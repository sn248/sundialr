//   Copyright (c) 2016-2026, Satyaprakash Nayak
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions are
//   met:
//
//   Redistributions of source code must retain the above copyright
//   notice, this list of conditions and the following disclaimer.
//
//   Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in
//   the documentation and/or other materials provided with the
//   distribution.
//
//   Neither sundialr nor the names of its
//   contributors may be used to endorse or promote products derived
//   from this software without specific prior written permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <Rcpp.h>

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

#include <check_retval.h>
#include <rhs_func.h>
#include <jac_func.h>
#include <sundials_scope_guard.h>

// CRAN fix: replace SUNDIALS' default abort()-based error handler with one that
// records the error for the solver to raise via stop() (see the header)
#include <sundials_err_handler.h>

#define Ith(v,i)    NV_Ith_S(v,i-1)  /* i-th vector component i=1..NEQ */

using namespace Rcpp;


//------------------------------------------------------------
// call the jac_eval from header file included
// optional arguments tmp1, tmp2, tmp3 need to be included
static int jac_cvode(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix JAC,
                     void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

  struct rhs_func *data = (struct rhs_func*)user_data;
  if (!data) { return -1; }
  return sundials_callback_guard(data->err, [&]() -> int {
    return jac_eval(t, y, JAC, data->jac_eqn, data->params);
  });
}


//'cvode
//'
//' CVODE solver to solve stiff ODEs
//'@param time_vector time vector
//'@param IC Initial Conditions
//'@param input_function Right Hand Side function of ODEs
//'@param Parameters Parameters input to ODEs
//'@param reltolerance Relative Tolerance (a scalar, default value  = 1e-04)
//'@param abstolerance Absolute Tolerance (a scalar or vector with length equal to ydot (dy/dx), default = 1e-04)
//'@param jacobian (Optional) Jacobian of the RHS with signature \code{function(t, y, p)} returning an n-by-n matrix where entry [i,j] is d(ydot_i)/d(y_j). Default is NULL and SUNDIALS uses internal finite-difference approximation.
//'@returns A Matrix. First column is the time-vector, the other columns are values of y in order they are provided.
//'@example /inst/examples/cv_Roberts_dns.r
// [[Rcpp::export]]
NumericMatrix cvode(NumericVector time_vector, NumericVector IC,
                     SEXP input_function,
                     NumericVector Parameters,
                     double reltolerance = 0.0001,
                     NumericVector abstolerance = 0.0001,
                     Nullable<Function> jacobian = R_NilValue){

   int flag;

   int time_vec_len = time_vector.length();
   double time;
   int NOUT = time_vec_len;
   sunrealtype T0 = SUN_RCONST(time_vector[0]);  //RCONST(0.0); // Initial Time

   // Initial Conditions
   int y_len = IC.length();

   // Relative tolerance
   sunrealtype reltol = reltolerance;

   // Receives SUNDIALS errors. Declared before the guard so that it is
   // destroyed after it - the SUNContext freed there holds a pointer to it.
   sundials_err_record sun_err;

   // SUNDIALS objects, released by the guard below on every exit path
   SUNContext sunctx      = NULL;
   void *cvode_mem        = NULL;
   N_Vector y0            = NULL;
   N_Vector abstol        = NULL;
   SUNMatrix SM           = NULL;
   SUNLinearSolver LS     = NULL;

   // Free in the same order the trailing free calls used to; runs on the
   // normal return and on the exceptions thrown by stop() below.
   auto sundials_cleanup = make_scope_guard([&]{
     if (y0)        N_VDestroy(y0);
     if (abstol)    N_VDestroy(abstol);
     if (cvode_mem) CVodeFree(&cvode_mem);
     if (LS)        SUNLinSolFree(LS);
     if (SM)        SUNMatDestroy(SM);
     if (sunctx)    SUNContext_Free(&sunctx);
   });

   // Set Sundials context
   SUNContext_Create(SUN_COMM_NULL, &sunctx);
   // CRAN fix: redirect SUNDIALS fatal errors to R instead of calling abort()
   SUNContext_PushErrHandler(sunctx, sundials_r_err_handler, &sun_err);
   sundials_check(sun_err);   // context creation is not otherwise checked

   // Absolute tolerance
   // Set the vector absolute tolerance -----------------------------------------
   // abstol must be same length as IC
   int abstol_len = abstolerance.length();

   // absolute tolerance is either length == 1 or equal to length of IC
   // If abstol is not equal to 1 and abstol is not equal to IC, then stop
   if(abstol_len != 1 && abstol_len != y_len){
     stop("Absolute tolerance must be a scalar or a vector of same length as IC\n");
   }

   abstol = N_VNew_Serial(y_len, sunctx);
   sunrealtype *abstol_ptr = N_VGetArrayPointer(abstol);
   if(abstol_len == 1){
     // if a scalar is provided - use it to make a vector with same values
     for (int i = 0; i<y_len; i++){
       abstol_ptr[i] = abstolerance[0];
     }
   }
   if (abstol_len == y_len){
     for (int i = 0; i<y_len; i++){
       abstol_ptr[i] = abstolerance[i];
     }
   }

   // Set the initial conditions-------------------------------------------------
   y0 = N_VNew_Serial(y_len, sunctx);
   sundials_check(sun_err);   // vector allocations are not otherwise checked
   sunrealtype *y0_ptr = N_VGetArrayPointer(y0);
   for (int i = 0; i<y_len; i++){
     y0_ptr[i] = IC[i]; // NV_Ith_S(y0, i)
   }

   // Call CVodeCreate to create the solver memory and specify the
   // Backward Differentiation Formula (BDF)
   cvode_mem = CVodeCreate(CV_BDF, sunctx);
   if (check_retval(cvode_mem, "CVodeCreate")) {
     sundials_stop(sun_err, "CVodeCreate", "Something went wrong in assigning memory, stopping cvode!");
   }

   //-- assign user input to the struct based on SEXP type of input_function
   if (!input_function){ stop("There is no input function, stopping!"); }

   if(TYPEOF(input_function) != CLOSXP) { stop("Incorrect input function type - input function can be an R or Rcpp function"); }

   // Check if jacobian is not NULL, fill the jac_sexp with input value
   SEXP jac_sexp = R_NilValue;
   if (jacobian.isNotNull()) jac_sexp = as<SEXP>(jacobian);
   struct rhs_func my_rhs_function = {input_function, Parameters, jac_sexp, &sun_err};

   // setting the user_data in rhs function
   flag = CVodeSetUserData(cvode_mem, (void*)&my_rhs_function);
   if (check_retval(flag, "CVodeSetUserData")) { sundials_stop(sun_err, "CVodeSetUserData", "Stopping cvode, something went wrong in setting user data!"); }

   flag = CVodeInit(cvode_mem, rhs_function, T0, y0);
   if (check_retval(flag, "CVodeInit")) { sundials_stop(sun_err, "CVodeInit", "Stopping cvode, something went wrong in initializing CVODE!"); }

   // Call CVodeSVtolerances to specify the scalar relative tolerance and vector absolute tol
   flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
   if (check_retval(flag, "CVodeSVtolerances")) { sundials_stop(sun_err, "CVodeSVtolerances", "Stopping cvode, something went wrong in setting solver tolerances!"); }

   // Create dense SUNMatrix for use in linear solves
   sunindextype y_len_M = y_len;
   SM = SUNDenseMatrix(y_len_M, y_len_M, sunctx);
   if(check_retval(SM, "SUNDenseMatrix")) { sundials_stop(sun_err, "SUNDenseMatrix", "Stopping cvode, something went wrong in setting the dense matrix!"); }

   // Create dense SUNLinearSolver object for use by CVode
   LS = SUNLinSol_Dense(y0, SM, sunctx);
   if(check_retval(LS, "SUNLinSol_Dense")) { sundials_stop(sun_err, "SUNLinSol_Dense", "Stopping cvode, something went wrong in setting the linear solver!"); }

   // Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode
   flag = CVodeSetLinearSolver(cvode_mem, LS, SM);
   if(check_retval(flag, "CVodeSetLinearSolver")) { sundials_stop(sun_err, "CVodeSetLinearSolver", "Stopping cvode, something went wrong in setting the linear solver!"); }

   if (jacobian.isNotNull()) {
     flag = CVodeSetJacFn(cvode_mem, jac_cvode);
     if(check_retval(flag, "CVodeSetJacFn")) { sundials_stop(sun_err, "CVodeSetJacFn", "Stopping cvode, something went wrong in setting the Jacobian function!"); }
   }
   // NumericMatrix to store results - filled with 0.0

   // Call CVodeInit to initialize the integrator memory and specify the
   // user's right hand side function in y'=f(time,y),
   // // the inital time T0, and the initial dependent variable vector y.
   sunrealtype tout;  // For output times

   int y_len_1 = y_len + 1; // remove later
   NumericMatrix soln(Dimension(time_vec_len,y_len_1));  // remove later

   // fill the first row of soln matrix with Initial Conditions
   soln(0,0) = time_vector[0];   // get the first time value
   for(int i = 0; i<y_len; i++){
     soln(0,i+1) = IC[i];
   }

   for(int iout = 0; iout < NOUT-1; iout++) {

     // output times start from the index after initial time
     tout = time_vector[iout+1];

     flag = CVode(cvode_mem, tout, y0, &time, CV_NORMAL);

     // If something went wrong in solving it!
     if (check_retval(flag, "CVode")) {
       sundials_stop(sun_err, "CVode", "Stopping CVODE, something went wrong in solving the system of ODEs!");
     }

     if (flag == CV_SUCCESS) {
       // store results in soln matrix
       soln(iout+1, 0) = time;           // first column is for time
       for (int i = 0; i<y_len; i++){
         soln(iout+1, i+1) = y0_ptr[i];
       }
     }
   }

   // SUNDIALS objects are released by sundials_cleanup on scope exit

   return soln;

}
 //--- cvode definition ends ----------------------------------------------------
