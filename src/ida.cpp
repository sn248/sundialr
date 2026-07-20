//   Copyright (c) 2024, Satyaprakash Nayak
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

#include <ida/ida.h>                          /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>           /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h>        /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h>        /* access to dense SUNLinearSolver      */
#include <sunnonlinsol/sunnonlinsol_newton.h> /* access to Newton SUNNonlinearSolver  */
#include <sundials/sundials_types.h>          /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>           /* defs. of SUNRabs, SUNRexp, etc.      */

#include <check_retval.h>
#include <jac_func.h>
#include <sundials_scope_guard.h>
// CRAN fix: replace SUNDIALS' default abort()-based error handler with one that
// records the error for the solver to raise via stop() (see the header)
#include <sundials_err_handler.h>

using namespace Rcpp;

// struct to use if R or Rcpp function is input as RHS residual function
struct res_func{
  Function res_eqn;
  NumericVector params;
  SEXP jac_eqn;
  sundials_err_record *err;   // collects errors raised inside the callbacks
};

// function called by IDAInit if user inputs R function
// Called by SUNDIALS from its own C code, so the body runs under
// sundials_callback_guard; see the note on rhs_function in rhs_func.cpp.
int res_function(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data){

  // cast void pointer to pointer to struct and assign rhs residual to a Function
  struct res_func *my_res_fun = (struct res_func*)user_data;

  // nothing to record against, so just report the failure to SUNDIALS
  if(!my_res_fun){ return(-1); }

  return sundials_callback_guard(my_res_fun->err, [&]() -> int {

    // convert y to NumericVector y1
    int yy_len = NV_LENGTH_S(yy);

    // NumericVector analog of yy
    NumericVector yy1(yy_len);                      // filled with zeros
    sunrealtype *yy_ptr = N_VGetArrayPointer(yy);
    for (int i = 0; i < yy_len; i++){
      yy1[i] = yy_ptr[i];                           // Ith(y,i+1);
    }

    int yp_len = NV_LENGTH_S(yp);
    // NumericVector analog of yp
    NumericVector yp1(yp_len);                      // filled with zeros
    sunrealtype *yp_ptr = N_VGetArrayPointer(yp);
    for (int i = 0; i < yp_len; i++){
      yp1[i] = yp_ptr[i];                           // Ith(y,i+1);
    }

    // NumericVector analog of rr
    int rr_len = NV_LENGTH_S(rr);
    NumericVector rr1(rr_len);                      // filled with zeros

    Function res_fun = (*my_res_fun).res_eqn;           // function
    NumericVector p_values = (*my_res_fun).params;      // rate parameters

    if (!res_fun || (TYPEOF(res_fun) != CLOSXP)){
      stop("The residual function is not an R or Rcpp function, stopping!");
    }

    // use the function to calculate value of RHS ----
    rr1 = res_fun(t, yy1, yp1, p_values);

    // guards the element-by-element copy below against a short return
    if (rr1.length() != rr_len){
      stop("The residual function must return a vector of the same length as the state vector: expected %d, got %d",
           rr_len, rr1.length());
    }

    // convert NumericVector rr1 to N_Vector rr
    sunrealtype *rr_ptr = N_VGetArrayPointer(rr);
    for (int i = 0; i <  rr_len; i++){
      rr_ptr[i] = rr1[i];
    }

    // everything went smoothly
    return(0);
  });
}
//---RHS residual function definition ends -------------------------------------

//-- Manual Jacobian -----------------------------------------------------------
/*IDA's Jacobian callback has a different SUNDIALS signature than CVODE's:
 * it receives cj (a scalar IDA manages internally) and yp (current y').
 * The matrix to fill is dF/dy + cj * dF/dy' - the user must return this
 * combined matrix from their R function
 * R function signature for IDA Jacobian:
 * f(t, y, ydot, cj, p) -> n-by-n matrix of dF/dy + cj * dF/dy'
 * Note the IDA callback type is IDALsJacFn,
 * which has cj and yp as extra arguments vs. CVLsJacFn:
*/

static int jac_ida(sunrealtype t, sunrealtype cj,
                     N_Vector yy, N_Vector yp, N_Vector rr, SUNMatrix JAC,
                     void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    struct res_func *data = (struct res_func*)user_data;
    if (!data) { return -1; }
    return sundials_callback_guard(data->err, [&]() -> int {
      return jac_eval_ida(t, cj, yy, yp, JAC, data->jac_eqn, data->params);
    });
}
//------------------------------------------------------------------------------

//'ida
//'
//' IDA solver to solve stiff DAEs
//'@param time_vector time vector
//'@param IC Initial Value of y
//'@param IRes Inital Value of ydot
//'@param input_function Right Hand Side function of DAEs
//'@param Parameters Parameters input to ODEs
//'@param reltolerance Relative Tolerance (a scalar, default value  = 1e-04)
//'@param abstolerance Absolute Tolerance (a scalar or vector with length equal to ydot, default = 1e-04)
//'@param jacobian (Optional) Jacobian with signature \code{function(t, y, ydot, cj, p)} returning an n-by-n matrix of \code{dF/dy + cj*dF/dydot}. Default NULL.
//'@returns A Matrix. First column is the time-vector, the other columns are values of y in order they are provided.
//'@example /inst/examples/ida_Roberts_dns.r
// [[Rcpp::export]]
NumericMatrix ida(NumericVector time_vector, NumericVector IC,
                  NumericVector IRes, SEXP input_function,
                  NumericVector Parameters,
                  double reltolerance = 0.0001,
                  NumericVector abstolerance = 0.0001,
                  Nullable<Function> jacobian = R_NilValue){

  int time_vec_len = time_vector.length();
  int y_len = IC.length();
  // Receives SUNDIALS errors. Declared before the guard so that it is
  // destroyed after it - the SUNContext freed there holds a pointer to it.
  sundials_err_record sun_err;

  // SUNDIALS objects, released by the guard below on every exit path
  SUNContext sunctx        = NULL;
  void *ida_mem            = NULL;
  N_Vector yy0             = NULL;
  N_Vector yp0             = NULL;
  N_Vector abstol          = NULL;
  SUNMatrix SM             = NULL;
  SUNLinearSolver LS       = NULL;
  SUNNonlinearSolver NLS   = NULL;

  // Free in the same order the trailing free calls used to; runs on the
  // normal return and on the exceptions thrown by stop() below - including
  // the two input checks immediately after the context is created.
  auto sundials_cleanup = make_scope_guard([&]{
    if (ida_mem) IDAFree(&ida_mem);
    if (NLS)     SUNNonlinSolFree(NLS);
    if (LS)      SUNLinSolFree(LS);
    if (SM)      SUNMatDestroy(SM);
    if (abstol)  N_VDestroy(abstol);
    if (yy0)     N_VDestroy(yy0);
    if (yp0)     N_VDestroy(yp0);
    if (sunctx)  SUNContext_Free(&sunctx);
  });

  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  if(y_len != IRes.length()){ stop("IC (Initial Conditions) and IRes (Residuals) should be of same length"); }

  int abstol_len = abstolerance.length();
  // absolute tolerance is either length == 1 or equal to length of IC
  // If abstol is not equal to 1 and abstol is not equal to IC, then stop
  if(abstol_len != 1 && abstol_len != y_len){
    stop("Absolute tolerance must be a scalar or a vector of same length as IC \n");
  }

  // CRAN fix: redirect SUNDIALS fatal errors to R instead of calling abort()
  SUNContext_PushErrHandler(sunctx, sundials_r_err_handler, &sun_err);
  sundials_check(sun_err);   // context creation is not otherwise checked

  int flag;
  sunrealtype reltol = reltolerance;
  sunrealtype T0 = SUN_RCONST(time_vector[0]);     //RCONST(0.0);  // Initial Time

  double time;
  int NOUT = time_vec_len;

  // Set the vector absolute tolerance -----------------------------------------
  // abstol must be same length as IC
  // NOTE: sized by y_len (not abstol_len) because the scalar-abstolerance
  // branch below writes y_len entries into abstol_ptr regardless of
  // abstol_len; sizing by abstol_len (== 1 in that branch) overflows the heap.
  abstol = N_VNew_Serial(y_len, sunctx);
  sunrealtype *abstol_ptr = N_VGetArrayPointer(abstol);
  if(abstol_len == 1){
    // if a scalar is provided - use it to make a vector with same values
    for (int i = 0; i<y_len; i++){
      abstol_ptr[i] = abstolerance[0];
    }
  }
  else if (abstol_len == y_len){

    for (int i = 0; i<abstol_len; i++){
      abstol_ptr[i] = abstolerance[i];
    }
  }

  //----------------------------------------------------------------------------
  // // Set the initial values of y -----------------------------------------------
  yy0 = N_VNew_Serial(y_len, sunctx);          // declared as yy0 to be consistent with example C code
  sunrealtype *yy0_ptr = N_VGetArrayPointer(yy0);
  for (int i = 0; i<y_len; i++){
    yy0_ptr[i] = IC[i];
  }

  // // Set the initial values of ydot --------------------------------------------
  yp0 = N_VNew_Serial(y_len, sunctx);         // declared as yp0 to be consistent with example C code
  sundials_check(sun_err);   // vector allocations are not otherwise checked
  sunrealtype *yp0_ptr = N_VGetArrayPointer(yp0);
  for (int i = 0; i<y_len; i++){
    yp0_ptr[i] = IRes[i];
  }

  // //----------------------------------------------------------------------------
  /* Call IDACreate and IDAInit to initialize IDA memory */
  ida_mem = IDACreate(sunctx);

  if(check_retval((void *)ida_mem, "IDACreate", 0)) { sundials_stop(sun_err, "IDACreate", "Stopping IDA, something went wrong in allocating memory!"); }

  // -- assign user input to the struct based on SEXP type of input_function
  if(!input_function){ stop("Something is wrong with the input function, stopping!"); }

  if(TYPEOF(input_function) != CLOSXP) { stop("Incorrect input function type - input function can be an R or Rcpp function"); }

  SEXP jac_sexp = R_NilValue;  // for manual jacobian, if provided
  if (jacobian.isNotNull()) jac_sexp = as<SEXP>(jacobian);
  struct res_func my_res_function = {input_function, Parameters, jac_sexp, &sun_err};

  // setting the user data in the rhs residual function
  flag = IDASetUserData(ida_mem, (void*)&my_res_function);
  if (check_retval(&flag, "IDASetUserData", 1)) { sundials_stop(sun_err, "IDASetUserData", "Stopping IDA, something went wrong in setting user data!"); }

  flag = IDAInit(ida_mem, res_function, T0, yy0, yp0);
  if(check_retval(&flag, "IDAInit", 1)) { sundials_stop(sun_err, "IDAInit", "Stopping, something went wrong in initializing IDA!"); };

  /* Call IDASVtolerances to set tolerances */
  flag = IDASVtolerances(ida_mem, reltol, abstol);
  if(check_retval(&flag, "IDASVtolerances", 1)) { sundials_stop(sun_err, "IDASVtolerances", "Stopping, something went wrong in setting tolerances!"); };

  /* Create dense SUNMatrix for use in linear solves */
  sunindextype y_len_M = y_len;
  SM = SUNDenseMatrix(y_len_M, y_len_M, sunctx);
  if(check_retval((void *)SM, "SUNDenseMatrix", 0)) { sundials_stop(sun_err, "SUNDenseMatrix", "Stopping IDA, something went wrong in setting the dense matrix!"); }


  // Create dense SUNLinearSolver object for use by IDA
  LS = SUNLinSol_Dense(yy0, SM, sunctx);
  if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) { sundials_stop(sun_err, "SUNLinSol_Dense", "Stopping IDA, something went wrong in setting the linear solver!"); }

  /* Attach the matrix and linear solver */
  flag = IDASetLinearSolver(ida_mem, LS, SM);
  if(check_retval(&flag, "IDASetLinearSolver", 1))  { sundials_stop(sun_err, "IDASetLinearSolver", "Stopping IDA, something went wrong in setting the linear solver!"); }

  // Add user-provided Jacobian, if not NULL
  if (jacobian.isNotNull()) {
    flag = IDASetJacFn(ida_mem, jac_ida);
    if(check_retval(&flag, "IDASetJacFn", 1)) { sundials_stop(sun_err, "IDASetJacFn", "Stopping IDA, something went wrong in setting the Jacobian function!"); }
  }

  /* Create Newton SUNNonlinearSolver object. IDA uses a
   * Newton SUNNonlinearSolver by default, so it is unecessary
   * to create it and attach it. */
  NLS = SUNNonlinSol_Newton(yy0, sunctx);
  if(check_retval((void *)NLS, "SUNNonlinSol_Newton", 0))  { sundials_stop(sun_err, "SUNNonlinSol_Newton", "Stopping IDA, something went wrong in creating the Non-linear Solver in IDA!"); }

  /* Attach the nonlinear solver */
  flag = IDASetNonlinearSolver(ida_mem, NLS);
  if(check_retval(&flag, "IDASetNonlinearSolver", 1)) { sundials_stop(sun_err, "IDASetNonlinearSolver", "Stopping IDA, something went wrong in attaching the Non-linear Solver in IDA!"); }

  /* In loop, call IDASolve, print results, and test for error.
   Break out of loop when NOUT preset output times have been reached. */

  sunrealtype tout;  // For output times

  int y_len_1 = y_len + 1; // remove later
  NumericMatrix soln(Dimension(time_vec_len,y_len_1));  // remove later


  // nothing to do with Events - single initialization of the ODE system
  // First row for initial conditions, First column is for time
  // int y_len_1 = y_len + 1;
  // NumericMatrix soln(Dimension(time_vec_len,y_len_1));

  // fill the first row of soln matrix with Initial Conditions
  soln(0,0) = time_vector[0];   // get the first time value
  for(int i = 0; i<y_len; i++){
    soln(0,i+1) = IC[i];
  }

  for(int iout = 0; iout < NOUT-1; iout++) {

    // output times start from the index after initial time
    tout = time_vector[iout+1];

    flag = IDASolve(ida_mem, tout, &time, yy0, yp0, IDA_NORMAL);

    if (check_retval(&flag, "IDASolve", 1)) {
      sundials_stop(sun_err, "IDASolve", "Stopping IDA, something went wrong in solving the system of DAEs!");
    } // Something went wrong in solving it!

    if (flag == IDA_SUCCESS) {

      // store results in soln matrix
      soln(iout+1, 0) = time;           // first column is for time
      for (int i = 0; i<y_len; i++){
        soln(iout+1, i+1) = yy0_ptr[i];
      }
    }
  }

  /* SUNDIALS objects are released by sundials_cleanup on scope exit */

  return soln;

}

