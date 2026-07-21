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
#include <algorithm>                   // to convert SensType to upper case - input cleaning

#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sundials/sundials_math.h>    /* definition of ABS */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h>

#include <check_retval.h>
#include <jac_func.h>
#include <sundials_scope_guard.h>
// CRAN fix: replace SUNDIALS' default abort()-based error handler with one that
// records the error for the solver to raise via stop() (see the header)
#include <sundials_err_handler.h>


// #define Ith(v,i)    NV_Ith_S(v,i-1)         /* i-th vector component i=1..NEQ */

using namespace Rcpp;

// struct to use if R or Rcpp function is input as RHS function----------------
struct rhs_func_sens{
  Function rhs_eqn;
  NumericVector params;
  double rtol;
  NumericVector atol;
  SEXP jac_eqn;
  SEXP sens_eqn;              // user sensitivity RHS, or R_NilValue
  sundials_err_record *err;   // collects errors raised inside the callbacks
};

// function called by CVodeInit if user inputs R function
// Called by SUNDIALS from its own C code, so the body runs under
// sundials_callback_guard; see the note on rhs_function in rhs_func.cpp.
int rhs_function_sens(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data){

  // cast void pointer to pointer to struct and assign rhs to a Function
  struct rhs_func_sens *my_rhs_fun = (struct rhs_func_sens*)user_data;

  // nothing to record against, so just report the failure to SUNDIALS
  if(!my_rhs_fun){ return(-1); }

  return sundials_callback_guard(my_rhs_fun->err, [&]() -> int {

    // convert y (N_Vector) to NumericVector y1
    int y_len = NV_LENGTH_S(y);

    NumericVector y1(y_len);    // filled with zeros
    sunrealtype *y_ptr = N_VGetArrayPointer(y);
    for (int i = 0; i < y_len; i++){
      y1[i] = y_ptr[i];
    }

    NumericVector ydot1(y_len);    // filled with zeros

    Function rhs_fun = (*my_rhs_fun).rhs_eqn;                // function
    NumericVector p_values = (*my_rhs_fun).params;           // rate parameters

    if (!rhs_fun || (TYPEOF(rhs_fun) != CLOSXP)){
      stop("The RHS function is not an R or Rcpp function, stopping!");
    }

    // use the function to calculate value of RHS ----
    ydot1 = rhs_fun(t, y1, p_values);

    // guards the element-by-element copy below against a short return
    if (ydot1.length() != y_len){
      stop("The RHS function must return a vector of the same length as the state vector: expected %d, got %d",
           y_len, ydot1.length());
    }

    // convert NumericVector ydot1 to N_Vector ydot
    sunrealtype *ydot_ptr = N_VGetArrayPointer(ydot);
    for (int i = 0; i<  y_len; i++){
      ydot_ptr[i] = ydot1[i];
    }

    // If everything went smoothly, return 0
    return(0);
  });
}

// EwtSet Function. Computes the error weights at the current solution.
int ewt(N_Vector y, N_Vector w, void *user_data)
{
  sunrealtype yy, ww, rtol;

  struct rhs_func_sens *my_rhs_fun = (struct rhs_func_sens*)user_data;
  if(!my_rhs_fun){ return(-1); }

  rtol = my_rhs_fun->rtol;

  // Read the tolerances through a pointer rather than copying the vector.
  // CVODES calls this on every internal step, far more often than the RHS,
  // so a copy here is on the hottest path in the solve.
  const double *atol = my_rhs_fun->atol.begin();

  sunindextype y_len = N_VGetLength_Serial(y);

  for (int i=1; i<=y_len; i++) {

    yy = NV_Ith_S(y,i-1);
    ww = rtol * SUNRabs(yy) + atol[i-1];
    if (ww <= 0.0) return (-1);
    NV_Ith_S(w,i-1) = 1.0/ww;

    // Rcpp::Rcout << NV_Ith_S(w,i-1) << "\n";
  }

  return(0);
}

// Call the jacobian function if manual jacobian is provided
static int jac_cvodes(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix JAC,
                      void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    struct rhs_func_sens *data = (struct rhs_func_sens*)user_data;
    if (!data) { return -1; }
    return sundials_callback_guard(data->err, [&]() -> int {
      return jac_eval(t, y, JAC, data->jac_eqn, data->params);
    });
}

// Sensitivity right-hand side, used when the caller supplies `sensitivity`.
// CVODES calls this once per parameter (a CVSensRhs1Fn); it delegates to the R
// function stored in data->sens_eqn, whose signature is
//   sens_rhs(t, y, ydot, iS, yS, p)  ->  numeric vector of length(y)
// returning d(yS_iS)/dt = J %*% yS_iS + df/dp_iS. iS is passed 1-based to R.
// Called by SUNDIALS from its own C code, so the body runs under
// sundials_callback_guard; see the note on rhs_function in rhs_func.cpp.
static int sens_rhs1_cvodes(int Ns, sunrealtype t, N_Vector y, N_Vector ydot,
                            int iS, N_Vector yS, N_Vector ySdot,
                            void *user_data, N_Vector tmp1, N_Vector tmp2) {
    struct rhs_func_sens *data = (struct rhs_func_sens*)user_data;
    if (!data) { return -1; }
    return sundials_callback_guard(data->err, [&]() -> int {

      int n = NV_LENGTH_S(y);

      NumericVector y1(n), ydot1(n), yS1(n);
      sunrealtype *y_ptr    = N_VGetArrayPointer(y);
      sunrealtype *ydot_ptr = N_VGetArrayPointer(ydot);
      sunrealtype *yS_ptr   = N_VGetArrayPointer(yS);
      for (int i = 0; i < n; i++) { y1[i] = y_ptr[i]; ydot1[i] = ydot_ptr[i]; yS1[i] = yS_ptr[i]; }

      Function sens_fun(data->sens_eqn);
      // iS arrives 0-based from CVODES; hand R the 1-based parameter index
      NumericVector ySdot1 = sens_fun(t, y1, ydot1, iS + 1, yS1, data->params);

      // guards the element-by-element copy below against a short return
      if (ySdot1.length() != n) {
        stop("The sensitivity function must return a vector of the same length as the state vector: expected %d, got %d",
             n, ySdot1.length());
      }

      sunrealtype *ySdot_ptr = N_VGetArrayPointer(ySdot);
      for (int i = 0; i < n; i++) ySdot_ptr[i] = ySdot1[i];

      return 0;
    });
}
//------------------------------------------------------------------------------
//' cvodes
//'
//' CVODES solver to solve ODEs and calculate sensitivities
//'@param time_vector time vector
//'@param IC Initial Conditions
//'@param input_function Right Hand Side function of ODEs
//'@param Parameters Parameters input to ODEs
//'@param reltolerance Relative Tolerance (a scalar, default value  = 1e-04)
//'@param abstolerance Absolute Tolerance (a scalar or vector with length equal to ydot, default = 1e-04)
//'@param SensType Sensitivity Type - allowed values are "STG" (for Staggered, default) or "SIM" (for Simultaneous)
//'@param ErrCon Error Control - allowed values are TRUE or FALSE (default)
//'@param jacobian (Optional) Jacobian of the RHS with signature \code{function(t, y, p)}. Default is NULL
//'@param sensitivity (Optional) Sensitivity right-hand side with signature \code{function(t, y, ydot, iS, yS, p)} returning the derivative \code{d(yS_iS)/dt = J \%*\% yS_iS + df/dp_iS} as a numeric vector of \code{length(y)}, where \code{iS} is the 1-based parameter index. Default is NULL, in which case the sensitivity equations are approximated by finite differences of the RHS
//'@returns A Matrix. First column is the time-vector, the next y * p columns are sensitivities of y1 w.r.t all parameters, then y2 w.r.t all parameters etc. y is the state vector, p is the parameter vector
//'@example /inst/examples/cvs_Roberts_dns.r
// [[Rcpp::export]]
NumericMatrix cvodes(NumericVector time_vector, NumericVector IC,
                      SEXP input_function,
                      NumericVector Parameters,
                      double reltolerance = 0.0001,
                      NumericVector abstolerance = 0.0001,
                      std::string SensType = "STG",
                      bool ErrCon = 'F',
                      Nullable<Function> jacobian = R_NilValue,
                      Nullable<Function> sensitivity = R_NilValue){

  int flag;

  int time_vec_len = time_vector.length();
  double time;
  int NOUT = time_vec_len;
  sunrealtype T0 = SUN_RCONST(time_vector[0]);     // Initial Time

  // Initial Conditions
  int y_len = IC.length();

  // Relative Tolerance
  sunrealtype reltol = reltolerance;

  // absolute tolerance is either length == 1 or equal to length of IC
  // If abstol is not equal to 1 or abstol is not equal to IC, then stop
  int abstol_len = abstolerance.length();
  if(abstol_len != 1 && abstol_len != y_len){
     stop("Absolute tolerance must be a scalar or a vector of same length as IC \n");
  }
  // Create an abstolerance vector equal to IC length
  Rcpp::NumericVector abstol(y_len);
  if(abstol_len == 1){
    for (int i = 0; i<y_len; i++){
      abstol[i] = abstolerance[0];
    }
  }
  if (abstol_len == y_len){
    for (int i = 0; i<y_len; i++){
      abstol[i] = abstolerance[i];
    }
  }

  // Number of sensitivity parameters; needed by the guard to size the yS free
  int NP = Parameters.length();

  // Receives SUNDIALS errors. Declared before the guard so that it is
  // destroyed after it - the SUNContext freed there holds a pointer to it.
  sundials_err_record sun_err;

  // SUNDIALS objects, released by the guard below on every exit path
  SUNContext sunctx      = NULL;
  void *cvode_mem        = NULL;
  N_Vector y0            = NULL;
  N_Vector *yS           = NULL;
  SUNMatrix SM           = NULL;
  SUNLinearSolver LS     = NULL;

  // Free in the same order the trailing free calls used to; runs on the
  // normal return and on the exceptions thrown by stop() below. yS must be
  // freed with the same count it was cloned with (NP).
  auto sundials_cleanup = make_scope_guard([&]{
    if (y0)        N_VDestroy(y0);
    if (yS)        N_VDestroyVectorArray(yS, NP);
    if (cvode_mem) CVodeFree(&cvode_mem);
    if (LS)        SUNLinSolFree(LS);
    if (SM)        SUNMatDestroy(SM);
    if (sunctx)    SUNContext_Free(&sunctx);
  });

  // Set context for sundials
  SUNContext_Create(SUN_COMM_NULL, &sunctx);
  // CRAN fix: redirect SUNDIALS fatal errors to R instead of calling abort()
  SUNContext_PushErrHandler(sunctx, sundials_r_err_handler, &sun_err);
  sundials_check(sun_err);   // context creation is not otherwise checked

  // Set the initial conditions-------------------------------------------------
  y0 = N_VNew_Serial(y_len, sunctx);
  sundials_check(sun_err);   // vector allocations are not otherwise checked
  sunrealtype *y0_ptr = N_VGetArrayPointer(y0);
  for (int i = 0; i<y_len; i++){
    y0_ptr[i] = IC[i];
  }

  // Check SensType input ---------------------------------------------------------
  // convert the SensType argument to upper case first
  std::transform(SensType.begin(), SensType.end(), SensType.begin(), ::toupper);
  if(SensType.compare("STG") != 0 && SensType.compare("SIM") != 0){
    stop("SensType argument can be STG (for Staggered) or SIM (for simulated) \nsee cvode help for example usage");
  }

  //-------------------------------------------------------------------------------

  // Check Error Control input ----------------------------------------------------
  if(ErrCon != TRUE && ErrCon != FALSE){
    stop("ErrCon can only be either TRUE or FALSE");
  }

  // //----------------------------------------------------------------------------
  /* Set sensitivity initial conditions */
  yS = N_VCloneVectorArray(NP, y0);   // vector of sensitivities (yS)

  // Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval(cvode_mem, "CVodeCreate")) { sundials_stop(sun_err, "CVodeCreate", "Stopping cvodes, cannot allocate memory for CVODES!"); }

  //-- assign user input to the struct based on SEXP type of input_function
  if (!input_function){
    stop("Something is wrong with input function, stopping!");
  }

  // Rcout << TYPEOF(input_function) << "\n";

  if (!input_function){ stop("There is no input function, stopping!"); }

  if(TYPEOF(input_function) != CLOSXP) { stop("Incorrect input function type - input function can be an R or Rcpp function"); }

  // realtype *params = Parameters.begin();
  // Initialize the struct for user data
  // Jacobian and sensitivity RHS have initial value of NULL
  SEXP jac_sexp = R_NilValue;
  if (jacobian.isNotNull()) jac_sexp = as<SEXP>(jacobian);
  SEXP sens_sexp = R_NilValue;
  if (sensitivity.isNotNull()) sens_sexp = as<SEXP>(sensitivity);
  struct rhs_func_sens my_rhs_function = {input_function,
                                          Parameters,
                                          reltol,
                                          abstol,
                                          jac_sexp,
                                          sens_sexp,
                                          &sun_err};

  // setting the user_data in rhs function
  flag = CVodeSetUserData(cvode_mem, (void*)&my_rhs_function);
  if (check_retval(flag, "CVodeSetUserData")) { sundials_stop(sun_err, "CVodeSetUserData", "Stopping cvodes, something went wrong in setting user data!"); }

  /* Allocate space for CVODES */
  flag = CVodeInit(cvode_mem, rhs_function_sens, T0, y0);
  if (check_retval(flag, "CVodeInit")) { sundials_stop(sun_err, "CVodeInit", "Stopping cvodes, something went wrong in allocating space for CVODES!"); }

  /* Use private function to compute error weights */
  flag = CVodeWFtolerances(cvode_mem, ewt);
  if (check_retval(flag, "CVodeWFtolerances")) { sundials_stop(sun_err, "CVodeWFtolerances", "Stopping cvodes, something went wrong in computing error weights!"); }

  /* Create dense SUNMatrix */
  sunindextype y_len_M = y_len;
  SM = SUNDenseMatrix(y_len_M, y_len_M, sunctx);
  if (check_retval(SM, "SUNDenseMatrix")) { sundials_stop(sun_err, "SUNDenseMatrix", "Stopping cvodes, something went wrong in setting SUNDenseMatrix!"); }

  /* Create dense SUNLinearSolver */
  LS = SUNLinSol_Dense(y0, SM, sunctx);
  if (check_retval(LS, "SUNLinSol_Dense")) { sundials_stop(sun_err, "SUNLinSol_Dense", "Stopping cvodes, something went wrong in setting Linear Solver!"); }

  /* Attach the matrix and linear solver */
  flag = CVodeSetLinearSolver(cvode_mem, LS, SM);
  if (check_retval(flag, "CVodeSetLinearSolver")) { sundials_stop(sun_err, "CVodeSetLinearSolver", "Stopping cvodes, something went wrong in attaching SUNDenseMatrix and Linear Solver!"); }

  /* If manual Jacobian is provided, use the jacobian */
  if (jacobian.isNotNull()){
    flag = CVodeSetJacFn(cvode_mem, jac_cvodes);
    if(check_retval(flag, "CVodeSetJacFn")) { sundials_stop(sun_err, "CVodeSetJacFn", "Stopping cvodes, something went wrong in setting the Jacobian function!"); }
  }

  if (check_retval(yS, "N_VCloneVectorArray")) { sundials_stop(sun_err, "N_VCloneVectorArray", "Stopping cvodes, something went wrong in setting Sensitivity Array!"); }
  for (int is=0;is<NP;is++) N_VConst(SUN_RCONST(0.0), yS[is]);

  /* Call CVodeSensInit1 to activate forward sensitivity computations
   and allocate internal memory for CVODES related to sensitivity
   calculations. Computes the right-hand sides of the sensitivity
   ODE, one at a time */
  //
  // The fourth argument is the sensitivity right-hand side. When the caller
  // supplies `sensitivity`, sens_rhs1_cvodes is passed and CVODES uses that R
  // function to evaluate the sensitivity equations analytically. Passing NULL
  // (the default) instead tells CVODES to approximate it by finite differences
  // of the state right-hand side, which needs no extra function from the caller
  // but evaluates the state RHS an extra NP times per internal step - each an R
  // callback, the most expensive thing in the solve - and is only as accurate
  // as the difference quotients.
  //
  // (int) flag to set sensitivity solution method - see CVodeSensInit1
  int ism = CV_STAGGERED;
  if (SensType.compare("SIM") == 0) ism = CV_SIMULTANEOUS;
  CVSensRhs1Fn fS1 = sensitivity.isNotNull() ? sens_rhs1_cvodes : NULL;
  flag = CVodeSensInit1(cvode_mem, NP, ism, fS1, yS);
  if(check_retval(flag, "CVodeSensInit1")) { sundials_stop(sun_err, "CVodeSensInit1", "Stopping cvodes, something went wrong in calculating Sensitivities!"); }

  /* Call CVodeSensEEtolerances to estimate tolerances for sensitivity
   variables based on the rolerances supplied for states variables and
   the scaling factor pbar */
  flag = CVodeSensEEtolerances(cvode_mem);
  if(check_retval(flag, "CVodeSensEEtolerances")) { sundials_stop(sun_err, "CVodeSensEEtolerances", "Stopping cvodes, something went wrong in estimating tolerances for sensitivities!"); }

  /* Call CVodeSetSensParams to specify problem parameter information for
   sensitivity calculations */
  // struct rhs_func_sens *ptr = &my_rhs_function;        // struct UserData storing my data
  // p = (my_rhs_function.params).begin();
  // Rcout << (my_rhs_function.params).begin() << "\n";
  // Rcout << &(ptr->params[0]);
  flag = CVodeSetSensParams(cvode_mem, (my_rhs_function.params).begin(), Parameters.begin(), NULL);  // double *y = x.begin()
  if (check_retval(flag, "CVodeSetSensParams")) { sundials_stop(sun_err, "CVodeSetSensParams", "Stopping cvodes, something went wrong in setting Sensitivity Parameters!"); }

  // First row for initial conditions, First column is for time
  int y_len_1 = y_len + 1;
  NumericMatrix soln(Dimension(time_vec_len,y_len_1));

  // fill the first row of soln matrix with Initial Conditions
  soln(0,0) = time_vector[0];   // get the first time value
  for(int i = 0; i<y_len; i++){
    soln(0,i+1) = IC[i];
  }

  NumericMatrix sens(Dimension(time_vec_len, (y_len * NP) + 1));
  // Store the Sensitivity Results
  // Sensitivity of each entitiy w.r.t. parameters is zero at initial time
  sens(0,0) = time_vector[0];           // get the first time value
  for(int i = 0; i<y_len; i++){
    sens(0,i+1) = 0.0;
  }

  sunrealtype tout;  // For output times
  sunrealtype *yS_ptr = NULL;
  for(int iout = 0; iout < NOUT-1; iout++) {

    // output times start from the index after initial time
    tout = time_vector[iout+1];

    flag = CVode(cvode_mem, tout, y0, &time, CV_NORMAL);
    if (check_retval(flag, "CVode")) { sundials_stop(sun_err, "CVode", "Stopping cvodes, something went wrong in solving the system using CVODE!"); } // Something went wrong in solving it!
    if (flag == CV_SUCCESS) {

      // store results in soln matrix
      soln(iout+1, 0) = time;           // first column is for time
      for (int i = 0; i<y_len; i++){
        soln(iout+1, i+1) = y0_ptr[i];
      }
    }

    flag = CVodeGetSens(cvode_mem, &time, yS);
    if (check_retval(flag, "CVodeGetSens")) { sundials_stop(sun_err, "CVodeGetSens", "Stopping cvodes, something went wrong in calculating Sensitivities!"); }

    if (flag == CV_SUCCESS) {
      sens(iout+1, 0) = time;               // first column is for time
      for (int i = 0; i < NP; i++){
        yS_ptr = N_VGetArrayPointer(yS[i]);   // sensitivities w.r.t. param i
        for(int j = 0; j < y_len; j++){
          // store results in soln matrix
          sens(iout+1, y_len*i+j+1) = yS_ptr[j]; // sensitivities of y1 ... yJ w.r.t param i
        }
      }
    }

  }

  /* SUNDIALS objects are released by sundials_cleanup on scope exit */

  return sens;



}

