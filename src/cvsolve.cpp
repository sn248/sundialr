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

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>    /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h>   /* definition of type realtype */
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

#include <check_retval.h>
#include <rhs_func.h>
#include <jac_func.h>
#include <sundials_scope_guard.h>

// CRAN fix: replace SUNDIALS' default abort()-based error handler with one that
// records the error for the solver to raise via stop() (see the header)
#include <sundials_err_handler.h>
#include "sortTimes.h"

using namespace Rcpp;
using namespace arma;

//------------------------------------------------------------
// call the jac_eval from header file included
// optional arguments tmp1, tmp2, tmp3 need to be included
static int jac_cvsolve(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix JAC,
                       void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

  struct rhs_func *data = (struct rhs_func*)user_data;
  if (!data) { return -1; }
  return sundials_callback_guard(data->err, [&]() -> int {
    return jac_eval(t, y, JAC, data->jac_eqn, data->params);
  });
}


//------------------------------------------------------------------------------
//'cvsolve
//'
//'CVSOLVE solver to solve stiff ODEs with discontinuties
//'@param time_vector time vector
//'@param IC Initial Conditions
//'@param input_function Right Hand Side function of ODEs
//'@param Parameters Parameters input to ODEs
//'@param Events Discontinuities in the solution (a DataFrame, default value is NULL). Three columns, names ignored: the 1-based index of the state, the time of the discontinuity, and the value to add to that state at that time. The value is always added to the current value of the state, including at the initial time, so the initial conditions in \code{IC} are the starting point and an event at t = 0 adds to them.
//'@param reltolerance Relative Tolerance (a scalar, default value  = 1e-04)
//'@param abstolerance Absolute Tolerance (a scalar or vector with length equal to ydot, default = 1e-04)
//'@param jacobian (Optional) Jacobian of the RHS with signature \code{function(t, y, p)} returning an n-by-n matrix where entry [i,j] is d(ydot_i)/d(y_j). Default is NULL and SUNDIALS uses internal finite-difference approximation.
//'@returns A Matrix. First column is the time-vector, the other columns are values of y in order they are provided.
//'@example /inst/examples/cvsolve_1D.r
// [[Rcpp::export]]
NumericMatrix cvsolve(NumericVector time_vector, NumericVector IC,
                      SEXP input_function,
                      NumericVector Parameters,
                      Nullable<DataFrame> Events = R_NilValue,
                      double reltolerance = 0.0001,
                      NumericVector abstolerance = 0.0001,
                      Nullable<Function> jacobian = R_NilValue){

  int y_len = IC.length();
  int NSTATES = IC.length();
  int abstol_len = abstolerance.length();

  // absolute tolerance is either length == 1 or equal to length of IC
  // If abstol is not equal to 1 and abstol is not equal to IC, then stop
  if(abstol_len != 1 && abstol_len != y_len){
    stop("Absolute tolerance must be a scalar or a vector of same length as IC \n");
  }

  // Receives SUNDIALS errors. Declared before the guard so that it is
  // destroyed after it - the SUNContext freed there holds a pointer to it.
  sundials_err_record sun_err;

  // SUNDIALS objects, released by the guard below on every exit path
  SUNContext sunctx      = NULL;
  void *cvode_mem        = NULL;
  N_Vector y0            = NULL;
  N_Vector abstol        = NULL;
  N_Vector constraints   = NULL;
  SUNMatrix SM           = NULL;
  SUNLinearSolver LS     = NULL;

  // Free in the same order the trailing free calls used to; runs on the
  // normal return and on the exceptions thrown by stop() below. constraints
  // was previously never freed at all, on any path.
  auto sundials_cleanup = make_scope_guard([&]{
    if (y0)          N_VDestroy(y0);
    if (abstol)      N_VDestroy(abstol);
    if (constraints) N_VDestroy(constraints);
    if (cvode_mem)   CVodeFree(&cvode_mem);
    if (LS)          SUNLinSolFree(LS);
    if (SM)          SUNMatDestroy(SM);
    if (sunctx)      SUNContext_Free(&sunctx);
  });

  SUNContext_Create(SUN_COMM_NULL, &sunctx);
  // CRAN fix: redirect SUNDIALS fatal errors to R instead of calling abort()
  SUNContext_PushErrHandler(sunctx, sundials_r_err_handler, &sun_err);
  sundials_check(sun_err);   // context creation is not otherwise checked

  int flag;
  sunrealtype reltol = reltolerance;

  sunrealtype T0 = SUN_RCONST(time_vector[0]);     //RCONST(0.0);  // Initial Time

  double time;

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
  if (abstol_len == y_len){
    for (int i = 0; i<abstol_len; i++){
      abstol_ptr[i] = abstolerance[i];
    }
  }

  SEXP jac_sexp = R_NilValue;
  if (jacobian.isNotNull()) jac_sexp = as<SEXP>(jacobian);


  // If Events is not NULL, change IC and generate a combined dataset-----------
  // Combine the time vector and Events vector into a single dataframe
  // If Events is null, then TCOMB is same as time_vec (i.e, is a vector)
  // If Events is not null, then TCOMB is a matrix with 4 columns
  // Column 1 - cpp index of state with discontinuity, -1 for all sampling time points
  // Column 2 - Time of discontinuity or sampling
  // Column 3 - Value at the time of discontinuity, 0 for sampling
  // Column 4 - 0 for sampling and 1 for discontinuity time point
  // An event at t = 0 adds to the IC, the same way events at later times add
  // to the current state
  NumericMatrix TCOMB(time_vector.length(), 4); // a matrix with 1 column
  NumericMatrix::Column col1 = TCOMB(_, 1);     /// reference to the second column (containing time points)
  col1 = time_vector;                           // set to be equal to initial time vector
  NumericMatrix::Column col0 = TCOMB(_,0);      // set state index to -1
  NumericVector stateSamplingTimes(TCOMB.nrow(), -1.0);
  col0 = stateSamplingTimes;

  NumericVector ICchanged(NSTATES);              // ICchanged created to store updated IC - if changed by Events
  for (int i = 0; i < NSTATES; i++){
    ICchanged[i] = IC[i];
  }

  if(Events.isNotNull()){

    Rcpp::DataFrame Events_DF(Events);

    // Check the state index in the first column.
    // It is converted to a 0-based index below and then used to subscript the
    // state vector directly (ICchanged(...) and y0_ptr[...]), so a value that
    // is out of range, fractional or NA is not merely wrong - it reads and
    // writes outside the state vector. Reject anything that is not a whole
    // number in 1..y_len.
    Rcpp::NumericVector EventsDF_state_col = Events_DF[0];
    for(int i = 0; i < EventsDF_state_col.length(); i++){
      double state_i = EventsDF_state_col[i];
      if(ISNAN(state_i) || state_i < 1 || state_i > y_len ||
         state_i != std::floor(state_i)){
        stop("The state index in the first column of the Events dataframe must be a whole number between 1 and the number of states");
      }
    }

    // Check the event times in the second column lie inside the output window.
    // An event before the first output time cannot be applied - the solve
    // starts at time_vector[0] - and previously produced a row of zeros while
    // the event itself was silently dropped.
    Rcpp::NumericVector EventsDF_time_col = Events_DF[1];
    for(int i = 0; i < EventsDF_time_col.length(); i++){
      double time_i = EventsDF_time_col[i];
      if(ISNAN(time_i) || time_i < time_vector[0] || time_i > Rcpp::max(time_vector)){
        stop("The event time in the second column of the Events dataframe must lie between the first and last time points");
      }
    }

    // Sort the Event DataFrame to combine discontinuities and sampling data points
    TCOMB = sorted_times(Events_DF, time_vector, NSTATES);  // get the sorted  Combined time matrix
    // Decrease State index by 1 to convert R state index to internal cpp index, all sampling times get state index of -1
    for(int i = 0; i < TCOMB.nrow(); i++){
      TCOMB(i,0) = TCOMB(i,0) - 1;
    }
    // Apply the events at t = 0 to the initial conditions. These ADD to the
    // IC, exactly as events at every later time add to the current state, so
    // the IC the user supplied is always the starting point and several events
    // at t = 0 on one state accumulate. Before 0.1.8 an event at t = 0
    // replaced the IC instead, which silently discarded it and made t = 0
    // behave differently from every other time.
    for(int i = 0; i < TCOMB.nrow();  i++){
      if(TCOMB(i,1) == time_vector[0] && TCOMB(i,3) == 1){
        int disc_index = static_cast<int>(TCOMB(i,0));
        ICchanged(disc_index) = ICchanged(disc_index) + TCOMB(i,2);
      }
    }
  }

  int NOUT = TCOMB.nrow();

  // Set the initial conditions-------------------------------------------------
  y0 = N_VNew_Serial(y_len, sunctx);
  sunrealtype *y0_ptr = N_VGetArrayPointer(y0);
  for (int i = 0; i<y_len; i++){
    y0_ptr[i] = ICchanged[i];
  }

  // Set constraints to all 1's for non-negative solution values.----------------
  constraints = N_VNew_Serial(y_len, sunctx);
  sundials_check(sun_err);   // vector allocations are not otherwise checked
  N_VConst(1, constraints);

  // Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void *) cvode_mem, "CVodeCreate", 0)) { sundials_stop(sun_err, "CVodeCreate", "Stopping cvsolve!"); }

  //-- assign user input to the struct based on SEXP type of input_function
  if (!input_function){ stop("Something is wrong with input function, stopping!"); }

  //-- assign user input to the struct based on SEXP type of input_function
  if (!input_function){ stop("There is no input function, stopping!"); }

  // order of input is rhs input function, Parameters and User-supplied Jacobian (optional)
  struct rhs_func my_rhs_function = {input_function, Parameters, jac_sexp, &sun_err};

  // setting the user_data in rhs function
  flag = CVodeSetUserData(cvode_mem, (void*)&my_rhs_function);
  if (check_retval(&flag, "CVodeSetUserData", 1)) { sundials_stop(sun_err, "CVodeSetUserData", "Stopping cvsolve, something went wrong in setting user data!"); }

  flag = CVodeInit(cvode_mem, rhs_function, T0, y0);
  if (check_retval(&flag, "CVodeInit", 1)) { sundials_stop(sun_err, "CVodeInit", "Stopping cvsolve, something went wrong in initializing CVODE!"); }

  // Call CVodeSVtolerances to specify the scalar relative tolerance and vector absolute tol
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&flag, "CVodeSVtolerances", 1)) { sundials_stop(sun_err, "CVodeSVtolerances", "Stopping cvsolve, something went wrong in setting solver tolerances!"); }

  // Create dense SUNMatrix for use in linear solves
  sunindextype y_len_M = y_len;
  SM = SUNDenseMatrix(y_len_M, y_len_M, sunctx);
  if(check_retval((void *)SM, "SUNDenseMatrix", 0)) { sundials_stop(sun_err, "SUNDenseMatrix", "Stopping cvsolve, something went wrong in setting the dense matrix!"); }

  // Create dense SUNLinearSolver object for use by CVode
  LS = SUNLinSol_Dense(y0, SM, sunctx);
  if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) { sundials_stop(sun_err, "SUNLinSol_Dense", "Stopping cvsolve, something went wrong in setting the linear solver!"); }

  // Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode
  flag = CVodeSetLinearSolver(cvode_mem, LS, SM);
  if(check_retval(&flag, "CVodeSetLinearSolver", 1)) { sundials_stop(sun_err, "CVodeSetLinearSolver", "Stopping cvsolve, something went wrong in setting the linear solver!"); }

  if (jacobian.isNotNull()) {
    flag = CVodeSetJacFn(cvode_mem, jac_cvsolve);
    if(check_retval(&flag, "CVodeSetJacFn", 1)) { sundials_stop(sun_err, "CVodeSetJacFn", "Stopping cvsolve, something went wrong in setting the Jacobian function!"); }
  }

  // Call CVodeSetConstraints to initialize constraints
  flag = CVodeSetConstraints(cvode_mem, constraints);
  if(check_retval(&flag, "CVodeSetConstraints", 1)) { sundials_stop(sun_err, "CVodeSetConstraints", "Stopping cvsolve, something went wrong in setting constraints!"); }

  sunrealtype tout;  // For output times

  // Solution vector has length equal to number of rows in TCOMB
  // Solution vector has width equal to number of IC + 1 (first column for time)
  int soln_rows = TCOMB.nrow();
  NumericMatrix soln(Dimension(soln_rows,y_len + 1));

  // fill the first row of soln matrix with Initial Conditions
  soln(0,0) = TCOMB(0, 1);   // get the first time value
  for(int i = 0; i<y_len; i++){
    soln(0,i+1) = ICchanged[i];
  }

  for(int iout = 0; iout < NOUT-1; iout++) {

    sunrealtype tprev = TCOMB(iout, 1);
    // output times start from the index after initial time
    tout = TCOMB(iout + 1, 1);

    // A record at the initial time is already accounted for by the initial
    // conditions, and one at the same time as the previous record was already
    // applied there. Neither advances the solution, so carry the current state
    // into the row rather than leaving it zero-filled.
    if (tout == time_vector[0] || tout == tprev){
      soln(iout+1, 0) = tout;
      for (int i = 0; i<y_len; i++){
        soln(iout+1, i+1) = y0_ptr[i];
      }
      continue;
    }
    else {

      // integrate upto the next time point (whether a sampling time or a discontinuity)
      flag = CVode(cvode_mem, tout, y0, &time, CV_NORMAL);
      if (check_retval(&flag, "CVode", 1)) { sundials_stop(sun_err, "CVode", "Stopping cvsolve, something went wrong in solving the system of ODEs!"); }

      // check whether the this records is sampling or discontinuity using the
      // fourth column of the TCOMB matrix to confirm discontinuity
      if(TCOMB(iout+1,3) == 1){

        // include the discontinuity, i.e. ADD to the solution
        // add the event value to the current value of the state
        for(int i = 0; i < TCOMB.nrow();  i++){

          if(TCOMB(i,1) == tout && TCOMB(i,3) == 1){

            // update y0 - index for y0 needs to be an integer
            int disc_index = static_cast<int>(TCOMB(i,0));
            y0_ptr[disc_index] = y0_ptr[disc_index] + TCOMB(i,2);
          }
        }

        if (flag == CV_SUCCESS) {
          // store results in soln matrix
          soln(iout+1, 0) = time;           // first column is for time
          for (int i = 0; i<y_len; i++){
            soln(iout+1, i+1) = y0_ptr[i];
          }
        }

        // re-initialize the solver
        flag = CVodeReInit(cvode_mem, tout, y0);
        if (check_retval((void *)&flag, "CVodeReInit", 1)) { sundials_stop(sun_err, "CVodeReInit", "Stopping cvsolve, something went wrong in reinitializing the ODE system!"); }

      } else {                                     // store results for the sampling record

        if (flag == CV_SUCCESS) {
          // store results in soln matrix
          soln(iout+1, 0) = time;           // first column is for time
          for (int i = 0; i<y_len; i++){
            soln(iout+1, i+1) = y0_ptr[i];
          }
        }
      }
    }


  }

  // SUNDIALS objects are released by sundials_cleanup on scope exit

  return soln;


}
// cvsolve definition ends ----------------------------------------------------


