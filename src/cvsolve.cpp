//   Copyright (c) 2020, Satyaprakash Nayak
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
//   Neither the name of the <ORGANIZATION> nor the names of its
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
#include "sortTimes.cpp"

using namespace Rcpp;
using namespace arma;

//------------------------------------------------------------------------------
//'cvsolve
//'
//'CVSOLVE solver to solve stiff ODEs with discontinuties
//'@param time_vector time vector
//'@param IC Initial Conditions
//'@param input_function Right Hand Side function of ODEs
//'@param Parameters Parameters input to ODEs
//'@param Events Discontinuities in the solution (a DataFrame, default value is NULL)
//'@param reltolerance Relative Tolerance (a scalar, default value  = 1e-04)
//'@param abstolerance Absolute Tolerance (a scalar or vector with length equal to ydot, default = 1e-04)
//'@example /inst/examples/cvsolve_1D.r
// [[Rcpp::export]]
NumericMatrix cvsolve(NumericVector time_vector, NumericVector IC,
                      SEXP input_function,
                      NumericVector Parameters,
                      Nullable<DataFrame> Events = R_NilValue,
                      double reltolerance = 0.0001,
                      NumericVector abstolerance = 0.0001){

  int y_len = IC.length();
  int NSTATES = IC.length();
  int abstol_len = abstolerance.length();

  int flag;
  realtype reltol = reltolerance;

  realtype T0 = RCONST(time_vector[0]);     //RCONST(0.0);  // Initial Time

  double time;

  // Set the vector absolute tolerance -----------------------------------------
  // abstol must be same length as IC
  N_Vector abstol = N_VNew_Serial(abstol_len);
  realtype *abstol_ptr = N_VGetArrayPointer(abstol);
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
  else if(abstol_len != 1 || abstol_len != y_len){
    stop("Absolute tolerance must be a scalar or a vector of same length as IC \n");
  }

  // If Events is not NULL, change IC and generate a combined dataset-----------
  // Combine the time vector and Events vector into a single dataframe
  // If Events is null, then TCOMB is same as time_vec (i.e, is a vector)
  // If Events is not null, then TCOMB is a matrix with 4 columns
  // Column 1 - cpp index of state with discontinuity, -1 for all sampling time points
  // Column 2 - Time of discontinuity or sampling
  // Column 3 - Value at the time of discontinuity, 0 for sampling
  // Column 4 - 0 for sampling and 1 for discontinuity time point
  // If IC of a state is specified by both IC and Events, then the IC value is
  // overwritten by the Events value
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

    DataFrame Events_DF(Events);
    TCOMB = sorted_times(Events_DF, time_vector, NSTATES);  // get the sorted  Combined time matrix
    // Decrease State index by 1 to convert R state index to internal cpp index, all sampling times get state index of -1
    for(int i = 0; i < TCOMB.nrow(); i++){
      TCOMB(i,0) = TCOMB(i,0) - 1;
    }
    // Change IC such that IC at time = 0 are overwritten by values in TCOMB at t = 0
    for(int i = 0; i < TCOMB.nrow();  i++){
      if(TCOMB(i,1) == 0 && TCOMB(i,3) == 1){
        ICchanged(TCOMB(i,0)) = TCOMB(i,2);
      }
    }
  }

  int NOUT = TCOMB.nrow();

  // Set the initial conditions-------------------------------------------------
  N_Vector y0 = N_VNew_Serial(y_len);
  realtype *y0_ptr = N_VGetArrayPointer(y0);
  for (int i = 0; i<y_len; i++){
    y0_ptr[i] = ICchanged[i];
  }

  // Set constraints to all 1's for non-negative solution values.----------------
  N_Vector constraints = N_VNew_Serial(y_len);
  N_VConst(1, constraints);

  // Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula
  void *cvode_mem;
  cvode_mem = NULL;

  cvode_mem = CVodeCreate(CV_BDF);
  if (check_retval((void *) cvode_mem, "CVodeCreate", 0)) { stop("Stopping cvsolve!"); }

  //-- assign user input to the struct based on SEXP type of input_function
  if (!input_function){
    stop("Something is wrong with input function, stopping!");
  }
  switch(TYPEOF(input_function)){

  case CLOSXP:{

    struct rhs_func my_rhs_function = {input_function, Parameters};

    // setting the user_data in rhs function
    flag = CVodeSetUserData(cvode_mem, (void*)&my_rhs_function);
    if (check_retval(&flag, "CVodeSetUserData", 1)) { stop("Stopping cvsolve, something went wrong in setting user data!"); }

    flag = CVodeInit(cvode_mem, rhs_function, T0, y0);
    if (check_retval(&flag, "CVodeInit", 1)) { stop("Stopping cvsolve, something went wrong in initializing CVODE!"); }

    // Call CVodeSVtolerances to specify the scalar relative tolerance and vector absolute tol
    flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
    if (check_retval(&flag, "CVodeSVtolerances", 1)) { stop("Stopping cvsolve, something went wrong in setting solver tolerances!"); }

    // Create dense SUNMatrix for use in linear solves
    sunindextype y_len_M = y_len;
    SUNMatrix SM = SUNDenseMatrix(y_len_M, y_len_M);
    if(check_retval((void *)SM, "SUNDenseMatrix", 0)) { stop("Stopping cvsolve, something went wrong in setting the dense matrix!"); }

    // Create dense SUNLinearSolver object for use by CVode
    SUNLinearSolver LS = SUNLinSol_Dense(y0, SM);
    if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) { stop("Stopping cvsolve, something went wrong in setting the linear solver!"); }

    // Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode
    flag = CVodeSetLinearSolver(cvode_mem, LS, SM);
    if(check_retval(&flag, "CVDlsSetLinearSolver", 1)) { stop("Stopping cvsolve, something went wrong in setting the linear solver!"); }

    // Call CVodeSetConstraints to initialize constraints
    flag = CVodeSetConstraints(cvode_mem, constraints);
    if(check_retval(&flag, "CVodeSetConstraints", 1)) { stop("Stopping cvsolve, something went wrong in setting constraints!"); }

    realtype tout;  // For output times

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

      realtype tprev = TCOMB(iout, 1);
      // output times start from the index after initial time
      tout = TCOMB(iout + 1, 1);

      // tout 0 is taken care of by initial conditions, if tout equal to previous time
      // then is taken care by the dosing record at the same time point
      if (tout == 0.0 || tout == tprev){
        continue;

      } else {

        // integrate upto the next time point (whether a sampling time or a discontinuity)
        flag = CVode(cvode_mem, tout, y0, &time, CV_NORMAL);
        if (check_retval(&flag, "CVode", 1)) { stop("Stopping cvsolve, something went wrong in solving the system of ODEs!"); break; }

        // check whether the this records is sampling or discontinuity using the
        // fourth column of the TCOMB matrix to confirm discontinuity
        if(TCOMB(iout+1,3) == 1){

          // include the discontinuity, i.e. ADD to the solution
          // Change y0 such that IC at time = 0 are overwritten by values in TCOMB at t = 0
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
          if (check_retval((void *)&flag, "CVodeReInit", 1)) { stop("Stopping cvsolve, something went wrong in reinitializing the ODE system!"); break; }

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

    // free the vectors
    N_VDestroy(y0);
    N_VDestroy(abstol);

    // free integrator memory
    CVodeFree(&cvode_mem);
    // free the linear solver memory
    SUNLinSolFree(LS);
    // Free the matrix memory
    SUNMatDestroy(SM);

    return soln;
    break;
  }

  default: {
    stop("Incorrect input function type - input function can be an R or Rcpp function");
  }

  }

  // return TCOMB;

}
// cvsolve definition ends ----------------------------------------------------


