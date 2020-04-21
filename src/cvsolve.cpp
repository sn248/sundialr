// sundialr is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// any later version.
//
// sundialr is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with sundialr.  If not, see <http://www.gnu.org/licenses/>.


#include <Rcpp.h>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>    /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h>   /* definition of type realtype */
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

#include <check_retval.h>
#include <rhs_func.h>

using namespace Rcpp;

//------------------------------------------------------------------------------
//'cvsolve
//'
//'solver to solve stiff ODEs with discontinuties
//'@param time_vector time vector
//'@param IC Initial Conditions
//'@param input_function Right Hand Side function of ODEs
//'@param Parameters Parameters input to ODEs
//'@param events Discontinuities in the solution (a DataFrame, default value is NULL)
//'@param Jacobian User-supplied Jacobian(a matrix, default value is NULL)
//'@param reltolerance Relative Tolerance (a scalar, default value  = 1e-04)
//'@param abstolerance Absolute Tolerance (a scalar or vector with length equal to ydot, default = 1e-04)
//'@param constraints By default all the solution values are constrained to be non-negative
//'@example /inst/examples/cv_Roberts_dns.r
// [[Rcpp::export]]
NumericMatrix cvsolve(NumericVector time_vector, NumericVector IC, SEXP input_function,
                    NumericVector Parameters,
                    NumericMatrix events = R_NilValue,
                    NumericMatrix Jacobian = R_NilValue,
                    double reltolerance = 0.0001,
                    NumericVector abstolerance = 0.0001){

  int time_vec_len = time_vector.length();
  int y_len = IC.length();
  int abstol_len = abstolerance.length();

  int flag;
  realtype reltol = reltolerance;

  realtype T0 = RCONST(time_vector[0]);     //RCONST(0.0);  // Initial Time

  double time;
  int NOUT = time_vec_len;

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

  // Set the initial conditions-------------------------------------------------
  N_Vector y0 = N_VNew_Serial(y_len);
  realtype *y0_ptr = N_VGetArrayPointer(y0);
  for (int i = 0; i<y_len; i++){
    y0_ptr[i] = IC[i]; // NV_Ith_S(y0, i)
  }

  // Set constraints to all 1's for nonnegative solution values.----------------
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

    // // User-Supplied Jacobian-----------------------------------------------------
    // if(Jacobian != R_NilValue){
    //   // Jacobian should be size y_len * y_len
    //   if(nrow(Jacobian) != y_len || ncol(Jacobian) != y_len){
    //     stop("Jacobian should be a square matrix of size equal to length of IC \n")
    //   } else {
    //     // Set the user-supplied Jacobian
    //     flag = CVodeSetJacFn(cvode_mem, Jacobian);
    //     if(check_retval(&flag, "CVodeSetJacFn", 1)) { stop("Stopping cvsolve, something went wrong in setting Jacobian!"); }
    //   }
    // }

    // Call CVodeSetConstraints to initialize constraints
    flag = CVodeSetConstraints(cvode_mem, constraints);
    if(check_retval(&flag, "CVodeSetConstraints", 1)) { stop("Stopping cvsolve, something went wrong in setting constraints!"); }

    // NumericMatrix to store results - filled with 0.0
    // Call CVodeInit to initialize the integrator memory and specify the
    // user's right hand side function in y'=f(time,y),
    // // the inital time T0, and the initial dependent variable vector y.
    realtype tout;  // For output times

    int y_len_1 = y_len + 1; // remove later
    NumericMatrix soln(Dimension(time_vec_len,y_len_1));  // remove later

   // fill the first row of soln matrix with Initial Conditions
    soln(0,0) = time_vector[0];   // get the first time value
    for(int i = 0; i<y_len; i++){
      soln(0,i+1) = IC[i];
    }

    // use the events matrix to get discontinuities in solution
    // if(events != R_NilValue){
    //   if(ncol(events) != 3) { stop("events should have 3 columns, first column is for ID of the discontinuous species, second column for Time\nand third column is for value of solution at the discontinuity"); }
    //
    //   // Combine the sampling and dosing times
    // }


    for(int iout = 0; iout < NOUT-1; iout++) {

      // output times start from the index after initial time
      tout = time_vector[iout+1];

      flag = CVode(cvode_mem, tout, y0, &time, CV_NORMAL);

      if (check_retval(&flag, "CVode", 1)) { stop("Stopping cvsolve, something went wrong in solving the system of ODEs!"); break; } // Something went wrong in solving it!
      if (flag == CV_SUCCESS) {

        // store results in soln matrix
        soln(iout+1, 0) = time;           // first column is for time
        for (int i = 0; i<y_len; i++){
          soln(iout+1, i+1) = y0_ptr[i];
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

}
// cvsolve definition ends ----------------------------------------------------
