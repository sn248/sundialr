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
//   Neither Satyaprakash Nayak nor the names of other
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
#include <iostream>                    // to compare SensType strings

#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <cvodes/cvodes_direct.h>      /* access to CVDls interface            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* definition of ABS */

#include <check_retval.h>
// #include "check_retval.cpp"


// #define Ith(v,i)    NV_Ith_S(v,i-1)         /* i-th vector component i=1..NEQ */

using namespace Rcpp;

// struct to use if R or Rcpp function is input as RHS function----------------
struct rhs_func_sens{
  Function rhs_eqn;
  NumericVector params;
  // std::vector<double> params;
  double rtol;
  NumericVector atol;
};

//struct for tolerance passed for error weight calculations----------------------
// typedef struct {
//   double rtol;
//   NumericVector atol;
// } errorStruct;

// int ewt(N_Vector y, N_Vector w, void *user_data);

// function called by CVodeInit if user inputs R function
int rhs_function_sens(realtype t, N_Vector y, N_Vector ydot, void* user_data){

  // convert y to NumericVector y1
  int y_len = NV_LENGTH_S(y);

  NumericVector y1(y_len);    // filled with zeros
  realtype *y_ptr = N_VGetArrayPointer(y);
  for (int i = 0; i < y_len; i++){
    y1[i] = y_ptr[i];
  }

  NumericVector ydot1(y_len);    // filled with zeros

  if(!user_data){
    stop("Something went wrong, stopping!");
  }

  // cast void pointer to pointer to struct and assign rhs to a Function
  struct rhs_func_sens *my_rhs_fun = (struct rhs_func_sens*)user_data;

  if(my_rhs_fun){
    Function rhs_fun = (*my_rhs_fun).rhs_eqn;                // function
    NumericVector p_values = (*my_rhs_fun).params;           // rate parameters
    // std::vector<double> p_values = (*my_rhs_fun).params;

    if (rhs_fun && (TYPEOF(rhs_fun) == CLOSXP)){
      // use the function to calculate value of RHS ----
      ydot1 = rhs_fun(t, y1, p_values);
    }
    else{
      stop("Something went wrong, stopping!");
    }
  }
  else {
    stop("Something went wrong, stopping!");
  }

  // convert NumericVector ydot1 to N_Vector ydot
  realtype *ydot_ptr = N_VGetArrayPointer(ydot);
  for (int i = 0; i<  y_len; i++){
    ydot_ptr[i] = ydot1[i];
    // Rcout << ydot_ptr[i] << "\n";
  }


  // everything went smoothly
  return(0);
}

// EwtSet Function. Computes the error weights at the current solution.
int ewt(N_Vector y, N_Vector w, void *user_data)
{
  int i;
  realtype yy, ww, rtol;
  NumericVector atol;

  struct rhs_func_sens *my_rhs_fun = (struct rhs_func_sens*)user_data;

  rtol = my_rhs_fun->rtol;
  atol = my_rhs_fun->atol;

  for (i=1; i<=3; i++) {
    yy = NV_Ith_S(y,i-1);
    ww = rtol * SUNRabs(yy) + atol[i-1];
    if (ww <= 0.0) return (-1);
    NV_Ith_S(w,i-1) = 1.0/ww;
  }

  // Rcout << Ith(w,1) << "\n";

  return(0);
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
//'@param SensType Sensitivity Type - allowed values are Staggered (default)", "STG" (for Staggered) or "SIM" (for Simultaneous)
//'@param ErrCon Error Control - allowed values are TRUE or FALSE (default)
//'@example /inst/examples/cvs_Roberts_dns.r
// [[Rcpp::export]]
NumericMatrix cvodes(NumericVector time_vector, NumericVector IC, SEXP input_function,
                     NumericVector Parameters,
                     double reltolerance = 0.0001, NumericVector abstolerance = 0.0001,
                     std::string SensType = "STG", bool ErrCon = 'F'){

  int time_vec_len = time_vector.length();
  int y_len = IC.length();
  int abstol_len = abstolerance.length();

  int flag;

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
  //----------------------------------------------------------------------------
  // // Set the initial conditions-------------------------------------------------
  N_Vector y0 = N_VNew_Serial(y_len);
  realtype *y0_ptr = N_VGetArrayPointer(y0);
  for (int i = 0; i<y_len; i++){
    y0_ptr[i] = IC[i]; // NV_Ith_S(y0, i)
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
  void *cvode_mem;
  cvode_mem = NULL;

  /* Set sensitivity initial conditions */
  int NS = Parameters.length();
  N_Vector *yS;
  yS = N_VCloneVectorArray(NS, y0);

  // Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula
  cvode_mem = CVodeCreate(CV_BDF);
  if (check_retval((void *) cvode_mem, "CVodeCreate", 0)) { stop("Stopping cvodes, cannot allocate memory for CVODES!"); }

  //-- assign user input to the struct based on SEXP type of input_function
  if (!input_function){
    stop("Something is wrong with input function, stopping!");
  }

  // Rcout << TYPEOF(input_function) << "\n";

  switch(TYPEOF(input_function)){

  case CLOSXP:{

    // realtype *params = Parameters.begin();

    struct rhs_func_sens my_rhs_function = {input_function,
                                            Parameters,
                                            // as<std::vector<double> >(Parameters),
                                            reltolerance, abstolerance};

    // setting the user_data in rhs function
    flag = CVodeSetUserData(cvode_mem, (void*)&my_rhs_function);
    if (check_retval(&flag, "CVodeSetUserData", 1)) { stop("Stopping cvodes, something went wrong in setting user data!"); }

    /* Allocate space for CVODES */
    flag = CVodeInit(cvode_mem, rhs_function_sens, T0, y0);
    if (check_retval(&flag, "CVodeInit", 1)) { stop("Stopping cvodes, something went wrong in allocating space for CVODES!"); }

    /* Use private function to compute error weights */
    flag = CVodeWFtolerances(cvode_mem, ewt);
    if (check_retval(&flag, "CVodeSetEwtFn", 1)) { stop("Stopping cvodes, something went wrong in computing error weights!"); }

    /* Create dense SUNMatrix */
    sunindextype y_len_M = y_len;
    SUNMatrix SM = SUNDenseMatrix(y_len_M, y_len_M);
    if (check_retval((void *)SM, "SUNDenseMatrix", 0)) { stop("Stopping cvodes, something went wrong in setting SUNDenseMatrix!"); }

    /* Create dense SUNLinearSolver */
    SUNLinearSolver LS = SUNLinSol_Dense(y0, SM);
    if (check_retval((void *)LS, "SUNLinSol_Dense", 0)) { stop("Stopping cvodes, something went wrong in setting Linear Solver!"); }

    /* Attach the matrix and linear solver */
    flag = CVodeSetLinearSolver(cvode_mem, LS, SM);
    if (check_retval(&flag, "CVodeSetLinearSolver", 1)) { stop("Stopping cvodes, something went wrong in attaching SUNDenseMatrix and Linear Solver!"); }


    if (check_retval((void *)yS, "N_VCloneVectorArray", 0)) { stop("Stopping cvodes, something went wrong in setting Sensitivity Array!"); }
    for (int is=0;is<NS;is++) N_VConst(RCONST(0.0), yS[is]);

    /* Call CVodeSensInit1 to activate forward sensitivity computations
     and allocate internal memory for COVEDS related to sensitivity
     calculations. Computes the right-hand sides of the sensitivity
     ODE, one at a time */
    // (int) flag to set sensitivity solution method - see CVodeSensInit1
    int ism = CV_STAGGERED;
    if (SensType.compare("SIM") == 0) ism = CV_SIMULTANEOUS;
    flag = CVodeSensInit1(cvode_mem, NS, ism, NULL, yS);
    if(check_retval(&flag, "CVodeSensInit", 1)) { stop("Stopping cvodes, something went wrong in calculating Sensitivities!"); }

    /* Call CVodeSensEEtolerances to estimate tolerances for sensitivity
     variables based on the rolerances supplied for states variables and
     the scaling factor pbar */
    flag = CVodeSensEEtolerances(cvode_mem);
    if(check_retval(&flag, "CVodeSensEEtolerances", 1)) { stop("Stopping cvodes, something went wrong in estimating tolerances for sensitivities!"); }

    /* Call CVodeSetSensParams to specify problem parameter information for
     sensitivity calculations */
    // struct rhs_func_sens *ptr = &my_rhs_function;        // struct UserData storing my data
    // p = (my_rhs_function.params).begin();
    // Rcout << (my_rhs_function.params).begin() << "\n";
    // Rcout << &(ptr->params[0]);
    flag = CVodeSetSensParams(cvode_mem, (my_rhs_function.params).begin(), Parameters.begin(), NULL);  // double *y = x.begin()
    if (check_retval(&flag, "CVodeSetSensParams", 1)) { stop("Stopping cvodes, something went wrong in setting Sensitivity Parameters!"); }

    // First row for initial conditions, First column is for time
    int y_len_1 = y_len + 1;
    NumericMatrix soln(Dimension(time_vec_len,y_len_1));

    // fill the first row of soln matrix with Initial Conditions
    soln(0,0) = time_vector[0];   // get the first time value
    for(int i = 0; i<y_len; i++){
      soln(0,i+1) = IC[i];
    }

    NumericMatrix sens(Dimension(time_vec_len, (NS * y_len) + 1));
    // Store the Sensitivity Results
    // Sensitivity of each entitiy w.r.t. parameters is zero at initial time
    sens(0,0) = time_vector[0];           // get the first time value
    for(int i = 0; i<y_len; i++){
      sens(0,i+1) = 0.0;
    }

    realtype tout;  // For output times
    realtype *yS_ptr = NULL;
    for(int iout = 0; iout < NOUT-1; iout++) {

      // output times start from the index after initial time
      tout = time_vector[iout+1];

      flag = CVode(cvode_mem, tout, y0, &time, CV_NORMAL);
      if (check_retval(&flag, "CVode", 1)) { stop("Stopping cvodes, something went wrong in solving the system using CVODE!"); break; } // Something went wrong in solving it!
      if (flag == CV_SUCCESS) {

        // store results in soln matrix
        soln(iout+1, 0) = time;           // first column is for time
        for (int i = 0; i<y_len; i++){
          soln(iout+1, i+1) = y0_ptr[i];
        }
      }

      flag = CVodeGetSens(cvode_mem, &time, yS);
      if (check_retval(&flag, "CVodeGetSens", 1)) { stop("Stopping cvodes, something went wrong in calculating Sensitivities!"); break; }

      if (flag == CV_SUCCESS) {
        for (int i = 0; i < Parameters.length(); i++){
          yS_ptr = N_VGetArrayPointer(yS[i]);
          for(int j = 0; j < y_len; j++){
            // store results in soln matrix
            sens(iout+1, 0) = time;           // first column is for time
            sens(iout+1, (NS*i)+j+1) = yS_ptr[j];
          }
        }
      }

    }

    N_VDestroy(y0);                          /* Free y vector */
    N_VDestroyVectorArray(yS, NS);           /* Free yS vector */
    // free(data);                           /* Free user data */
    CVodeFree(&cvode_mem);                   /* Free CVODES memory */
    SUNLinSolFree(LS);                       /* Free the linear solver memory */
    SUNMatDestroy(SM);                       /* Free the matrix memory */

    return sens;
    break;

  }

  default: {
    stop("Incorrect input function type - input function can be an R or Rcpp function");
  }

  }

}
