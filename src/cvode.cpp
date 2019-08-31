
#include <Rcpp.h>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>    /* serial N_Vector types, fcts., macros */
// #include <cvode/cvode_direct.h>         /* prototype for CVDense */
// #include <sundials/sundials_dense.h>   /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h>   /* definition of type realtype */
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

#include <check_retval.h>
// #include "check_retval.cpp"

#define Ith(v,i)    NV_Ith_S(v,i-1)         /* i-th vector component i=1..NEQ */

using namespace Rcpp;

// struct to use if R or Rcpp function is input as RHS function
struct rhs_func{
  Function rhs_eqn;
  NumericVector params;
};

// function called by CVodeInit if user inputs R function
int rhs_function(realtype t, N_Vector y, N_Vector ydot, void* user_data){

  // convert y to NumericVector y1
  int y_len = NV_LENGTH_S(y);

  NumericVector y1(y_len);    // filled with zeros
  // realtype *y_ptr = N_VGetArrayPointer(y);
  for (int i = 0; i < y_len; i++){
    y1[i] = Ith(y,i+1); // y_ptr[i];
  }

  NumericVector ydot1(y_len);    // filled with zeros

  if(!user_data){
    stop("Something went wrong, stopping!");
  }

  // cast void pointer to pointer to struct and assign rhs to a Function
  struct rhs_func *my_rhs_fun = (struct rhs_func*)user_data;

  if(my_rhs_fun){
    Function rhs_fun = (*my_rhs_fun).rhs_eqn;           // function
    NumericVector p_values = (*my_rhs_fun).params;      // rate parameters

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
  }

  // everything went smoothly
  return(0);
}
//---RHS function definition ends ----------------------------------------------

//------------------------------------------------------------------------------
//'cvode
//'
//' CVODE solver to solve stiff ODEs
//'@param time_vector time vector
//'@param IC Initial Conditions
//'@param input_function Right Hand Side function of ODEs
//'@param Parameters Parameters input to ODEs
//'@param reltolerance Relative Tolerance (a scalar)
//'@param abstolerance Absolute Tolerance (a vector with length equal to ydot)
//'@example /inst/examples/cv_Roberts_dns.r
// [[Rcpp::export]]
NumericMatrix cvode(NumericVector time_vector, NumericVector IC, SEXP input_function,
                    NumericVector Parameters, NumericMatrix Events = R_NilValue,
                    double reltolerance = 1e-04, NumericVector abstolerance = 1e-04){

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
  //----------------------------------------------------------------------------
  // // Set the initial conditions-------------------------------------------------
  N_Vector y0 = N_VNew_Serial(y_len);
  realtype *y0_ptr = N_VGetArrayPointer(y0);
  for (int i = 0; i<y_len; i++){
    y0_ptr[i] = IC[i]; // NV_Ith_S(y0, i)
  }
  // //----------------------------------------------------------------------------
  void *cvode_mem;
  cvode_mem = NULL;

  // // Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula
  cvode_mem = CVodeCreate(CV_BDF);
  if (check_retval((void *) cvode_mem, "CVodeCreate", 0)) { stop("Stopping cvode!"); }

  //-- assign user input to the struct based on SEXP type of input_function

  if (!input_function){
    stop("Something is wrong with input function, stopping!");
  }
  switch(TYPEOF(input_function)){

  case CLOSXP:{

    struct rhs_func my_rhs_function = {input_function, Parameters};

    // setting the user_data in rhs function
    flag = CVodeSetUserData(cvode_mem, (void*)&my_rhs_function);
    if (check_retval(&flag, "CVodeSetUserData", 1)) { stop("Stopping cvode, something went wrong in setting user data!"); }

    flag = CVodeInit(cvode_mem, rhs_function, T0, y0);
    if (check_retval(&flag, "CVodeInit", 1)) { stop("Stopping cvode, something went wrong in initializing CVODE!"); }

    // Call CVodeSVtolerances to specify the scalar relative tolerance and vector absolute tol
    flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
    if (check_retval(&flag, "CVodeSVtolerances", 1)) { stop("Stopping cvode, something went wrong in setting solver tolerances!"); }

    //--required in the new version ----------------------------------------------
    // Create dense SUNMatrix for use in linear solves
    sunindextype y_len_M = y_len;
    SUNMatrix SM = SUNDenseMatrix(y_len_M, y_len_M);
    if(check_retval((void *)SM, "SUNDenseMatrix", 0)) { stop("Stopping cvode, something went wrong in setting the dense matrix!"); }

    // Create dense SUNLinearSolver object for use by CVode
    SUNLinearSolver LS = SUNLinSol_Dense(y0, SM);
    if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) { stop("Stopping cvode, something went wrong in setting the linear solver!"); }

    // Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode
    flag = CVodeSetLinearSolver(cvode_mem, LS, SM);
    if(check_retval(&flag, "CVDlsSetLinearSolver", 1)) { stop("Stopping cvode, something went wrong in setting the linear solver!"); }
    // NumericMatrix to store results - filled with 0.0



    // Call CVodeInit to initialize the integrator memory and specify the
    // user's right hand side function in y'=f(time,y),
    // // the inital time T0, and the initial dependent variable vector y.
    realtype tout;  // For output times

    int y_len_1 = y_len + 1; // remove later
    NumericMatrix soln(Dimension(time_vec_len,y_len_1));  // remove later

    if(Events == R_NilValue){

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

        flag = CVode(cvode_mem, tout, y0, &time, CV_NORMAL);

        if (check_retval(&flag, "CVode", 1)) { stop("Stopping CVODE, something went wrong in solving the system of ODEs!"); break; } // Something went wrong in solving it!
        if (flag == CV_SUCCESS) {

          // store results in soln matrix
          soln(iout+1, 0) = time;           // first column is for time
          for (int i = 0; i<y_len; i++){
            soln(iout+1, i+1) = y0_ptr[i];
          }
        }
      }
    } else {                                 // When Events is not R_NilValue

      double T1 = 4.0;
      double T2 = 8.0;

      Rcout << "\nDiscontinuity in solution\n\n";

      // for(int iout = 0; iout < NOUT-1; iout++) {
      //
      //
      // }
      //
      // /* set TSTOP (max time solution proceeds to) - this is not required */
      // flag = CVodeSetStopTime(cvode_mem, T1);
      // if (check_retval((void *)&flag, "CVodeSetStopTime", 1)) {stop("Stopping CVODE, something went wrong!");};
      //
      // while(time < T1){
      //
      //   /* set TSTOP (max time solution proceeds to) - this is not required */
      //   flag = CVodeSetStopTime(cvode_mem, t1);
      //   if (check_retval((void *)&ret, "CVodeSetStopTime", 1)) {stop("Stopping CVODE, something went wrong!");};
      // }



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

// //--- CVODE definition ends ----------------------------------------------------
