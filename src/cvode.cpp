
#include <Rcpp.h>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>    /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_direct.h>         /* prototype for CVDense */
// #include <sundials/sundials_dense.h>   /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h>   /* definition of type realtype */
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

#include <check_flag.h>
#include "check_flag.cpp"

using namespace Rcpp;

//-- typedefs for RHS function pointer input in from R -------------------------
// typedef NumericVector (*funcPtr) (double t, NumericVector y);
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-- user data is passed to CVODE in a struct, one component is the pointer to
//-- RHS function---------------------------------------------------------------
// struct to use is pointer to Rcpp is input as RHS function
// struct rhs_xptr{
//   SEXP rhs_eqn;
// };

// struct to use if R or Rcpp function is input as RHS function
struct rhs_func{
  Function rhs_eqn;
};
//------------------------------------------------------------------------------

//---RHS function that is an input to CVODE, takes user input of RHS in user_data

// function called by CVodeInit if user inputs Rcpp function pointer
// int rhs_pointer(realtype t, N_Vector y, N_Vector ydot, void* user_data){
//
//   // convert y to NumericVector y1
//   int y_len = NV_LENGTH_S(y);
//
//   NumericVector y1(y_len);    // filled with zeros
//   realtype *y_ptr = N_VGetArrayPointer(y);
//   for (int i = 0; i < y_len; i++){
//     y1[i] = y_ptr[i];
//   }
//
//   // convert ydot to NumericVector ydot1
//   int ydot_len = NV_LENGTH_S(ydot);
//
//   NumericVector ydot1(ydot_len);    // filled with zeros
//
//   // cast void pointer to pointer to struct and assign rhs to a SEXP
//   Rprintf("Reached in function\n");
//   struct rhs_xptr *my_rhs_ptr = NULL;
//   my_rhs_ptr = (struct rhs_xptr*)user_data;
//   SEXP rhs_fun_sexp = (*my_rhs_ptr).rhs_eqn;
//
//   Rprintf("type of rhs_fun_sexp is %d\n", TYPEOF(rhs_fun_sexp));
//
//   // use function pointer to get the derivatives
//   XPtr<funcPtr> rhs_fun_xptr(rhs_fun_sexp);
//   funcPtr rhs_fun = *rhs_fun_xptr;
//
//   // Rprintf("type of rhs_fun is %d\n", TYPEOF(rhs_fun));
//   // use the function to calculate value of RHS ----
//   ydot1 = rhs_fun(t, y1);
//
//   // convert NumericVector ydot1 to N_Vector ydot
//   realtype *ydot_ptr = N_VGetArrayPointer(ydot);
//   for (int i = 0; i<ydot1.length(); i++){
//     ydot_ptr[i] = ydot1[i];
//   }
//
//   // everything went smoothly
//   return(0);
// }

// function called by CVodeInit if user inputs R function
int rhs_function(realtype t, N_Vector y, N_Vector ydot, void* user_data){

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
  struct rhs_func *my_rhs_fun = (struct rhs_func*)user_data;

  if(my_rhs_fun){
    Function rhs_fun = (*my_rhs_fun).rhs_eqn;
    if (rhs_fun && (TYPEOF(rhs_fun) == CLOSXP)){
      // use the function to calculate value of RHS ----
      ydot1 = rhs_fun(t, y1);
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
//'@param reltolerance Relative Tolerance (a scalar)
//'@param abstolerance Absolute Tolerance (a vector with length equal to ydot)
//'@example /inst/examples/cv_Roberts_dns.r
// [[Rcpp::export]]
NumericMatrix cvode(NumericVector time_vector, NumericVector IC, SEXP input_function,
                    double reltolerance, NumericVector abstolerance){

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
  // // and the use of a Newton iteration
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *) cvode_mem, "CVodeCreate", 0)) { stop("Stopping cvode!"); }

  //-- assign user input to the struct based on SEXP type of input_function

  if (!input_function){
    stop("Something is wrong with input function, stopping!");
  }
  switch(TYPEOF(input_function)){

  case CLOSXP:{

    struct rhs_func my_rhs_function = {input_function};
    // setting the user_data in rhs function
    flag = CVodeSetUserData(cvode_mem, (void*)&my_rhs_function);
    if (check_flag(&flag, "CVodeSetUserData", 1)) { stop("Stopping cvode!"); }

    flag = CVodeInit(cvode_mem, rhs_function, T0, y0);
    if (check_flag(&flag, "CVodeInit", 1)) { stop("Stopping cvode!"); }

    // Call CVodeSVtolerances to specify the scalar relative tolerance and vector absolute tol
    flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
    if (check_flag(&flag, "CVodeSVtolerances", 1)) { stop("Stopping cvode!"); }

    //--required in the new version ----------------------------------------------
    // Create dense SUNMatrix for use in linear solves
    sunindextype y_len_M = y_len;
    SUNMatrix SM = SUNDenseMatrix(y_len_M, y_len_M);
    if(check_flag((void *)SM, "SUNDenseMatrix", 0)) { stop("Stopping cvode!"); }

    // Create dense SUNLinearSolver object for use by CVode
    SUNLinearSolver LS = SUNDenseLinearSolver(y0, SM);
    if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) { stop("Stopping cvode!"); }

    // Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode
    flag = CVDlsSetLinearSolver(cvode_mem, LS, SM);
    if(check_flag(&flag, "CVDlsSetLinearSolver", 1)) { stop("Stopping cvode!"); }
    // NumericMatrix to store results - filled with 0.0

    // First row for initial conditions, First column is for time
    int y_len_1 = y_len + 1;
    NumericMatrix soln(Dimension(time_vec_len,y_len_1));

    // fill the first row of soln matrix with Initial Conditions
    soln(0,0) = time_vector[0];   // get the first time value
    for(int i = 0; i<y_len; i++){
      soln(0,i+1) = IC[i];
    }

    // Call CVodeInit to initialize the integrator memory and specify the
    // user's right hand side function in y'=f(time,y),
    // // the inital time T0, and the initial dependent variable vector y.
    realtype tout;  // For output times
    for(int iout = 0; iout < NOUT-1; iout++) {

      // output times start from the index after initial time
      tout = time_vector[iout+1];

      flag = CVode(cvode_mem, tout, y0, &time, CV_NORMAL);

      if (check_flag(&flag, "CVode", 1)) { break; } // Something went wrong in solving it!
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

  // case EXTPTRSXP: {
  //   Rprintf("Reached in EXTPTRSXP\n");
  //   Rprintf("Type of input_function is %d\n", TYPEOF(input_function));
  //   struct rhs_xptr my_rhs_xptr = {input_function};
  //   // setting the user_data in rhs function
  //   flag = CVodeSetUserData(cvode_mem, (void*)&my_rhs_xptr);
  //   if (check_flag(&flag, "CVodeSetUserData", 1)) { stop("Stopping cvode!"); }
  //
  //   flag = CVodeInit(cvode_mem, rhs_pointer, T0, y0);
  //   if (check_flag(&flag, "CVodeInit", 1)) { stop("Stopping cvode!"); }
  //
  //   break;
  // }

  default: {
    stop("Incorrect input function type - input function can be an R or Rcpp function");
    }

  }



}

// //--- CVODE definition ends ----------------------------------------------------
