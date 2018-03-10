
#include <Rcpp.h>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>    /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_direct.h>         /* prototype for CVDense */
// #include <sundials/sundials_dense.h>   /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h>   /* definition of type realtype */
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

#include <check_flag.h>


using namespace Rcpp;

// #define T0    RCONST(0.0)             // initial time
// #define T1    RCONST(0.4)                // first output time
#define TMULT RCONST(10.0)               // output time factor
// #define NOUT  12

//------------------------------------------------------------------------------
typedef int (*funcPtr)(realtype time, N_Vector y, N_Vector ydot, void *user_data);
//------------------------------------------------------------------------------

//'cvode
//'
//' CVODE solver to solve stiff ODEs
//'@param time_vec time vector
//'@param IC Initial Conditions
//'@param xpsexp External pointer to RHS function
//'@param reltolerance Relative Tolerance (a scalar)
//'@param abstolerance Absolute Tolerance (a vector with length equal to ydot)
// [[Rcpp::export]]
NumericMatrix cvode (NumericVector time_vec, NumericVector IC, SEXP xpsexp,
           double reltolerance, NumericVector abstolerance){

  int flag;
  realtype reltol = reltolerance;
  N_Vector abstol, y0;
  y0 = abstol = NULL;

  // realtype T0 = RCONST(0.0);
  // realtype tout = T1;
  realtype T0 = RCONST(time_vec[0]); //RCONST(0.0);  // Initial Time
  realtype tout = RCONST(time_vec[1]); // T1;        // First output time
  realtype iout;
  SUNMatrix SM;
  SUNLinearSolver LS;
  SM = NULL;
  LS = NULL;

  double time;
  int NOUT = time_vec.length() - 1;

  // Set the vector absolute tolerance -----------------------------------------
  // abstol must be same length as IC
  abstol = N_VNew_Serial(IC.length());
  if(abstolerance.length() == 1){
    // if a scalar is provided - use it to make a vector with same values
    for (int i = 0; i<IC.length(); i++){
      NV_Ith_S(abstol, i) = abstolerance[0];
    }
  }
  else if (abstolerance.length() == IC.length()){

    for (int i = 0; i<abstolerance.length(); i++){
      NV_Ith_S(abstol, i) = abstolerance[i];
    }
  }
  else if(abstolerance.length() != 1 || abstolerance.length() != IC.length()){

    stop("Absolute tolerance must be a scalar or a vector of same length as IC \n");
  }
  //----------------------------------------------------------------------------

  // Set the initial conditions-------------------------------------------------
  y0 = N_VNew_Serial(IC.length());
  for (int i = 0; i<IC.length(); i++){
    NV_Ith_S(y0, i) = IC[i];
  }
  //----------------------------------------------------------------------------

  // void pointer, can point to any data type -  will have to be cast before use
  void *cvode_mem;
  cvode_mem = NULL;

  // Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula
  // and the use of a Newton iteration
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *) cvode_mem, "CVodeCreate", 0)) {
    return (1);
  }

  XPtr<funcPtr> xpfun(xpsexp);
  funcPtr fun = *xpfun;

  flag = CVodeInit(cvode_mem, fun, T0, y0);
  if (check_flag(&flag, "CVodeInit", 1)) {
    return (1);
  }

  // Call CVodeSVtolerances to specify the scalar relative tolerance and vector absolute tol
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_flag(&flag, "CVodeSVtolerances", 1)) {
    return (1);
  }

  //------not required now for the new version ---------------------------------
  // Call CVDense to specify the CVDENSE dense linear solver
  // flag = CVDense(cvode_mem, 3);
  // if (check_flag(&flag, "CVDense", 1)) {
  //   return (1);
  // }
  //----------------------------------------------------------------------------

  //--required in the new version ----------------------------------------------
  // Create dense SUNMatrix for use in linear solves
  SM = SUNDenseMatrix(IC.length(), IC.length());
  if(check_flag((void *)SM, "SUNDenseMatrix", 0)) return (1);

  // Create dense SUNLinearSolver object for use by CVode
  LS = SUNDenseLinearSolver(y0, SM);
  if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);

  // Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode
  flag = CVDlsSetLinearSolver(cvode_mem, LS, SM);
  if(check_flag(&flag, "CVDlsSetLinearSolver", 1)) return(1);

  // NumericMatrix to store results - filled with 0.0
  // First row for initial conditions, First column is for time
  NumericMatrix soln = NumericMatrix(time_vec.length()+1, IC.length() + 1);

  // fill the first row of soln matrix with Initial Conditions
  soln(0,0) = time_vec[0];
  for(int i = 0; i<IC.length(); i++){
    soln(0,i+1) = IC[i];
  }

  // Call CVodeInit to initialize the integrator memory and specify the
  // user's right hand side function in y'=f(time,y),
  // the inital time T0, and the initial dependent variable vector y.
  iout = 0;
  while (1) {

    flag = CVode(cvode_mem, tout, y0, &time, CV_NORMAL);

    if (check_flag(&flag, "CVode", 1)) { break; } // Something went wrong in solving it!
    if (flag == CV_SUCCESS) {

      // store results in soln matrix
      soln(iout+1, 0) = time;           // first column is for time
      for (int i = 0; i<IC.length(); i++){
        soln(iout+1, i+1) = NV_Ith_S(y0, i);
      }

      // Rcout << time <<  "\t" << NV_Ith_S(y0,0) << "\t" << NV_Ith_S(y0,1) << "\t" << NV_Ith_S(y0,2) << '\n';

      iout++;

      if (iout < NOUT){
        tout += time_vec[iout+1] - time_vec[iout];   //TMULT; // tout *= TMULT orginally
        }

    }

    if (iout > NOUT) { break; }
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
  // return(0); // everything went fine
}

int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL)
  {
    // fprintf(stderr, "\n SUNDIALS_ERROR: %s() failed - returned NULL pointer \n\n", funcname);
    Rprintf("SUNDIALS_ERROR: %s() failed - returned NULL pointer", funcname);
    return(1);
  }

  /* Check if flag < 0 */
  else if (opt == 1)
  {
    errflag = (int *) flagvalue;
    if (*errflag < 0)
    {
      //fprintf(stderr, "\n SUNDIALS_ERROR: %s() failed with flag = %d \n\n", funcname, *errflag);
      Rprintf("SUNDIALS_ERROR: %s() failed with flag = %d", funcname, *errflag);
      return(1);
    }
  }
  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    // fprintf(stderr, "\n MEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    Rprintf("MEMORY_ERROR: %s() failed - returned NULL pointer", funcname);
    return (1);
  }

  return (0); // so returns 0 if everything is okay

}
