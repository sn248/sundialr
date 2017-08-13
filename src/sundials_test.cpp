// [[Rcpp::depends(BH)]]
#include <Rcpp.h>

#include <iostream>
#include <vector>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>    /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>         /* prototype for CVDense */
#include <sundials/sundials_dense.h>   /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h>   /* definition of type realtype */

using namespace Rcpp;

// #define T0    RCONST(0.0)             // initial time
#define T1    RCONST(0.4)                // first output time
#define TMULT RCONST(10.0)               // output time factor
#define NOUT  12

// Check function return value...
//   opt == 0 means SUNDIALS function allocates memory so check if
//            returned NULL pointer
//   opt == 1 means SUNDIALS function returns a flag so check if
//            flag >= 0
//   opt == 2 means function allocates memory so check if returned
//            NULL pointer
//
int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL)
  {
    fprintf(stderr, "\n SUNDIALS_ERROR: %s() failed - returned NULL pointer \n\n", funcname);
    return(1);
  }

  /* Check if flag < 0 */
  else if (opt == 1)
  {
    errflag = (int *) flagvalue;
    if (*errflag < 0)
    {
      fprintf(stderr, "\n SUNDIALS_ERROR: %s() failed with flag = %d \n\n", funcname, *errflag);
      return(1);
    }
  }
  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\n MEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return (1);
  }

  return (0); // so returns 0 if everything is okay

}


static int MassBalances (realtype t, N_Vector y, N_Vector ydot, void *user_data){

  // realtype t1 = t;

  // N_Vector ydot = NULL;
  // ydot = N_VNew_Serial(3);

  NV_Ith_S(ydot,0) = -0.04 * NV_Ith_S(y,0) + 1e04 * NV_Ith_S(y,1) * NV_Ith_S(y,2);
  NV_Ith_S(ydot,2) = 3e07 * NV_Ith_S(y,1) * NV_Ith_S(y,1);
  NV_Ith_S(ydot,1) = -NV_Ith_S(ydot,0) - NV_Ith_S(ydot,2);

  // ydot[0] = NV_Ith_S(dydt,0);
  // ydot[1] = NV_Ith_S(dydt,1);
  // ydot[2] = NV_Ith_S(dydt,2);

  return(0);

}



// [[Rcpp::export]]
int Rcppcvode (double t, NumericVector IC, NumericVector ydot, double reltolerance, NumericVector abstolerance){

  int flag;
  realtype reltol = reltolerance;
  N_Vector abstol, y0;
  realtype T0 = RCONST(0.0);
  realtype tout = T1;
  realtype iout;


  // Set the vector absolute tolerance
  abstol = N_VNew_Serial(3);
  for (int i = 0; i<3; i++){
    NV_Ith_S(abstol, i) = abstolerance[i];
  }

  // Set the initial conditions
  y0 = N_VNew_Serial(3);
  // for (int i = 0; i<3; i++){
    NV_Ith_S(y0, 0) = 1.0; //IC[i];
    NV_Ith_S(y0, 1) = 0.0;
    NV_Ith_S(y0, 2) = 0.0;
  // }

  // void pointer, can point to any data type -  will have to be cast before use
  void *cvode_mem;
  cvode_mem = NULL;

  // Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula
  // and the use of a Newton iteration
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *) cvode_mem, "CVodeCreate", 0)) {
    return (1);
  }

  flag = CVodeInit(cvode_mem, MassBalances, T0, y0);
  if (check_flag(&flag, "CVodeInit", 1)) {
    return (1);
  }

  // Call CVodeSVtolerances to specify the scalar relative tolerance and vector absolute tol
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_flag(&flag, "CVodeSVtolerances", 1)) {
    return (1);
  }

  // Call CVDense to specify the CVDENSE dense linear solver
  flag = CVDense(cvode_mem, 3);
  if (check_flag(&flag, "CVDense", 1)) {
    return (1);
  }

  // Call CVodeInit to initialize the integrator memory and specify the user's right hand side function in y'=f(t,y),
  // the inital time T0, and the initial dependent variable vector y.
  iout = 0;
  while (1) {

    flag = CVode(cvode_mem, tout, y0, &t, CV_NORMAL);

    if (check_flag(&flag, "CVode", 1)) { break; } // Something went wrong in solving it!
    if (flag == CV_SUCCESS) {

      Rcout << t <<  "\t" << NV_Ith_S(y0,0) << "\t" << NV_Ith_S(y0,1) << "\t" << NV_Ith_S(y0,2) << '\n';
      iout++;
      tout *= TMULT;
    }

    if (iout > NOUT) { break; }
  }

  return(0); // everything went fine
}


typedef int (*funcPtr)(realtype t, N_Vector y, N_Vector ydot, void *user_data);

// [[Rcpp::export]]
int Rcppcvode_str (double t, NumericVector IC, SEXP xpsexp, double reltolerance, NumericVector abstolerance){

  int flag;
  realtype reltol = reltolerance;
  N_Vector abstol, y0;
  realtype T0 = RCONST(0.0);
  realtype tout = T1;
  realtype iout;


  // Set the vector absolute tolerance
  abstol = N_VNew_Serial(3);
  for (int i = 0; i<3; i++){
    NV_Ith_S(abstol, i) = abstolerance[i];
  }

  // Set the initial conditions
  y0 = N_VNew_Serial(3);
  // for (int i = 0; i<3; i++){
  NV_Ith_S(y0, 0) = 1.0; //IC[i];
  NV_Ith_S(y0, 1) = 0.0;
  NV_Ith_S(y0, 2) = 0.0;
  // }

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

  // Call CVDense to specify the CVDENSE dense linear solver
  flag = CVDense(cvode_mem, 3);
  if (check_flag(&flag, "CVDense", 1)) {
    return (1);
  }

  // Call CVodeInit to initialize the integrator memory and specify the user's right hand side function in y'=f(t,y),
  // the inital time T0, and the initial dependent variable vector y.
  iout = 0;
  while (1) {

    flag = CVode(cvode_mem, tout, y0, &t, CV_NORMAL);

    if (check_flag(&flag, "CVode", 1)) { break; } // Something went wrong in solving it!
    if (flag == CV_SUCCESS) {

      Rcout << t <<  "\t" << NV_Ith_S(y0,0) << "\t" << NV_Ith_S(y0,1) << "\t" << NV_Ith_S(y0,2) << '\n';
      iout++;
      tout *= TMULT;
    }

    if (iout > NOUT) { break; }
  }

  return(0); // everything went fine
}
