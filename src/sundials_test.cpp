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

// [[Rcpp::export]]
NumericVector MassBalances (double t, NumericVector y, NumericVector ydot){

  realtype t1 = t;

  N_Vector dydt = NULL;

  dydt = N_VNew_Serial(3);

  NV_Ith_S(dydt,0) = -0.04 * y[0] + 1e04 * y[1] * y[2];
  NV_Ith_S(dydt,2) = 3e07 * y[1] *y[2];
  NV_Ith_S(dydt,1) = -NV_Ith_S(dydt,0) - NV_Ith_S(dydt,2);

  ydot[0] = NV_Ith_S(dydt,0);
  ydot[1] = NV_Ith_S(dydt,1);
  ydot[2] = NV_Ith_S(dydt,2);

  return ydot;

}

// [[Rcpp::export]]
int Rcppcvode (double t, NumericVector IC, NumericVector ydot, double reltolerance, double abstolerance){

  int flag;
  realtype reltol = reltolerance;
  N_Vector abstol, y;


  // Set the vector absolute tolerance
  abstol = N_VNew_Serial(3);
  for (int i = 0; i<3; i++){
    NV_Ith_S(abstol, 0) = abstolerance;
  }

  // Set the initial conditions
  y = N_VNew_Serial(3);
  for (int i = 0; i<3; i++){
    NV_Ith_S(y, 0) = IC[i];
  }

  // void pointer, can point to any data type -  will have to be cast before use
  void *cvode_mem;
  cvode_mem = NULL;

  // Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula and the use of a Newton iteration
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *) cvode_mem, "CVodeCreate", 0)) {
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



  return(0); // everything went fine
}

// [[Rcpp::export]]
Rcpp::NumericVector timesTwo(Rcpp::NumericVector x) {
  return x * 2;
}
