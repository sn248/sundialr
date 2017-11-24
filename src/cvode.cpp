
#include <Rcpp.h>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>    /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>         /* prototype for CVDense */
#include <sundials/sundials_dense.h>   /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h>   /* definition of type realtype */

#include "checkflag.h"

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
int cvode (NumericVector time_vec, NumericVector IC, SEXP xpsexp,
           double reltolerance, NumericVector abstolerance){

  int flag;
  realtype reltol = reltolerance;
  N_Vector abstol, y0;
  // realtype T0 = RCONST(0.0);
  // realtype tout = T1;
  realtype T0 = RCONST(time_vec[0]); //RCONST(0.0);  // Initial Time
  realtype tout = RCONST(time_vec[1]); // T1;        // First output time
  realtype iout;

  double time;
  int NOUT = time_vec.length() - 1;

  // Set the vector absolute tolerance
  abstol = N_VNew_Serial(abstolerance.length());
  for (int i = 0; i<abstolerance.length(); i++){
    NV_Ith_S(abstol, i) = abstolerance[i];
  }

  // Set the initial conditions
  y0 = N_VNew_Serial(IC.length());
  for (int i = 0; i<IC.length(); i++){
    NV_Ith_S(y0, i) = IC[i];
  }

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

  // Call CVodeInit to initialize the integrator memory and specify the user's right hand side function in y'=f(time,y),
  // the inital time T0, and the initial dependent variable vector y.
  iout = 0;
  while (1) {

    flag = CVode(cvode_mem, tout, y0, &time, CV_NORMAL);

    if (check_flag(&flag, "CVode", 1)) { break; } // Something went wrong in solving it!
    if (flag == CV_SUCCESS) {

      Rcout << time <<  "\t" << NV_Ith_S(y0,0) << "\t" << NV_Ith_S(y0,1) << "\t" << NV_Ith_S(y0,2) << '\n';
      iout++;
      tout += TMULT; // tout *= TMULT orginally
    }

    if (iout > NOUT) { break; }
  }

  return(0); // everything went fine
}
