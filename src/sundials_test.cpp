#include <Rcpp.h>

// #include <iostream.h>
#include "../inst/include/sundials/cvode/cvode.h"               /* prototypes for CVODE fcts., consts. */
#include "../inst/include/sundials/nvector/nvector_serial.h"    /* serial N_Vector types, fcts., macros */
#include "../inst/include/sundials/cvode/cvode_dense.h"         /* prototype for CVDense */
#include "../inst/include/sundials/sundials/sundials_dense.h"   /* definitions DlsMat DENSE_ELEM */
#include "../inst/include/sundials/sundials/sundials_types.h"   /* definition of type realtype */

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

const int NEQ = 3;
N_Vector y = NULL;
N_Vector abstol = NULL;

realtype reltol, t, tout;
void *cvode_mem = NULL;     // void pointer, can point to any data type -  will have to be cast before use
int flag, iout;

y = N_VNew_Serial(NEQ);
// if (check_flag((void *) y, "N_VNew_Serial", 0)) {
//   return (1);
// }
//
// abstol = N_VNew_Serial(NEQ);
// if (check_flag((void *) abstol, "N_VNew_Serial", 0)) {
//   return (1);
// }

// [[Rcpp::export]]
NumericVector IC(NumericVector x) {
  return x;
}

// [[Rcpp::export]]
NumericVector ode_fun(NumericVector x) {

  NumericVector dxdt(x.length(), 0.0);  //initialize vector of length x with zeros

  dxdt[0] = -0.04 * x[0] + 1.0e4 * x[1] * x[2];
  dxdt[2] = 3.0e7* x[1] * x[1];
  dxdt[1] = -dxdt[0] - dxdt[2];

  return dxdt;
}



// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
