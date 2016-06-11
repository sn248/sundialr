#include <Rcpp.h>

// #include <iostream.h>
#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>    /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>         /* prototype for CVDense */
#include <sundials/sundials_dense.h>   /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h>   /* definition of type realtype */


// [[Rcpp::export]]
Rcpp::NumericVector timesTwo(Rcpp::NumericVector x) {
  return x * 2;
}