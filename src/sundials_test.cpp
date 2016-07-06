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

// [[Rcpp::export]]
void InitialConditions (NumericVector x){

  N_Vector y = NULL;
  // y = N_VNew_Serial(3);

  // NV_Ith_S(y,0) = 1;
  // NV_Ith_S(y,1) = 2;
  // NV_Ith_S(y,2) = 3;


}

// [[Rcpp::export]]
Rcpp::NumericVector timesTwo(Rcpp::NumericVector x) {
  return x * 2;
}
