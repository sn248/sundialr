#include <Rcpp.h>

#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <cvodes/cvodes_direct.h>      /* access to CVDls interface            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* definition of ABS */

using namespace Rcpp;



//------------------------------------------------------------------------------
//' cvodes
//'
//' CVODES solver to solve ODEs and calculate sensitivities
//' @param time_vector time vector
//' @param IC Initial Conditions

