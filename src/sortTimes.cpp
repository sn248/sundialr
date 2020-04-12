/* This function sorts the dosing and sampling times
 * in the increasing order of sampling times with
 * dosing records inserted in between in the order of
 * increasing times, e.g.,
 * > TDOSE - first column is time, second column is 1
 * indicating dose times
 *      [,1] [,2]
 * [1,]    5    1
 * [2,]    2    1
 * [3,]    9    1
 * [4,]    7    1
 * [5,]    1    1
 * > TSAMP - first column is time, second column is 0
 * indicating sampling times
 *       [,1] [,2]
 * [1,]    2    0
 * [2,]    6    0
 * [3,]    4    0
 * [4,]   10    0
 * [5,]    8    0
 * The output is
 * sorted_mat(TDOSE,TSAMP)
 *       [,1] [,2]
 * [2,]    2    1
 * [3,]    2    0
 * [4,]    4    0
 * [1,]    1    1
 * [5,]    5    1
 * [6,]    6    0
 * [7,]    7    1
 * [8,]    8    0
 * [9,]    9    1
 * [10,]   10   0
 */
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

////[[Rcpp::export]] - not exported
Rcpp::NumericMatrix sorted_mat(Rcpp::NumericMatrix TDOSE1, Rcpp::NumericMatrix TSAMP1){

  arma::mat TDOSE(TDOSE1.begin(), TDOSE1.nrow(), TDOSE1.ncol(), false);
  arma::mat TSAMP(TSAMP1.begin(), TSAMP1.nrow(), TSAMP1.ncol(), false);

  arma::mat TCOMB = arma::join_vert(TDOSE, TSAMP);   // rbind in R
  arma::mat TOUT = arma::zeros(TCOMB.n_rows,2);
  // TOUT.col(1) = -1;

  if(sum(TSAMP.col(1)) != 0){
    Rcpp::stop("Something wrong in Sampling Time Matrix \nAll elements in second column should be 0!");
  }

  if(sum(TDOSE.col(1)) != TDOSE.n_rows){
    Rcpp::stop("Something is wrong in the Dosing Matrix \nAll elements in second column should be 1");
  }

  // sort based on the first column (TIME)
  arma::uvec sorted_index = arma::stable_sort_index(TCOMB.col(0));

  // rearrange TOUT based on sorted_index
  for (unsigned int index = 0; index < TOUT.n_rows; index++){

    TOUT(index,0) = TCOMB(sorted_index(index),0);
    TOUT(index,1) = TCOMB(sorted_index(index),1);
  }
  // return TOUT; // sorted_index;

  Rcpp::NumericMatrix TOUT1(TOUT.n_rows, TOUT.n_cols, TOUT.begin());
  return TOUT1;

}

/*** R
TDOSE <- matrix(c(5,2,9,7,1,1,1,1,1,1), ncol = 2)
TSAMP <- matrix(c(2,6,4,10,8,0,0,0,0,0), ncol = 2)
TDOSE
TSAMP
sorted_mat(TDOSE,TSAMP)
*/


