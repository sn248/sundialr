#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

Rcpp::NumericMatrix sorted_mat(Rcpp::NumericMatrix TDOSE1, Rcpp::NumericMatrix TSAMP1){

  arma::mat TDOSE(TDOSE1.begin(), TDOSE1.nrow(), TDOSE1.ncol(), false);
  arma::mat TSAMP(TSAMP1.begin(), TSAMP1.nrow(), TSAMP1.ncol(), false);

  arma::mat TCOMB = arma::join_vert(TDOSE, TSAMP);   // rbind in R
  arma::mat TOUT = arma::zeros(TCOMB.n_rows,2);
  // TOUT.col(1) = -1;

  if(sum(TSAMP.col(1)) != 0){
    Rcpp::stop("Something wrong in T Sample!");
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
sorted_mat(TDOSE,TSAMP)
*/

