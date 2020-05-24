/* This function sorts the dosing dataframe and sampling times
 * in the increasing order of sampling times with
 * dosing records inserted in between in the order of
 * increasing times, e.g., dosing second species
 > TDOSE <- data.frame(ID = 2, TIMES = c(0,10, 20, 30), VAL = 100)
 > TDOSE
 ID TIMES VAL
 1  2     0 100
 2  2    10 100
 3  2    20 100
 4  2    30 100
 > TSAMP <- seq(from = 0, to = 50, by = 10)
 > TSAMP
 [1]  0 10 20 30 40 50
 */

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
// not exported
////[[Rcpp::export]]
Rcpp::NumericMatrix sorted_times(Rcpp::DataFrame TDOSE, Rcpp::NumericVector TSAMP){

  // Dosing dataframe - 1st column - index of species being dosed + 1
  // Dosing dataframe - 2nd column - time of dosing
  // Dosing dataframe - 3rd column - Value of dosing

  // convert Dosing dataframe to matrix first and then an arma matrix
  // add another column (4th column) of all ones to indicate dosing
  Rcpp::NumericMatrix TDOSE1 = Rcpp::internal::convert_using_rfunction(TDOSE, "as.matrix");
  arma::mat TDOSE2(TDOSE1.begin(), TDOSE1.nrow(), TDOSE1.ncol(), false);
  TDOSE2.insert_cols(3, 1);  // add 1 column at column position 4
  TDOSE2.col(3).ones();      // set column 4 to all ones T0 INDICATE DOSING

  // convert Sampling vector to an arma matrix
  // all columns are zero except 2nd column for time
  // 1st column - zero indicates no species in being dosed
  // 2nd column - value indicates sampling time
  // 3rd column - zero indicates dose is zero
  // 4th column - zero indicates sampling time
  Rcpp::NumericMatrix TSAMP1(TSAMP.length(), 4); // a matrix of zeros
  // second column contains the sampling times
  for(int i = 0; i < TSAMP1.nrow(); i++){
    TSAMP1(i,1) = TSAMP(i);      // fill second column with sampling times
  }

  arma::mat TSAMP2(TSAMP1.begin(), TSAMP1.nrow(), TSAMP1.ncol(), false);
  if(sum(TSAMP2.col(0)) != 0 & sum(TSAMP2.col(3)) != 0 ){
    Rcpp::stop("Something wrong in Sampling Time Matrix \nAll elements in first and third columns should be 0!");
  }

  arma::mat TCOMB = arma::join_vert(TDOSE2, TSAMP2);        // rbind in R
  arma::mat TOUT = arma::zeros(TCOMB.n_rows, TCOMB.n_cols); // output

  // sort based on the second column (TIME)
  arma::uvec sorted_index = arma::stable_sort_index(TCOMB.col(1));

  // rearrange TOUT based on sorted_index
  for (unsigned int index = 0; index < TOUT.n_rows; index++){
    TOUT(index,0) = TCOMB(sorted_index(index),0);
    TOUT(index,1) = TCOMB(sorted_index(index),1);
    TOUT(index,2) = TCOMB(sorted_index(index),2);
    TOUT(index,3) = TCOMB(sorted_index(index),3);
  }

  Rcpp::NumericMatrix TOUT1(TOUT.n_rows, TOUT.n_cols, TOUT.begin());
  return TOUT1;

}

// /*** R
// TDOSE <- DOSE <- data.frame(ID = 2, TIMES = c(0,10, 20, 30), VAL = 100)
// TSAMP <- seq(from = 0, to = 50, by = 10)
// TDOSE
// TSAMP
// sorted_mat(TDOSE,TSAMP)
// */


