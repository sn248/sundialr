// File: sortTimes.h
//
// Declares sorted_times(), defined in sortTimes.cpp and used by cvsolve().
// Compiled as its own translation unit: including the .cpp directly would
// define the function again in any other file that included it.

#ifndef SORT_TIMES_H
#define SORT_TIMES_H

#include <Rcpp.h>

// Merges the dosing records in TDOSE with the sampling times in TSAMP into one
// matrix ordered by time, with a fourth column flagging dosing rows.
Rcpp::NumericMatrix sorted_times(Rcpp::DataFrame TDOSE,
                                 Rcpp::NumericVector TSAMP,
                                 int NSTATES);

#endif /* SORT_TIMES_H */
