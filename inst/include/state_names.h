// File: state_names.h

#ifndef STATE_NAMES_H
#define STATE_NAMES_H

// Prerequisites: Rcpp.h

// Every callback copies the SUNDIALS state into a fresh NumericVector before
// handing it to the user's R function, so the names on the IC the caller
// supplied would otherwise never reach that function: y arrives unnamed and a
// right-hand side written by name (y["A"]) cannot work without an adapter.
// Each solver keeps names(IC) in its callback data and stamps it onto every
// state-shaped vector it builds (y, ydot, yS) with this helper. names is
// R_NilValue when IC is unnamed, and then nothing changes.
//
// The parameter vector needs no such treatment: the callbacks pass the
// caller's own Parameters object (or, in cvodes with IC sensitivities, a copy
// that carries the same names), so its names already arrive intact.
static inline void set_state_names(Rcpp::NumericVector &v, SEXP names) {
  if (names != R_NilValue) v.attr("names") = names;
}

#endif /* STATE_NAMES_H */
