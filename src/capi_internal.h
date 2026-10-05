//   Copyright (c) 2016-2026, Satyaprakash Nayak
//   Distributed under the BSD-3 licence; see the header of sundialr_capi.h.

// Pieces shared by the CVODE and ARKODE halves of the C API
// (sundialr_capi.cpp, sundialr_capi_arkode.cpp). Private to the package: it is
// in src/, not inst/include/, so a consumer never sees it, and it deliberately
// includes neither Rcpp nor the R API, since no C API entry point may touch R.

#ifndef SUNDIALR_CAPI_INTERNAL_H
#define SUNDIALR_CAPI_INTERNAL_H

#include <cmath>
#include <exception>
#include <string>

#include <nvector/nvector_serial.h>

// Copy the current failure into last_err_msg: the SUNDIALS-recorded message if
// there is one, otherwise a description naming the solver and numeric code.
// H is a handle type with `err` (sundials_err_record) and `last_err_msg`.
template <class H>
static void capi_capture_err(H* h, int flag, const char* solver) {
  if (h->err.has_error) h->last_err_msg = h->err.message;
  else h->last_err_msg = std::string(solver) + " returned error code " + std::to_string(flag);
}

// Clear the per-call error state before an operation that can record one.
template <class H>
static void capi_clear_err(H* h) {
  h->err.has_error = false;
  h->err.message.clear();
}

// Error weights realising per-equation rtol/atol: w[i] = 1/(rtol[i]*|y[i]| +
// atol[i]). Used for CVODE's itol == 4 and ARKODE's equivalent, neither of which
// has a vector-rtol setter. H has `neq`, `rtol_v` and `atol_v`.
template <class H>
static int capi_fill_ewt(const H* h, N_Vector y, N_Vector w) {
  int n = h->neq;
  double* yp = N_VGetArrayPointer(y);
  double* wp = N_VGetArrayPointer(w);
  for (int i = 0; i < n; i++) {
    double ww = h->rtol_v[i] * std::fabs(yp[i]) + h->atol_v[i];
    if (ww <= 0.0) return -1;
    wp[i] = 1.0 / ww;
  }
  return 0;
}

// Wrap every entry-point body so no C++ exception can escape to the caller (the
// whole reason the C API exists in this form). A thrown std::exception is
// recorded and turned into an error code; anything else becomes a generic
// memory error.
#define CAPI_GUARD(handle, failcode, body)                                     \
  try { body }                                                                 \
  catch (std::exception& e) {                                                  \
    if (handle) { (handle)->last_err_msg = e.what(); }                         \
    return failcode;                                                           \
  }                                                                            \
  catch (...) {                                                                \
    if (handle) { (handle)->last_err_msg = "unidentified C++ exception in sundialr C API"; } \
    return failcode;                                                           \
  }

#endif  // SUNDIALR_CAPI_INTERNAL_H
