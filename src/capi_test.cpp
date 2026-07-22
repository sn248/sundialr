//   Copyright (c) 2016-2026, Satyaprakash Nayak
//   Distributed under the BSD-3 licence; see the header of sundialr_capi.h.

// Test drivers for the sundialr C API. Each one exercises the plain-C surface
// exactly as an external consumer (mrgsolve) would - create / configure /
// reinit / solve / free - and returns results to R so testthat can compare them
// with closed-form solutions (see tests/testthat/test-capi.r). They call the C
// API directly because they are linked into the same package.
//
// The R wrappers are named ".capi_test_*": the leading dot keeps them out of the
// package namespace (exportPattern only exports alpha-initial names), so they
// ship without being user-visible and need no documentation for R CMD check.

#include <Rcpp.h>
#include <string>
#include <sundialr_capi.h>

using namespace Rcpp;

// --- Model callbacks (C linkage, so their type matches the C API typedefs) ---
extern "C" {

// Scalar exponential decay:  y' = -k y   ->  y(t) = y0 exp(-k t)
static int decay_rhs(double t, const double* y, double* ydot, void* udata) {
  (void) t;
  double k = *((double*) udata);
  ydot[0] = -k * y[0];
  return 0;
}

// Two-compartment linear:  y1' = -k1 y1 ;  y2' = k1 y1 - k2 y2
// With y1(0)=A, y2(0)=0:  y1 = A e^{-k1 t},
//   y2 = A k1/(k2-k1) (e^{-k1 t} - e^{-k2 t}).
static int twocmt_rhs(double t, const double* y, double* ydot, void* udata) {
  (void) t;
  double* p = (double*) udata;
  ydot[0] = -p[0] * y[0];
  ydot[1] =  p[0] * y[0] - p[1] * y[1];
  return 0;
}

// Column-major Jacobian: J[i + j*2] = d(ydot_i)/d(y_j).
static int twocmt_jac(double t, const double* y, double* J, void* udata) {
  (void) t; (void) y;
  double* p = (double*) udata;
  J[0] = -p[0];  // (0,0)
  J[1] =  p[0];  // (1,0)
  J[2] =  0.0;   // (0,1)
  J[3] = -p[1];  // (1,1)
  return 0;
}

}  // extern "C"

// Fetch the recorded message (or a fallback) and free the handle before raising.
static void capi_test_fail(void* m, const char* what) {
  const char* e = sundialr_cvode_last_err(m);
  std::string msg = e ? std::string(e) : std::string(what);
  sundialr_cvode_free(m);
  stop(msg);
}

// Scalar decay. tol_mode: 1 = scalar tolerances, 2 = per-equation vectors.
// reinit_each = true reinitialises at every output time with the current state
// (state continuity across reinit), the reinit-per-segment pattern mrgsolve uses.
// [[Rcpp::export(".capi_test_decay")]]
NumericVector capi_test_decay(NumericVector times, double y0, double k,
                              int tol_mode, bool reinit_each) {
  int n = times.size();
  double kk = k;
  void* m = sundialr_cvode_create(1, &kk);
  if (!m) stop("sundialr_cvode_create returned NULL");

  sundialr_cvode_set_rhs(m, decay_rhs);
  sundialr_cvode_set_max_steps(m, 100000);
  if (tol_mode == 2) {
    double rt = 1e-10, at = 1e-12;
    sundialr_cvode_set_tol_vector(m, &rt, &at);
  } else {
    sundialr_cvode_set_tol_scalar(m, 1e-10, 1e-12);
  }

  NumericVector out(n);
  out[0] = y0;
  double ycur = y0;
  if (sundialr_cvode_reinit(m, times[0], &ycur) < 0) capi_test_fail(m, "reinit failed");

  for (int i = 1; i < n; i++) {
    double yo = 0.0, tr = 0.0;
    if (sundialr_cvode_solve(m, times[i], &yo, &tr) < 0) capi_test_fail(m, "solve failed");
    out[i] = yo;
    if (reinit_each) {
      if (sundialr_cvode_reinit(m, times[i], &yo) < 0) capi_test_fail(m, "reinit failed");
    }
  }

  sundialr_cvode_free(m);
  return out;
}

// Two-compartment system; use_jac toggles the analytic Jacobian.
// [[Rcpp::export(".capi_test_twocmt")]]
NumericMatrix capi_test_twocmt(NumericVector times, double y10,
                               double k1, double k2, bool use_jac) {
  int n = times.size();
  double p[2] = {k1, k2};
  void* m = sundialr_cvode_create(2, p);
  if (!m) stop("sundialr_cvode_create returned NULL");

  sundialr_cvode_set_rhs(m, twocmt_rhs);
  if (use_jac) sundialr_cvode_set_jac(m, twocmt_jac);
  sundialr_cvode_set_tol_scalar(m, 1e-10, 1e-12);
  sundialr_cvode_set_max_steps(m, 100000);

  NumericMatrix out(n, 2);
  double y[2] = {y10, 0.0};
  if (sundialr_cvode_reinit(m, times[0], y) < 0) capi_test_fail(m, "reinit failed");
  out(0, 0) = y[0];
  out(0, 1) = y[1];

  for (int i = 1; i < n; i++) {
    double yo[2] = {0.0, 0.0}, tr = 0.0;
    if (sundialr_cvode_solve(m, times[i], yo, &tr) < 0) capi_test_fail(m, "solve failed");
    out(i, 0) = yo[0];
    out(i, 1) = yo[1];
  }

  sundialr_cvode_free(m);
  return out;
}

// Force a failure (maxsteps = 1 over a long interval) and report the code and
// message, to confirm the C API returns CVODE's own negative code and records a
// non-empty message rather than throwing.
// [[Rcpp::export(".capi_test_force_error")]]
List capi_test_force_error() {
  double kk = 1.0;
  void* m = sundialr_cvode_create(1, &kk);
  if (!m) stop("sundialr_cvode_create returned NULL");

  sundialr_cvode_set_rhs(m, decay_rhs);
  sundialr_cvode_set_tol_scalar(m, 1e-12, 1e-14);
  sundialr_cvode_set_max_steps(m, 1);          // far too few

  double y0 = 1.0;
  sundialr_cvode_reinit(m, 0.0, &y0);
  double yo = 0.0, tr = 0.0;
  int flag = sundialr_cvode_solve(m, 1e6, &yo, &tr);

  const char* e = sundialr_cvode_last_err(m);
  std::string es = e ? std::string(e) : std::string("");
  sundialr_cvode_free(m);

  return List::create(_["flag"] = flag, _["err"] = es);
}

// Number of internal steps for a decay solve, to confirm get_num_steps works.
// [[Rcpp::export(".capi_test_num_steps")]]
double capi_test_num_steps(double tout, double y0, double k) {
  double kk = k;
  void* m = sundialr_cvode_create(1, &kk);
  if (!m) stop("sundialr_cvode_create returned NULL");

  sundialr_cvode_set_rhs(m, decay_rhs);
  sundialr_cvode_set_tol_scalar(m, 1e-8, 1e-10);
  sundialr_cvode_set_max_steps(m, 100000);

  double yc = y0;
  sundialr_cvode_reinit(m, 0.0, &yc);
  double yo = 0.0, tr = 0.0;
  sundialr_cvode_solve(m, tout, &yo, &tr);
  long ns = sundialr_cvode_get_num_steps(m);

  sundialr_cvode_free(m);
  return (double) ns;
}

// After a successful solve, last_err must be NULL (returned as TRUE here).
// [[Rcpp::export(".capi_test_clean_err")]]
bool capi_test_clean_err() {
  double kk = 1.0;
  void* m = sundialr_cvode_create(1, &kk);
  if (!m) stop("sundialr_cvode_create returned NULL");

  sundialr_cvode_set_rhs(m, decay_rhs);
  sundialr_cvode_set_tol_scalar(m, 1e-8, 1e-10);
  sundialr_cvode_set_max_steps(m, 100000);

  double yc = 1.0;
  sundialr_cvode_reinit(m, 0.0, &yc);
  double yo = 0.0, tr = 0.0;
  sundialr_cvode_solve(m, 1.0, &yo, &tr);
  bool clean = (sundialr_cvode_last_err(m) == NULL);

  sundialr_cvode_free(m);
  return clean;
}

// [[Rcpp::export(".capi_test_abi")]]
int capi_test_abi() {
  return sundialr_abi_version();
}
