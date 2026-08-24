//   Copyright (c) 2016-2026, Satyaprakash Nayak
//   Distributed under the BSD-3 licence; see the header of sundialr_capi.h.

// Test drivers for the sundialr C API. Each one exercises the plain-C surface
// exactly as an external consumer would - create / configure /
// reinit / solve / free - and returns results to R so testthat can compare them
// with closed-form solutions (see tests/testthat/test-capi.r). They call the C
// API directly because they are linked into the same package.
//
// The R wrappers are named ".capi_test_*": the leading dot keeps them out of the
// package namespace (exportPattern only exports alpha-initial names), so they
// ship without being user-visible and need no documentation for R CMD check.

#include <Rcpp.h>
#include <string>
#include <thread>
#include <vector>
#include <sundialr_capi.h>

// White-box hook from sundialr_capi.cpp (deliberately not in the public
// header): the address of the handle's state buffer, for the tight-loop
// no-reallocation check below.
extern "C" const void* sundialr_cvode_state_ptr_internal(void* m);

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
// (state continuity across reinit), the reinit-per-segment pattern a dose/event-
// driven host loop uses.
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

// One handle serving two "subjects": solve decay with k1, repoint the callback
// data at k2 with set_udata, reinit and solve again. Column j holds subject
// j's trajectory; each must match its own closed form, which fails if the old
// udata pointer is still being read after the switch.
// [[Rcpp::export(".capi_test_set_udata")]]
NumericMatrix capi_test_set_udata(NumericVector times, double y0,
                                  double k1, double k2) {
  int n = times.size();
  double ka = k1, kb = k2;
  void* m = sundialr_cvode_create(1, &ka);
  if (!m) stop("sundialr_cvode_create returned NULL");

  sundialr_cvode_set_rhs(m, decay_rhs);
  sundialr_cvode_set_tol_scalar(m, 1e-10, 1e-12);
  sundialr_cvode_set_max_steps(m, 100000);

  NumericMatrix out(n, 2);
  for (int subj = 0; subj < 2; subj++) {
    if (subj == 1) {
      if (sundialr_cvode_set_udata(m, &kb) != 0) capi_test_fail(m, "set_udata failed");
    }
    double ycur = y0;
    out(0, subj) = ycur;
    if (sundialr_cvode_reinit(m, times[0], &ycur) < 0) capi_test_fail(m, "reinit failed");
    for (int i = 1; i < n; i++) {
      double yo = 0.0, tr = 0.0;
      if (sundialr_cvode_solve(m, times[i], &yo, &tr) < 0) capi_test_fail(m, "solve failed");
      out(i, subj) = yo;
    }
  }

  sundialr_cvode_free(m);
  return out;
}

// The handle-reuse guarantee: nseg reinit+solve segments on one handle, the
// dose/event cadence of a host fitting loop at realistic volume. Checks that
// the state buffer's address never changes (direct evidence the state vector
// is not reallocated; CVodeReInit's no-malloc contract covers the rest), that
// the trajectory stays exact, and that get_num_steps accumulates across the
// reinits instead of resetting with CVODE's internal counter.
// [[Rcpp::export(".capi_test_tight_loop")]]
List capi_test_tight_loop(int nseg, double dt, double y0, double k) {
  double kk = k;
  void* m = sundialr_cvode_create(1, &kk);
  if (!m) stop("sundialr_cvode_create returned NULL");

  sundialr_cvode_set_rhs(m, decay_rhs);
  sundialr_cvode_set_tol_scalar(m, 1e-10, 1e-12);
  sundialr_cvode_set_max_steps(m, 100000);

  double t = 0.0, ycur = y0;
  if (sundialr_cvode_reinit(m, t, &ycur) < 0) capi_test_fail(m, "reinit failed");
  const void* p0 = sundialr_cvode_state_ptr_internal(m);
  bool ptr_stable = (p0 != NULL);

  long steps_mid = -1;
  for (int i = 0; i < nseg; i++) {
    double yo = 0.0, tr = 0.0;
    if (sundialr_cvode_solve(m, t + dt, &yo, &tr) < 0) capi_test_fail(m, "solve failed");
    t += dt;
    ycur = yo;
    if (sundialr_cvode_reinit(m, t, &ycur) < 0) capi_test_fail(m, "reinit failed");
    if (sundialr_cvode_state_ptr_internal(m) != p0) ptr_stable = false;
    if (i == nseg / 2) steps_mid = sundialr_cvode_get_num_steps(m);
  }
  long steps_end = sundialr_cvode_get_num_steps(m);

  sundialr_cvode_free(m);
  return List::create(_["y_end"]      = ycur,
                      _["t_end"]      = t,
                      _["ptr_stable"] = ptr_stable,
                      _["steps_mid"]  = (double) steps_mid,
                      _["steps_end"]  = (double) steps_end);
}

// reset_stats semantics: the count reads zero immediately after a reset (even
// mid-segment), then resumes accumulating, including across a later reinit.
// [[Rcpp::export(".capi_test_reset_stats")]]
List capi_test_reset_stats() {
  double kk = 0.6;
  void* m = sundialr_cvode_create(1, &kk);
  if (!m) stop("sundialr_cvode_create returned NULL");

  sundialr_cvode_set_rhs(m, decay_rhs);
  sundialr_cvode_set_tol_scalar(m, 1e-10, 1e-12);
  sundialr_cvode_set_max_steps(m, 100000);

  double ycur = 3.0, yo = 0.0, tr = 0.0;
  if (sundialr_cvode_reinit(m, 0.0, &ycur) < 0) capi_test_fail(m, "reinit failed");
  if (sundialr_cvode_solve(m, 1.0, &yo, &tr) < 0) capi_test_fail(m, "solve failed");
  long before = sundialr_cvode_get_num_steps(m);

  if (sundialr_cvode_reset_stats(m) != 0) capi_test_fail(m, "reset_stats failed");
  long at_reset = sundialr_cvode_get_num_steps(m);   // mid-segment: must be 0

  if (sundialr_cvode_solve(m, 2.0, &yo, &tr) < 0) capi_test_fail(m, "solve failed");
  long after_solve = sundialr_cvode_get_num_steps(m);

  if (sundialr_cvode_reinit(m, 2.0, &yo) < 0) capi_test_fail(m, "reinit failed");
  if (sundialr_cvode_solve(m, 3.0, &yo, &tr) < 0) capi_test_fail(m, "solve failed");
  long after_reinit = sundialr_cvode_get_num_steps(m);

  sundialr_cvode_free(m);
  return List::create(_["before"]       = (double) before,
                      _["at_reset"]     = (double) at_reset,
                      _["after_solve"]  = (double) after_solve,
                      _["after_reinit"] = (double) after_reinit);
}

// --- Concurrency: two handles driven from two std::threads -----------------
// Everything a thread touches is preallocated here and owned by its job; the
// thread bodies call ONLY the C API (which never throws and never calls the R
// API), so no R interaction happens off the main thread. Two threads, matching
// CRAN's limit on cores used by tests.

namespace {

struct capi_thread_job {
  double k = 0.0;
  double y0 = 0.0;
  std::vector<double> times;
  std::vector<double> out;
  int status = 0;                     // 0 = ok, < 0 = first failing code
};

void capi_thread_run(capi_thread_job* job) {
  void* m = sundialr_cvode_create(1, &job->k);
  if (!m) { job->status = -100; return; }
  sundialr_cvode_set_rhs(m, decay_rhs);
  sundialr_cvode_set_tol_scalar(m, 1e-10, 1e-12);
  sundialr_cvode_set_max_steps(m, 100000);

  double ycur = job->y0;
  job->out[0] = ycur;
  int flag = sundialr_cvode_reinit(m, job->times[0], &ycur);
  if (flag < 0) { job->status = flag; sundialr_cvode_free(m); return; }

  for (size_t i = 1; i < job->times.size(); i++) {
    double yo = 0.0, tr = 0.0;
    flag = sundialr_cvode_solve(m, job->times[i], &yo, &tr);
    if (flag < 0) { job->status = flag; break; }
    job->out[i] = yo;
    // Reinit every segment so the reinit path runs concurrently too.
    flag = sundialr_cvode_reinit(m, job->times[i], &yo);
    if (flag < 0) { job->status = flag; break; }
  }
  sundialr_cvode_free(m);
}

}  // namespace

// [[Rcpp::export(".capi_test_concurrent")]]
List capi_test_concurrent(NumericVector times, double y0, double k1, double k2) {
  int n = times.size();
  capi_thread_job j1, j2;
  j1.k = k1; j1.y0 = y0; j1.times.assign(times.begin(), times.end()); j1.out.assign(n, 0.0);
  j2.k = k2; j2.y0 = y0; j2.times.assign(times.begin(), times.end()); j2.out.assign(n, 0.0);

  std::thread t1(capi_thread_run, &j1);
  std::thread t2(capi_thread_run, &j2);
  t1.join();
  t2.join();

  if (j1.status < 0) stop("thread 1 failed with code " + std::to_string(j1.status));
  if (j2.status < 0) stop("thread 2 failed with code " + std::to_string(j2.status));

  return List::create(_["y1"] = NumericVector(j1.out.begin(), j1.out.end()),
                      _["y2"] = NumericVector(j2.out.begin(), j2.out.end()));
}
