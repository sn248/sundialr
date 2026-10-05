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

// White-box hooks from sundialr_capi.cpp and sundialr_capi_arkode.cpp
// (deliberately not in the public header): the address of a handle's state
// buffer, for the tight-loop no-reallocation check below.
extern "C" const void* sundialr_cvode_state_ptr_internal(void* m);
extern "C" const void* sundialr_arkode_state_ptr_internal(void* m);

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

// Robertson chemical kinetics, the standard stiff test problem:
//   y1' = -0.04 y1 + 1e4 y2 y3
//   y2' =  0.04 y1 - 1e4 y2 y3 - 3e7 y2^2
//   y3' =  3e7 y2^2
// Its rate constants span nine orders of magnitude; y1 + y2 + y3 stays 1.
static int robertson_rhs(double t, const double* y, double* ydot, void* udata) {
  (void) t; (void) udata;
  ydot[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
  ydot[2] =  3.0e7 * y[1] * y[1];
  ydot[1] = -ydot[0] - ydot[2];
  return 0;
}

static int robertson_jac(double t, const double* y, double* J, void* udata) {
  (void) t; (void) udata;
  // Column-major, J[i + 3*j] = d(ydot_i)/d(y_j).
  J[0] = -0.04;  J[3] =  1.0e4 * y[2];                  J[6] =  1.0e4 * y[1];
  J[2] =  0.0;   J[5] =  6.0e7 * y[1];                  J[8] =  0.0;
  J[1] =  0.04;  J[4] = -1.0e4 * y[2] - 6.0e7 * y[1];   J[7] = -1.0e4 * y[1];
  return 0;
}

}  // extern "C"

// --- The two integrators behind one table ----------------------------------
// The CVODE and ARKODE halves of the C API have identical signatures, so every
// driver below runs against either: solver = 0 selects CVODE, 1 ARKODE with its
// default ESDIRK.

struct capi_api {
  void* (*create)(int, void*);
  void  (*free)(void*);
  int   (*set_rhs)(void*, sundialr_rhs);
  int   (*set_jac)(void*, sundialr_jac);
  int   (*set_tol_scalar)(void*, double, double);
  int   (*set_tol_vector)(void*, const double*, const double*);
  int   (*set_max_steps)(void*, long);
  int   (*set_udata)(void*, void*);
  int   (*reinit)(void*, double, const double*);
  int   (*solve)(void*, double, double*, double*);
  long  (*get_num_steps)(void*);
  int   (*reset_stats)(void*);
  const char* (*last_err)(void*);
  const void* (*state_ptr_internal)(void*);
};

static const capi_api CVODE_API = {
  sundialr_cvode_create, sundialr_cvode_free, sundialr_cvode_set_rhs,
  sundialr_cvode_set_jac, sundialr_cvode_set_tol_scalar,
  sundialr_cvode_set_tol_vector, sundialr_cvode_set_max_steps,
  sundialr_cvode_set_udata, sundialr_cvode_reinit, sundialr_cvode_solve,
  sundialr_cvode_get_num_steps, sundialr_cvode_reset_stats,
  sundialr_cvode_last_err, sundialr_cvode_state_ptr_internal
};

static const capi_api ARKODE_API = {
  sundialr_arkode_create, sundialr_arkode_free, sundialr_arkode_set_rhs,
  sundialr_arkode_set_jac, sundialr_arkode_set_tol_scalar,
  sundialr_arkode_set_tol_vector, sundialr_arkode_set_max_steps,
  sundialr_arkode_set_udata, sundialr_arkode_reinit, sundialr_arkode_solve,
  sundialr_arkode_get_num_steps, sundialr_arkode_reset_stats,
  sundialr_arkode_last_err, sundialr_arkode_state_ptr_internal
};

static const capi_api& capi_select(int solver) {
  if (solver == 0) return CVODE_API;
  if (solver == 1) return ARKODE_API;
  stop("solver must be 0 (CVODE) or 1 (ARKODE)");
}

// Fetch the recorded message (or a fallback) and free the handle before raising.
static void capi_test_fail(const capi_api& A, void* m, const char* what) {
  const char* e = A.last_err(m);
  std::string msg = e ? std::string(e) : std::string(what);
  A.free(m);
  stop(msg);
}

// Scalar decay. tol_mode: 1 = scalar tolerances, 2 = per-equation vectors.
// reinit_each = true reinitialises at every output time with the current state
// (state continuity across reinit), the reinit-per-segment pattern a dose/event-
// driven host loop uses.
// [[Rcpp::export(".capi_test_decay")]]
NumericVector capi_test_decay(NumericVector times, double y0, double k,
                              int tol_mode, bool reinit_each, int solver = 0) {
  const capi_api& A = capi_select(solver);
  int n = times.size();
  double kk = k;
  void* m = A.create(1, &kk);
  if (!m) stop("sundialr_cvode_create returned NULL");

  A.set_rhs(m, decay_rhs);
  A.set_max_steps(m, 100000);
  if (tol_mode == 2) {
    double rt = 1e-10, at = 1e-12;
    A.set_tol_vector(m, &rt, &at);
  } else {
    A.set_tol_scalar(m, 1e-10, 1e-12);
  }

  NumericVector out(n);
  out[0] = y0;
  double ycur = y0;
  if (A.reinit(m, times[0], &ycur) < 0) capi_test_fail(A, m, "reinit failed");

  for (int i = 1; i < n; i++) {
    double yo = 0.0, tr = 0.0;
    if (A.solve(m, times[i], &yo, &tr) < 0) capi_test_fail(A, m, "solve failed");
    out[i] = yo;
    if (reinit_each) {
      if (A.reinit(m, times[i], &yo) < 0) capi_test_fail(A, m, "reinit failed");
    }
  }

  A.free(m);
  return out;
}

// Two-compartment system; use_jac toggles the analytic Jacobian.
// [[Rcpp::export(".capi_test_twocmt")]]
NumericMatrix capi_test_twocmt(NumericVector times, double y10,
                               double k1, double k2, bool use_jac, int solver = 0) {
  const capi_api& A = capi_select(solver);
  int n = times.size();
  double p[2] = {k1, k2};
  void* m = A.create(2, p);
  if (!m) stop("sundialr_cvode_create returned NULL");

  A.set_rhs(m, twocmt_rhs);
  if (use_jac) A.set_jac(m, twocmt_jac);
  A.set_tol_scalar(m, 1e-10, 1e-12);
  A.set_max_steps(m, 100000);

  NumericMatrix out(n, 2);
  double y[2] = {y10, 0.0};
  if (A.reinit(m, times[0], y) < 0) capi_test_fail(A, m, "reinit failed");
  out(0, 0) = y[0];
  out(0, 1) = y[1];

  for (int i = 1; i < n; i++) {
    double yo[2] = {0.0, 0.0}, tr = 0.0;
    if (A.solve(m, times[i], yo, &tr) < 0) capi_test_fail(A, m, "solve failed");
    out(i, 0) = yo[0];
    out(i, 1) = yo[1];
  }

  A.free(m);
  return out;
}

// Force a failure (maxsteps = 1 over a long interval) and report the code and
// message, to confirm the C API returns CVODE's own negative code and records a
// non-empty message rather than throwing.
// [[Rcpp::export(".capi_test_force_error")]]
List capi_test_force_error(int solver = 0) {
  const capi_api& A = capi_select(solver);
  double kk = 1.0;
  void* m = A.create(1, &kk);
  if (!m) stop("sundialr_cvode_create returned NULL");

  A.set_rhs(m, decay_rhs);
  A.set_tol_scalar(m, 1e-12, 1e-14);
  A.set_max_steps(m, 1);          // far too few

  double y0 = 1.0;
  A.reinit(m, 0.0, &y0);
  double yo = 0.0, tr = 0.0;
  int flag = A.solve(m, 1e6, &yo, &tr);

  const char* e = A.last_err(m);
  std::string es = e ? std::string(e) : std::string("");
  A.free(m);

  return List::create(_["flag"] = flag, _["err"] = es);
}

// Number of internal steps for a decay solve, to confirm get_num_steps works.
// [[Rcpp::export(".capi_test_num_steps")]]
double capi_test_num_steps(double tout, double y0, double k, int solver = 0) {
  const capi_api& A = capi_select(solver);
  double kk = k;
  void* m = A.create(1, &kk);
  if (!m) stop("sundialr_cvode_create returned NULL");

  A.set_rhs(m, decay_rhs);
  A.set_tol_scalar(m, 1e-8, 1e-10);
  A.set_max_steps(m, 100000);

  double yc = y0;
  A.reinit(m, 0.0, &yc);
  double yo = 0.0, tr = 0.0;
  A.solve(m, tout, &yo, &tr);
  long ns = A.get_num_steps(m);

  A.free(m);
  return (double) ns;
}

// After a successful solve, last_err must be NULL (returned as TRUE here).
// [[Rcpp::export(".capi_test_clean_err")]]
bool capi_test_clean_err(int solver = 0) {
  const capi_api& A = capi_select(solver);
  double kk = 1.0;
  void* m = A.create(1, &kk);
  if (!m) stop("sundialr_cvode_create returned NULL");

  A.set_rhs(m, decay_rhs);
  A.set_tol_scalar(m, 1e-8, 1e-10);
  A.set_max_steps(m, 100000);

  double yc = 1.0;
  A.reinit(m, 0.0, &yc);
  double yo = 0.0, tr = 0.0;
  A.solve(m, 1.0, &yo, &tr);
  bool clean = (A.last_err(m) == NULL);

  A.free(m);
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
                                  double k1, double k2, int solver = 0) {
  const capi_api& A = capi_select(solver);
  int n = times.size();
  double ka = k1, kb = k2;
  void* m = A.create(1, &ka);
  if (!m) stop("sundialr_cvode_create returned NULL");

  A.set_rhs(m, decay_rhs);
  A.set_tol_scalar(m, 1e-10, 1e-12);
  A.set_max_steps(m, 100000);

  NumericMatrix out(n, 2);
  for (int subj = 0; subj < 2; subj++) {
    if (subj == 1) {
      if (A.set_udata(m, &kb) != 0) capi_test_fail(A, m, "set_udata failed");
    }
    double ycur = y0;
    out(0, subj) = ycur;
    if (A.reinit(m, times[0], &ycur) < 0) capi_test_fail(A, m, "reinit failed");
    for (int i = 1; i < n; i++) {
      double yo = 0.0, tr = 0.0;
      if (A.solve(m, times[i], &yo, &tr) < 0) capi_test_fail(A, m, "solve failed");
      out(i, subj) = yo;
    }
  }

  A.free(m);
  return out;
}

// The handle-reuse guarantee: nseg reinit+solve segments on one handle, the
// dose/event cadence of a host fitting loop at realistic volume. Checks that
// the state buffer's address never changes (direct evidence the state vector
// is not reallocated; CVodeReInit's no-malloc contract covers the rest), that
// the trajectory stays exact, and that get_num_steps accumulates across the
// reinits instead of resetting with CVODE's internal counter.
// [[Rcpp::export(".capi_test_tight_loop")]]
List capi_test_tight_loop(int nseg, double dt, double y0, double k, int solver = 0) {
  const capi_api& A = capi_select(solver);
  double kk = k;
  void* m = A.create(1, &kk);
  if (!m) stop("sundialr_cvode_create returned NULL");

  A.set_rhs(m, decay_rhs);
  A.set_tol_scalar(m, 1e-10, 1e-12);
  A.set_max_steps(m, 100000);

  double t = 0.0, ycur = y0;
  if (A.reinit(m, t, &ycur) < 0) capi_test_fail(A, m, "reinit failed");
  const void* p0 = A.state_ptr_internal(m);
  bool ptr_stable = (p0 != NULL);

  long steps_mid = -1;
  for (int i = 0; i < nseg; i++) {
    double yo = 0.0, tr = 0.0;
    if (A.solve(m, t + dt, &yo, &tr) < 0) capi_test_fail(A, m, "solve failed");
    t += dt;
    ycur = yo;
    if (A.reinit(m, t, &ycur) < 0) capi_test_fail(A, m, "reinit failed");
    if (A.state_ptr_internal(m) != p0) ptr_stable = false;
    if (i == nseg / 2) steps_mid = A.get_num_steps(m);
  }
  long steps_end = A.get_num_steps(m);

  A.free(m);
  return List::create(_["y_end"]      = ycur,
                      _["t_end"]      = t,
                      _["ptr_stable"] = ptr_stable,
                      _["steps_mid"]  = (double) steps_mid,
                      _["steps_end"]  = (double) steps_end);
}

// reset_stats semantics: the count reads zero immediately after a reset (even
// mid-segment), then resumes accumulating, including across a later reinit.
// [[Rcpp::export(".capi_test_reset_stats")]]
List capi_test_reset_stats(int solver = 0) {
  const capi_api& A = capi_select(solver);
  double kk = 0.6;
  void* m = A.create(1, &kk);
  if (!m) stop("sundialr_cvode_create returned NULL");

  A.set_rhs(m, decay_rhs);
  A.set_tol_scalar(m, 1e-10, 1e-12);
  A.set_max_steps(m, 100000);

  double ycur = 3.0, yo = 0.0, tr = 0.0;
  if (A.reinit(m, 0.0, &ycur) < 0) capi_test_fail(A, m, "reinit failed");
  if (A.solve(m, 1.0, &yo, &tr) < 0) capi_test_fail(A, m, "solve failed");
  long before = A.get_num_steps(m);

  if (A.reset_stats(m) != 0) capi_test_fail(A, m, "reset_stats failed");
  long at_reset = A.get_num_steps(m);   // mid-segment: must be 0

  if (A.solve(m, 2.0, &yo, &tr) < 0) capi_test_fail(A, m, "solve failed");
  long after_solve = A.get_num_steps(m);

  if (A.reinit(m, 2.0, &yo) < 0) capi_test_fail(A, m, "reinit failed");
  if (A.solve(m, 3.0, &yo, &tr) < 0) capi_test_fail(A, m, "solve failed");
  long after_reinit = A.get_num_steps(m);

  A.free(m);
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
  const capi_api* api = NULL;
};

void capi_thread_run(capi_thread_job* job) {
  const capi_api& A = *job->api;
  void* m = A.create(1, &job->k);
  if (!m) { job->status = -100; return; }
  A.set_rhs(m, decay_rhs);
  A.set_tol_scalar(m, 1e-10, 1e-12);
  A.set_max_steps(m, 100000);

  double ycur = job->y0;
  job->out[0] = ycur;
  int flag = A.reinit(m, job->times[0], &ycur);
  if (flag < 0) { job->status = flag; A.free(m); return; }

  for (size_t i = 1; i < job->times.size(); i++) {
    double yo = 0.0, tr = 0.0;
    flag = A.solve(m, job->times[i], &yo, &tr);
    if (flag < 0) { job->status = flag; break; }
    job->out[i] = yo;
    // Reinit every segment so the reinit path runs concurrently too.
    flag = A.reinit(m, job->times[i], &yo);
    if (flag < 0) { job->status = flag; break; }
  }
  A.free(m);
}

}  // namespace

// [[Rcpp::export(".capi_test_concurrent")]]
List capi_test_concurrent(NumericVector times, double y0, double k1, double k2,
                          int solver = 0) {
  const capi_api& A = capi_select(solver);
  int n = times.size();
  capi_thread_job j1, j2;
  j1.api = &A; j2.api = &A;
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

// --- ARKODE method selection ------------------------------------------------
// Decay solved with an explicitly chosen ESDIRK: by order when order > 0,
// otherwise by table name. Returns the trajectory and the step count, so the
// tests can check both accuracy and that a higher order takes fewer steps.
// [[Rcpp::export(".capi_test_arkode_method")]]
List capi_test_arkode_method(NumericVector times, double y0, double k,
                             int order, std::string table) {
  int n = times.size();
  double kk = k;
  void* m = sundialr_arkode_create(1, &kk);
  if (!m) stop("sundialr_arkode_create returned NULL");

  sundialr_arkode_set_rhs(m, decay_rhs);
  sundialr_arkode_set_tol_scalar(m, 1e-8, 1e-10);
  sundialr_arkode_set_max_steps(m, 100000);
  int flag = order > 0 ? sundialr_arkode_set_order(m, order)
                       : sundialr_arkode_set_table_name(m, table.c_str());
  if (flag < 0) capi_test_fail(ARKODE_API, m, "method selection failed");

  NumericVector out(n);
  out[0] = y0;
  double ycur = y0;
  if (sundialr_arkode_reinit(m, times[0], &ycur) < 0) capi_test_fail(ARKODE_API, m, "reinit failed");
  for (int i = 1; i < n; i++) {
    double yo = 0.0, tr = 0.0;
    if (sundialr_arkode_solve(m, times[i], &yo, &tr) < 0) capi_test_fail(ARKODE_API, m, "solve failed");
    out[i] = yo;
  }
  long ns = sundialr_arkode_get_num_steps(m);
  sundialr_arkode_free(m);
  return List::create(_["y"] = out, _["steps"] = (double) ns);
}

// The method-selection error paths: an out-of-range order, an unknown table
// name, and a valid choice made after the integrator exists must each return
// SUNDIALR_ARK_ILL_INPUT with a message, and leave the handle usable.
// [[Rcpp::export(".capi_test_arkode_method_errors")]]
List capi_test_arkode_method_errors() {
  double kk = 1.0;
  void* m = sundialr_arkode_create(1, &kk);
  if (!m) stop("sundialr_arkode_create returned NULL");
  sundialr_arkode_set_rhs(m, decay_rhs);

  int bad_order = sundialr_arkode_set_order(m, 7);
  std::string bad_order_msg = sundialr_arkode_last_err(m) ? sundialr_arkode_last_err(m) : "";
  int bad_name = sundialr_arkode_set_table_name(m, "ARKODE_NO_SUCH_TABLE");
  std::string bad_name_msg = sundialr_arkode_last_err(m) ? sundialr_arkode_last_err(m) : "";

  double y = 1.0, yo = 0.0, tr = 0.0;
  int reinit = sundialr_arkode_reinit(m, 0.0, &y);
  int late = sundialr_arkode_set_order(m, 3);
  int solve = sundialr_arkode_solve(m, 1.0, &yo, &tr);
  sundialr_arkode_free(m);

  return List::create(_["bad_order"] = bad_order, _["bad_order_msg"] = bad_order_msg,
                      _["bad_name"] = bad_name,   _["bad_name_msg"] = bad_name_msg,
                      _["reinit"] = reinit, _["late"] = late,
                      _["solve"] = solve, _["y1"] = yo);
}

// --- ARKODE cold reinit ------------------------------------------------------
// One handle first solves an unrelated problem (a fast two-compartment system,
// p_a), then is pointed at a second one (p_b, from y_b at t0) by one of three
// routes: reinit_cold() on the used handle, warm reinit() on the used handle,
// and a NEW handle. reinit_cold() must reproduce the new handle bit for bit (no
// history crosses it); warm reinit() carries the first problem's step size and
// controller state into the second. Also returned: step counts (the cold handle
// keeps a running total) and whether 1000 further cold reinit+solve cycles kept
// the state buffer in place. tol_mode: 1 = scalar, 2 = per-equation vectors.
namespace {
struct cold_cfg { bool use_jac; int tol_mode; };

void* cold_handle(double* p, const cold_cfg& c) {
  void* m = sundialr_arkode_create(2, p);
  if (!m) stop("sundialr_arkode_create returned NULL");
  sundialr_arkode_set_rhs(m, twocmt_rhs);
  if (c.use_jac) sundialr_arkode_set_jac(m, twocmt_jac);
  if (c.tol_mode == 2) {
    double rt[2] = {1e-8, 1e-8}, at[2] = {1e-10, 1e-12};
    sundialr_arkode_set_tol_vector(m, rt, at);
  } else {
    sundialr_arkode_set_tol_scalar(m, 1e-8, 1e-10);
  }
  sundialr_arkode_set_max_steps(m, 100000);
  return m;
}

// Solve from the handle's current start to each of `times`; rows = times.
NumericMatrix cold_march(void* m, const NumericVector& times) {
  NumericMatrix out(times.size(), 2);
  for (int i = 0; i < times.size(); i++) {
    double yo[2], tr = 0.0;
    if (sundialr_arkode_solve(m, times[i], yo, &tr) < 0) capi_test_fail(ARKODE_API, m, "solve failed");
    out(i, 0) = yo[0]; out(i, 1) = yo[1];
  }
  return out;
}
}  // namespace

// [[Rcpp::export(".capi_test_arkode_cold")]]
List capi_test_arkode_cold(NumericVector times_a, NumericVector times_b,
                           NumericVector p_a, NumericVector p_b,
                           NumericVector y_a, NumericVector y_b,
                           bool use_jac, int tol_mode) {
  cold_cfg c{use_jac, tol_mode};
  double pa[2] = {p_a[0], p_a[1]}, pb[2] = {p_b[0], p_b[1]};
  double ya[2] = {y_a[0], y_a[1]}, yb[2] = {y_b[0], y_b[1]};
  const double t0b = times_b[0];
  NumericVector tb = times_b[Range(1, times_b.size() - 1)];

  // the problem-A history, then each route to problem B
  auto used = [&](void) {
    void* m = cold_handle(pa, c);
    if (sundialr_arkode_reinit(m, times_a[0], ya) < 0) capi_test_fail(ARKODE_API, m, "reinit failed");
    cold_march(m, times_a[Range(1, times_a.size() - 1)]);
    sundialr_arkode_set_udata(m, pb);
    return m;
  };

  void* mc = used();
  long steps_a = sundialr_arkode_get_num_steps(mc);
  if (sundialr_arkode_reinit_cold(mc, t0b, yb) < 0) capi_test_fail(ARKODE_API, mc, "reinit_cold failed");
  NumericMatrix cold = cold_march(mc, tb);
  long steps_cold_total = sundialr_arkode_get_num_steps(mc);
  // further cold cycles: the state buffer must stay where it is
  const void* p0 = sundialr_arkode_state_ptr_internal(mc);
  bool ptr_stable = true;
  for (int k = 0; k < 1000; k++) {
    double y[2] = {yb[0], yb[1]}, yo[2], tr = 0.0;
    if (sundialr_arkode_reinit_cold(mc, 0.0, y) < 0) capi_test_fail(ARKODE_API, mc, "reinit_cold failed");
    if (sundialr_arkode_solve(mc, 0.01, yo, &tr) < 0) capi_test_fail(ARKODE_API, mc, "solve failed");
    if (sundialr_arkode_state_ptr_internal(mc) != p0) ptr_stable = false;
  }
  sundialr_arkode_reset_stats(mc);
  long steps_after_reset = sundialr_arkode_get_num_steps(mc);
  sundialr_arkode_free(mc);

  void* mw = used();
  if (sundialr_arkode_reinit(mw, t0b, yb) < 0) capi_test_fail(ARKODE_API, mw, "reinit failed");
  NumericMatrix warm = cold_march(mw, tb);
  sundialr_arkode_free(mw);

  void* mf = cold_handle(pb, c);
  if (sundialr_arkode_reinit(mf, t0b, yb) < 0) capi_test_fail(ARKODE_API, mf, "reinit failed");
  NumericMatrix fresh = cold_march(mf, tb);
  long steps_fresh = sundialr_arkode_get_num_steps(mf);
  sundialr_arkode_free(mf);

  // reinit_cold before the integrator exists is an ordinary first reinit
  void* mn = cold_handle(pb, c);
  if (sundialr_arkode_reinit_cold(mn, t0b, yb) < 0) capi_test_fail(ARKODE_API, mn, "first reinit_cold failed");
  NumericMatrix first_cold = cold_march(mn, tb);
  sundialr_arkode_free(mn);

  return List::create(_["cold"] = cold, _["warm"] = warm, _["fresh"] = fresh,
                      _["first_cold"] = first_cold,
                      _["steps_a"] = (double) steps_a,
                      _["steps_cold_total"] = (double) steps_cold_total,
                      _["steps_fresh"] = (double) steps_fresh,
                      _["steps_after_reset"] = (double) steps_after_reset,
                      _["ptr_stable"] = ptr_stable);
}

// Robertson to each of `times`, from (1, 0, 0) at t = 0. Rows are output
// times, columns the three species.
// [[Rcpp::export(".capi_test_robertson")]]
NumericMatrix capi_test_robertson(NumericVector times, bool use_jac, int solver = 0) {
  const capi_api& A = capi_select(solver);
  int n = times.size();
  void* m = A.create(3, NULL);
  if (!m) stop("create returned NULL");

  A.set_rhs(m, robertson_rhs);
  if (use_jac) A.set_jac(m, robertson_jac);
  double rtol[3] = {1e-6, 1e-6, 1e-6};
  double atol[3] = {1e-10, 1e-14, 1e-8};
  A.set_tol_vector(m, rtol, atol);
  A.set_max_steps(m, 100000);

  double y[3] = {1.0, 0.0, 0.0};
  if (A.reinit(m, 0.0, y) < 0) capi_test_fail(A, m, "reinit failed");
  NumericMatrix out(n, 3);
  for (int i = 0; i < n; i++) {
    double yo[3], tr = 0.0;
    if (A.solve(m, times[i], yo, &tr) < 0) capi_test_fail(A, m, "solve failed");
    for (int j = 0; j < 3; j++) out(i, j) = yo[j];
  }
  A.free(m);
  return out;
}
