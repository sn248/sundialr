#ifndef SUNDIALS_SCOPE_GUARD_H
#define SUNDIALS_SCOPE_GUARD_H

// Runs a cleanup callable when the enclosing scope is left, including when it
// is left by a C++ exception.
//
// The solvers release their SUNDIALS objects (SUNContext, N_Vector, SUNMatrix,
// SUNLinearSolver, SUNNonlinearSolver, integrator memory) with explicit free
// calls. Rcpp's stop() throws, so any stop() between allocation and those free
// calls - input validation, a failed check_retval(), a failed solve - unwinds
// straight past them and leaks every object allocated so far. These are plain
// C allocations, so R's garbage collector never reclaims them. Wrapping the
// teardown in a guard makes it run on the normal return and on every throw.
//
// Usage: declare the resource variables first (initialised to NULL), then the
// guard, then allocate. The guard's body must therefore be NULL-tolerant, and
// should free in the same order the trailing free calls used to.
//
// Note: this covers C++ exceptions only. Rf_error() (used by
// sundials_r_err_handler) unwinds with longjmp, which does not run
// destructors, so SUNDIALS' own fatal-error path is not covered by this guard.

template <typename F>
class scope_guard {
public:
  explicit scope_guard(F f) : f_(f), active_(true) {}
  scope_guard(scope_guard&& other) : f_(other.f_), active_(other.active_) {
    other.active_ = false;
  }
  ~scope_guard() { if (active_) f_(); }

  scope_guard(const scope_guard&) = delete;
  scope_guard& operator=(const scope_guard&) = delete;

private:
  F f_;
  bool active_;
};

template <typename F>
scope_guard<F> make_scope_guard(F f) { return scope_guard<F>(f); }

#endif /* SUNDIALS_SCOPE_GUARD_H */
