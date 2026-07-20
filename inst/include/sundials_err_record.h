#ifndef SUNDIALS_ERR_RECORD_H
#define SUNDIALS_ERR_RECORD_H

// Where a solver call collects the first error it meets, whether that error
// came from SUNDIALS itself (see sundials_err_handler.h) or from a callback
// into R (see sundials_callback_guard below), and the helpers that turn it
// into an R error.
//
// This is kept separate from sundials_err_handler.h so that translation units
// needing only the record - rhs_func.cpp, for instance - do not also pull in
// the handler function and warn about it being unused.

#include <string>
#include <Rcpp.h>

// One record per solver call, reached through err_user_data or through the
// user-data struct, so no global mutable state is involved and concurrent
// solves stay independent.
struct sundials_err_record {
  bool has_error;
  std::string message;

  sundials_err_record() : has_error(false) {}

  // Keeps the first error: it is the cause, later ones are consequences.
  void record(const std::string& msg) {
    if (!has_error) { has_error = true; message = msg; }
  }
};

// Raise the recorded error if there is one, otherwise raise fallback prefixed
// with the SUNDIALS routine that failed. Throws, so enclosing scope guards run.
static inline void sundials_stop(const sundials_err_record& rec,
                                 const char* funcname,
                                 const char* fallback) {
  if (rec.has_error) { Rcpp::stop(rec.message); }
  Rcpp::stop(std::string(funcname) + "() failed: " + fallback);
}

// Raise only if an error was recorded. Used after calls whose return value is
// not inspected, to preserve the previous fail-fast behaviour.
static inline void sundials_check(const sundials_err_record& rec) {
  if (rec.has_error) { Rcpp::stop(rec.message); }
}

// Runs the body of a SUNDIALS callback and converts any C++ exception it
// raises into a recorded message plus SUNDIALS' unrecoverable-error return
// code, so the exception never propagates out through SUNDIALS' C frames.
// The failing routine then returns an error code as usual and the solver
// reports the recorded message through sundials_stop(). This also lets
// SUNDIALS unwind its own internal state instead of being abandoned mid-step.
//
// R-level errors and user interrupts are deliberately re-thrown. Rcpp raises
// those as LongjumpException / InterruptedException carrying an R
// continuation, and only Rcpp's handler at the .Call boundary can resume one.
// Swallowing it would discard the R condition and break tryCatch() on the R
// side, so those keep using Rcpp's own mechanism. They are C++ exceptions, so
// the solver's scope guards still run and nothing leaks.
template <typename F>
static inline int sundials_callback_guard(sundials_err_record* rec, F body) {
  try {
    return body();
  }
  catch (Rcpp::LongjumpException&) { throw; }
  catch (Rcpp::internal::InterruptedException&) { throw; }
  catch (std::exception& e) {
    if (rec) { rec->record(e.what()); }
    return -1;
  }
  catch (...) {
    if (rec) { rec->record("unidentified C++ exception raised in a SUNDIALS callback"); }
    return -1;
  }
}

#endif /* SUNDIALS_ERR_RECORD_H */
