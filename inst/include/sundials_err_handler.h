#ifndef SUNDIALS_ERR_HANDLER_H
#define SUNDIALS_ERR_HANDLER_H

// CRAN fix: CRAN policy requires that compiled code must not call abort() or
// write to stdout/stderr directly (see "Writing R Extensions" §1.6.4).
// SUNDIALS' default error handler calls abort() on fatal errors, which would
// terminate the R session. This header provides a replacement handler that
// routes fatal SUNDIALS errors through R's error mechanism instead.
//
// The handler does NOT raise the error itself. Rf_error() unwinds with
// longjmp, which does not run C++ destructors, so raising from inside the
// handler would skip the scope guards that release the SUNDIALS objects (see
// sundials_scope_guard.h) and leak them on every SUNDIALS error. Instead the
// handler records the error in a caller-supplied record and returns - SUNDIALS
// error handlers are expected to return, and the failing routine then reports
// the error through its return code as usual. The solver turns that into an
// R error with Rcpp's stop(), which throws a C++ exception and so unwinds
// through the guards.
//
// Usage: declare the record BEFORE the scope guard (so it outlives the
// SUNContext that refers to it), then after SUNContext_Create():
//   sundials_err_record sun_err;
//   ...
//   SUNContext_PushErrHandler(sunctx, sundials_r_err_handler, &sun_err);
// and raise failures with sundials_stop(sun_err, "fallback message") or
// sundials_check(sun_err) after calls whose return value is not inspected.

#include <string>
#include <Rcpp.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>

// Holds the first SUNDIALS error seen during one solver call. One record per
// call, passed through err_user_data, so no global mutable state is involved.
struct sundials_err_record {
  bool has_error;
  std::string message;

  sundials_err_record() : has_error(false) {}
};

static void sundials_r_err_handler(int line, const char* func, const char* file,
                                   const char* msg, SUNErrCode err_code,
                                   void* err_user_data, SUNContext sunctx) {
  sundials_err_record* rec = (sundials_err_record*) err_user_data;
  if (rec == NULL || rec->has_error) return;   // keep the first error

  rec->has_error = true;
  rec->message = std::string("SUNDIALS error in ") +
                 (func ? func : "unknown") + " (" +
                 (file ? file : "unknown") + ":" + std::to_string(line) + "): " +
                 (msg  ? msg  : "no message");
}

// Raise the error behind a failed SUNDIALS call. Throws, so the enclosing
// scope guards run. When the handler recorded a message it is used as is - it
// already names the function, file, line and cause. Otherwise (a routine that
// returned an error code without going through the handler) the call site's
// own description is raised, prefixed with the SUNDIALS function that failed,
// so the failing call is always identified.
static inline void sundials_stop(const sundials_err_record& rec,
                                 const char* funcname,
                                 const char* fallback) {
  if (rec.has_error) { Rcpp::stop(rec.message); }
  Rcpp::stop(std::string(funcname) + "() failed: " + fallback);
}

// Raise only if SUNDIALS recorded an error. Used after calls whose return
// value is not inspected, to preserve the previous fail-fast behaviour.
static inline void sundials_check(const sundials_err_record& rec) {
  if (rec.has_error) { Rcpp::stop(rec.message); }
}

#endif /* SUNDIALS_ERR_HANDLER_H */
