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
// and raise failures with sundials_stop(sun_err, "SUNDIALSRoutine", "fallback")
// or sundials_check(sun_err) after calls whose return value is not inspected.
//
// The record itself, and the helpers that raise from it, live in
// sundials_err_record.h so that files needing only those do not pull in this
// handler as well.

#include <sundials_err_record.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>

static void sundials_r_err_handler(int line, const char* func, const char* file,
                                   const char* msg, SUNErrCode err_code,
                                   void* err_user_data, SUNContext sunctx) {
  sundials_err_record* rec = (sundials_err_record*) err_user_data;
  if (rec == NULL) return;

  rec->record(std::string("SUNDIALS error in ") +
              (func ? func : "unknown") + " (" +
              (file ? file : "unknown") + ":" + std::to_string(line) + "): " +
              (msg  ? msg  : "no message"));
}

#endif /* SUNDIALS_ERR_HANDLER_H */
