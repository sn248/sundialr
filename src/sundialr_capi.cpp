//   Copyright (c) 2016-2026, Satyaprakash Nayak
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions are
//   met:
//
//   Redistributions of source code must retain the above copyright
//   notice, this list of conditions and the following disclaimer.
//
//   Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in
//   the documentation and/or other materials provided with the
//   distribution.
//
//   Neither sundialr nor the names of its
//   contributors may be used to endorse or promote products derived
//   from this software without specific prior written permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// Implementation of the plain-C CVODE wrapper declared in sundialr_capi.h. See
// that header for the contract. The overriding rule here is that NO entry point
// may throw: each one is called from a consumer's foreign C/C++ stack, so a C++
// exception crossing the boundary would unwind through frames that cannot
// handle it. Errors are recorded (reusing the package's recording error
// handler, never the raising sundials_stop()) and returned as codes; the
// message is fetched with sundialr_cvode_last_err(). This is the same class of
// defect as the Rf_error-longjmp problem fixed in 0.1.8/0.2.0 - see the note in
// sundials_err_handler.h - and is just as easy to reintroduce by copying a
// solver body that uses stop(), so keep every body inside CAPI_GUARD.

#include <cmath>
#include <string>
#include <vector>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>

// Reused, unchanged, from the R-facing solvers: the recording error handler and
// its record. We use ONLY the recording handler here, never sundials_stop().
#include <sundials_err_record.h>
#include <sundials_err_handler.h>

#include <sundialr_capi.h>

// The C API speaks `double`; N_VGetArrayPointer hands back sunrealtype*. The
// package builds SUNDIALS at its default double precision, so the two coincide
// and the state buffers can be shared without conversion. Fail loudly at compile
// time if a future build ever changes that, rather than silently reinterpreting.
static_assert(sizeof(sunrealtype) == sizeof(double),
              "sundialr C API assumes sunrealtype is double");

// ---------------------------------------------------------------------------
// Handle
// ---------------------------------------------------------------------------
struct sundialr_cvode_handle {
  SUNContext      sunctx    = NULL;
  void*           cvode_mem = NULL;
  N_Vector        y         = NULL;   // persistent state, length neq
  SUNMatrix       SM        = NULL;
  SUNLinearSolver LS        = NULL;

  int   neq   = 0;
  void* udata = NULL;

  sundialr_rhs rhs = NULL;
  sundialr_jac jac = NULL;

  bool initialized = false;           // has CVodeInit run yet?

  // Tolerance configuration, applied on the first reinit and immediately on any
  // later change. 0 = none set yet (a scalar default is used), 1 = scalar,
  // 2 = per-equation vectors (via the error-weight function).
  int    tol_mode = 0;
  double rtol_s   = 1e-6;
  double atol_s   = 1e-8;
  std::vector<double> rtol_v;
  std::vector<double> atol_v;

  // Step-size configuration; applied when set (if already initialized) and on
  // the first reinit.
  bool has_maxsteps = false; long   maxsteps = 0;
  bool has_hmax     = false; double hmax     = 0.0;
  bool has_hmin     = false; double hmin     = 0.0;

  sundials_err_record err;            // filled by the recording handler
  std::string last_err_msg;           // returned by last_err()
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Copy the current failure into last_err_msg: the SUNDIALS-recorded message if
// there is one, otherwise a description naming the numeric code.
static void capi_capture_err(sundialr_cvode_handle* h, int flag) {
  if (h->err.has_error) h->last_err_msg = h->err.message;
  else h->last_err_msg = std::string("CVODE returned error code ") + std::to_string(flag);
}

// Clear the per-call error state before an operation that can record one.
static void capi_clear_err(sundialr_cvode_handle* h) {
  h->err.has_error = false;
  h->err.message.clear();
}

// --- SUNDIALS callback thunks (C linkage, never throw) ---------------------

static int capi_rhs_thunk(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data) {
  sundialr_cvode_handle* h = (sundialr_cvode_handle*) user_data;
  double* yp   = N_VGetArrayPointer(y);
  double* ydp  = N_VGetArrayPointer(ydot);
  return h->rhs((double) t, yp, ydp, h->udata);
}

static int capi_jac_thunk(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                          void* user_data, N_Vector t1, N_Vector t2, N_Vector t3) {
  (void) fy; (void) t1; (void) t2; (void) t3;
  sundialr_cvode_handle* h = (sundialr_cvode_handle*) user_data;
  double* yp = N_VGetArrayPointer(y);
  // SM_DATA_D is the dense matrix's own column-major buffer; the consumer writes
  // J[i + j*neq] straight into it (see the sundialr_jac contract in the header).
  double* Jdata = SM_DATA_D(J);
  return h->jac((double) t, yp, Jdata, h->udata);
}

// Error-weight function realising per-equation rtol/atol (CVODE itol == 4).
static int capi_ewt_thunk(N_Vector y, N_Vector w, void* user_data) {
  sundialr_cvode_handle* h = (sundialr_cvode_handle*) user_data;
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

// Apply the stored tolerance configuration to an initialized handle.
static int capi_apply_tol(sundialr_cvode_handle* h) {
  if (h->tol_mode == 2) return CVodeWFtolerances(h->cvode_mem, capi_ewt_thunk);
  if (h->tol_mode == 1) return CVodeSStolerances(h->cvode_mem, h->rtol_s, h->atol_s);
  return CVodeSStolerances(h->cvode_mem, h->rtol_s, h->atol_s);  // defaults
}

// Apply the stored step-size configuration to an initialized handle.
static int capi_apply_steps(sundialr_cvode_handle* h) {
  int flag;
  if (h->has_maxsteps) { flag = CVodeSetMaxNumSteps(h->cvode_mem, h->maxsteps); if (flag < 0) return flag; }
  if (h->has_hmax)     { flag = CVodeSetMaxStep(h->cvode_mem, h->hmax);         if (flag < 0) return flag; }
  if (h->has_hmin)     { flag = CVodeSetMinStep(h->cvode_mem, h->hmin);         if (flag < 0) return flag; }
  return 0;
}

// Wrap every entry-point body so no C++ exception can escape to the caller (the
// whole reason this file exists). A thrown std::exception is recorded and turned
// into an error code; anything else becomes a generic memory error.
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

extern "C" {

// ---------------------------------------------------------------------------
// Lifecycle
// ---------------------------------------------------------------------------
void* sundialr_cvode_create(int neq, void* udata) {
  if (neq <= 0) return NULL;
  sundialr_cvode_handle* h = NULL;
  try {
    h = new sundialr_cvode_handle();
    h->neq   = neq;
    h->udata = udata;

    if (SUNContext_Create(SUN_COMM_NULL, &h->sunctx) < 0) { sundialr_cvode_free(h); return NULL; }
    // Record SUNDIALS errors into the handle instead of aborting (CRAN rule);
    // crucially the handler only records, so nothing longjmps or throws here.
    SUNContext_PushErrHandler(h->sunctx, sundials_r_err_handler, &h->err);

    h->cvode_mem = CVodeCreate(CV_BDF, h->sunctx);
    if (!h->cvode_mem) { sundialr_cvode_free(h); return NULL; }

    h->y = N_VNew_Serial(neq, h->sunctx);
    if (!h->y) { sundialr_cvode_free(h); return NULL; }

    h->SM = SUNDenseMatrix(neq, neq, h->sunctx);
    if (!h->SM) { sundialr_cvode_free(h); return NULL; }

    h->LS = SUNLinSol_Dense(h->y, h->SM, h->sunctx);
    if (!h->LS) { sundialr_cvode_free(h); return NULL; }
  }
  catch (...) {
    if (h) sundialr_cvode_free(h);
    return NULL;
  }
  return h;
}

void sundialr_cvode_free(void* m) {
  sundialr_cvode_handle* h = (sundialr_cvode_handle*) m;
  if (!h) return;
  // Same order as the R solvers' cleanup: dependents first, context last.
  if (h->y)         N_VDestroy(h->y);
  if (h->LS)        SUNLinSolFree(h->LS);
  if (h->SM)        SUNMatDestroy(h->SM);
  if (h->cvode_mem) CVodeFree(&h->cvode_mem);
  if (h->sunctx)    SUNContext_Free(&h->sunctx);
  delete h;
}

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------
int sundialr_cvode_set_rhs(void* m, sundialr_rhs f) {
  sundialr_cvode_handle* h = (sundialr_cvode_handle*) m;
  if (!h) return SUNDIALR_CV_MEM_NULL;
  h->rhs = f;
  return SUNDIALR_CV_SUCCESS;
}

int sundialr_cvode_set_jac(void* m, sundialr_jac J) {
  sundialr_cvode_handle* h = (sundialr_cvode_handle*) m;
  if (!h) return SUNDIALR_CV_MEM_NULL;
  CAPI_GUARD(h, SUNDIALR_CV_MEM_FAIL, {
    h->jac = J;
    if (h->initialized) {
      int flag = J ? CVodeSetJacFn(h->cvode_mem, capi_jac_thunk)
                   : CVodeSetJacFn(h->cvode_mem, NULL);
      if (flag < 0) { capi_capture_err(h, flag); return flag; }
    }
    return SUNDIALR_CV_SUCCESS;
  })
}

int sundialr_cvode_set_tol_scalar(void* m, double rtol, double atol) {
  sundialr_cvode_handle* h = (sundialr_cvode_handle*) m;
  if (!h) return SUNDIALR_CV_MEM_NULL;
  CAPI_GUARD(h, SUNDIALR_CV_MEM_FAIL, {
    h->tol_mode = 1;
    h->rtol_s = rtol;
    h->atol_s = atol;
    if (h->initialized) {
      int flag = CVodeSStolerances(h->cvode_mem, rtol, atol);
      if (flag < 0) { capi_capture_err(h, flag); return flag; }
    }
    return SUNDIALR_CV_SUCCESS;
  })
}

int sundialr_cvode_set_tol_vector(void* m, const double* rtol, const double* atol) {
  sundialr_cvode_handle* h = (sundialr_cvode_handle*) m;
  if (!h) return SUNDIALR_CV_MEM_NULL;
  if (!rtol || !atol) return SUNDIALR_CV_ILL_INPUT;
  CAPI_GUARD(h, SUNDIALR_CV_MEM_FAIL, {
    h->tol_mode = 2;
    h->rtol_v.assign(rtol, rtol + h->neq);
    h->atol_v.assign(atol, atol + h->neq);
    if (h->initialized) {
      int flag = CVodeWFtolerances(h->cvode_mem, capi_ewt_thunk);
      if (flag < 0) { capi_capture_err(h, flag); return flag; }
    }
    return SUNDIALR_CV_SUCCESS;
  })
}

int sundialr_cvode_set_max_steps(void* m, long mxsteps) {
  sundialr_cvode_handle* h = (sundialr_cvode_handle*) m;
  if (!h) return SUNDIALR_CV_MEM_NULL;
  CAPI_GUARD(h, SUNDIALR_CV_MEM_FAIL, {
    h->has_maxsteps = true;
    h->maxsteps = mxsteps;
    if (h->initialized) {
      int flag = CVodeSetMaxNumSteps(h->cvode_mem, mxsteps);
      if (flag < 0) { capi_capture_err(h, flag); return flag; }
    }
    return SUNDIALR_CV_SUCCESS;
  })
}

int sundialr_cvode_set_max_step(void* m, double hmax) {
  sundialr_cvode_handle* h = (sundialr_cvode_handle*) m;
  if (!h) return SUNDIALR_CV_MEM_NULL;
  CAPI_GUARD(h, SUNDIALR_CV_MEM_FAIL, {
    h->has_hmax = true;
    h->hmax = hmax;
    if (h->initialized) {
      int flag = CVodeSetMaxStep(h->cvode_mem, hmax);
      if (flag < 0) { capi_capture_err(h, flag); return flag; }
    }
    return SUNDIALR_CV_SUCCESS;
  })
}

int sundialr_cvode_set_min_step(void* m, double hmin) {
  sundialr_cvode_handle* h = (sundialr_cvode_handle*) m;
  if (!h) return SUNDIALR_CV_MEM_NULL;
  CAPI_GUARD(h, SUNDIALR_CV_MEM_FAIL, {
    h->has_hmin = true;
    h->hmin = hmin;
    if (h->initialized) {
      int flag = CVodeSetMinStep(h->cvode_mem, hmin);
      if (flag < 0) { capi_capture_err(h, flag); return flag; }
    }
    return SUNDIALR_CV_SUCCESS;
  })
}

// ---------------------------------------------------------------------------
// Integration
// ---------------------------------------------------------------------------
int sundialr_cvode_reinit(void* m, double t0, const double* y0) {
  sundialr_cvode_handle* h = (sundialr_cvode_handle*) m;
  if (!h) return SUNDIALR_CV_MEM_NULL;
  CAPI_GUARD(h, SUNDIALR_CV_MEM_FAIL, {
    capi_clear_err(h);
    double* yp = N_VGetArrayPointer(h->y);
    for (int i = 0; i < h->neq; i++) yp[i] = y0[i];

    int flag;
    if (!h->initialized) {
      if (!h->rhs) {
        h->last_err_msg = "RHS function not set (call sundialr_cvode_set_rhs before reinit)";
        return SUNDIALR_CV_ILL_INPUT;
      }
      // One-time setup. CVodeReInit does not reset any of these, so later
      // reinits keep the linear solver, tolerances, Jacobian and step limits.
      flag = CVodeInit(h->cvode_mem, capi_rhs_thunk, t0, h->y);
      if (flag < 0) { capi_capture_err(h, flag); return flag; }
      flag = CVodeSetUserData(h->cvode_mem, h);
      if (flag < 0) { capi_capture_err(h, flag); return flag; }
      flag = CVodeSetLinearSolver(h->cvode_mem, h->LS, h->SM);
      if (flag < 0) { capi_capture_err(h, flag); return flag; }
      if (h->jac) {
        flag = CVodeSetJacFn(h->cvode_mem, capi_jac_thunk);   // must follow SetLinearSolver
        if (flag < 0) { capi_capture_err(h, flag); return flag; }
      }
      flag = capi_apply_tol(h);
      if (flag < 0) { capi_capture_err(h, flag); return flag; }
      flag = capi_apply_steps(h);
      if (flag < 0) { capi_capture_err(h, flag); return flag; }
      h->initialized = true;
    } else {
      flag = CVodeReInit(h->cvode_mem, t0, h->y);
      if (flag < 0) { capi_capture_err(h, flag); return flag; }
    }
    return SUNDIALR_CV_SUCCESS;
  })
}

int sundialr_cvode_solve(void* m, double tout, double* y, double* treached) {
  sundialr_cvode_handle* h = (sundialr_cvode_handle*) m;
  if (!h) return SUNDIALR_CV_MEM_NULL;
  CAPI_GUARD(h, SUNDIALR_CV_MEM_FAIL, {
    if (!h->initialized) {
      h->last_err_msg = "solve called before reinit";
      return SUNDIALR_CV_NO_MALLOC;
    }
    capi_clear_err(h);
    sunrealtype treal = 0.0;
    int flag = CVode(h->cvode_mem, tout, h->y, &treal, CV_NORMAL);

    double* yp = N_VGetArrayPointer(h->y);
    for (int i = 0; i < h->neq; i++) y[i] = yp[i];
    if (treached) *treached = (double) treal;

    if (flag < 0) capi_capture_err(h, flag);
    return flag;
  })
}

// ---------------------------------------------------------------------------
// Introspection
// ---------------------------------------------------------------------------
long sundialr_cvode_get_num_steps(void* m) {
  sundialr_cvode_handle* h = (sundialr_cvode_handle*) m;
  if (!h || !h->initialized) return -1;
  long ns = 0;
  if (CVodeGetNumSteps(h->cvode_mem, &ns) < 0) return -1;
  return ns;
}

const char* sundialr_cvode_last_err(void* m) {
  sundialr_cvode_handle* h = (sundialr_cvode_handle*) m;
  if (!h) return NULL;
  if (h->last_err_msg.empty()) return NULL;
  return h->last_err_msg.c_str();
}

int sundialr_abi_version(void) {
  return SUNDIALR_ABI_VERSION;
}

}  // extern "C"
