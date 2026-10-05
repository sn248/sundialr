//   Copyright (c) 2016-2026, Satyaprakash Nayak
//   Distributed under the BSD-3 licence; see the header of sundialr_capi.h.

// Implementation of the plain-C ARKODE (ESDIRK, dense) wrapper declared in
// sundialr_capi.h. It mirrors the CVODE half in sundialr_capi.cpp entry point
// for entry point, and the same overriding rule applies: NO entry point may
// throw, since each is called from a consumer's foreign C/C++ stack. Every body
// sits inside CAPI_GUARD, errors are recorded by the recording handler (never
// the raising sundials_stop()) and returned as codes, and the message is
// fetched with sundialr_arkode_last_err().
//
// The integrator is ARKStep run fully implicit (no explicit RHS), so each step
// is a diagonally implicit Runge-Kutta method solved with Newton iteration and
// a dense direct linear solver. The default table and every table selected by
// sundialr_arkode_set_order() is an ESDIRK; set_table_name() accepts any ARKODE
// DIRK table by name.
//
// Two differences from the CVODE half, both forced by ARKODE's API:
//   - ARKStepCreate needs t0 and y0, so arkode_mem is created by the first
//     reinit rather than by create(). Configuration made before that is stored
//     in the handle and applied then, as the CVODE half does for CVodeInit.
//   - Later reinits use ARKodeReset, which, unlike CVodeReInit, keeps the step
//     counter running. get_num_steps() is therefore the counter minus a base
//     that reset_stats() moves, with no banking at reinit.

#include <cstring>
#include <string>
#include <vector>

#include <arkode/arkode.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_butcher_dirk.h>
#include <arkode/arkode_butcher_erk.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>

#include <sundials_err_record.h>
#include <sundials_err_handler.h>

#include <sundialr_capi.h>

#include "capi_internal.h"

static_assert(sizeof(sunrealtype) == sizeof(double),
              "sundialr C API assumes sunrealtype is double");

// The header promises these numeric values without including arkode.h; keep
// the two in step at compile time.
static_assert(SUNDIALR_ARK_SUCCESS     == ARK_SUCCESS,     "ARK code mismatch");
static_assert(SUNDIALR_ARK_TOO_MUCH_WORK == ARK_TOO_MUCH_WORK, "ARK code mismatch");
static_assert(SUNDIALR_ARK_CONV_FAILURE == ARK_CONV_FAILURE, "ARK code mismatch");
static_assert(SUNDIALR_ARK_CONSTR_FAIL == ARK_CONSTR_FAIL, "ARK code mismatch");
static_assert(SUNDIALR_ARK_MEM_FAIL    == ARK_MEM_FAIL,    "ARK code mismatch");
static_assert(SUNDIALR_ARK_MEM_NULL    == ARK_MEM_NULL,    "ARK code mismatch");
static_assert(SUNDIALR_ARK_ILL_INPUT   == ARK_ILL_INPUT,   "ARK code mismatch");
static_assert(SUNDIALR_ARK_NO_MALLOC   == ARK_NO_MALLOC,   "ARK code mismatch");
static_assert(SUNDIALR_ARK_NLS_OP_ERR  == ARK_NLS_OP_ERR,  "ARK code mismatch");
static_assert(SUNDIALR_ARK_INVALID_TABLE == ARK_INVALID_TABLE, "ARK code mismatch");
static_assert(SUNDIALR_ARK_CONTROLLER_ERR == ARK_CONTROLLER_ERR, "ARK code mismatch");

// ---------------------------------------------------------------------------
// Handle
// ---------------------------------------------------------------------------
struct sundialr_arkode_handle {
  SUNContext      sunctx     = NULL;
  void*           arkode_mem = NULL;  // created by the first reinit
  N_Vector        y          = NULL;  // persistent state, length neq
  SUNMatrix       SM         = NULL;
  SUNLinearSolver LS         = NULL;

  int   neq   = 0;
  void* udata = NULL;

  sundialr_rhs rhs = NULL;
  sundialr_jac jac = NULL;

  bool initialized = false;           // has ARKStepCreate run yet?

  // Method. order > 0 selects ARKODE's default DIRK of that order (each an
  // ESDIRK for orders 2-5); otherwise table_id names the table. Both are fixed
  // once the integrator exists.
  int                order    = 0;
  ARKODE_DIRKTableID table_id = ARKODE_ESDIRK436L2SA_6_3_4;

  // ARKodeReset keeps the step counter running across reinits, so the count
  // reported is the counter minus this base, which reset_stats() moves.
  long nst_base = 0;

  // Tolerance configuration, as in the CVODE handle: 0 = none set (scalar
  // defaults), 1 = scalar, 2 = per-equation vectors via the error weights.
  int    tol_mode = 0;
  double rtol_s   = 1e-6;
  double atol_s   = 1e-8;
  std::vector<double> rtol_v;
  std::vector<double> atol_v;

  bool has_maxsteps = false; long   maxsteps = 0;
  bool has_hmax     = false; double hmax     = 0.0;
  bool has_hmin     = false; double hmin     = 0.0;

  sundials_err_record err;
  std::string last_err_msg;
};

// --- SUNDIALS callback thunks (C linkage, never throw) ---------------------

static int ark_rhs_thunk(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data) {
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) user_data;
  return h->rhs((double) t, N_VGetArrayPointer(y), N_VGetArrayPointer(ydot), h->udata);
}

static int ark_jac_thunk(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                         void* user_data, N_Vector t1, N_Vector t2, N_Vector t3) {
  (void) fy; (void) t1; (void) t2; (void) t3;
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) user_data;
  // Column-major dense storage, written in place (see sundialr_jac).
  return h->jac((double) t, N_VGetArrayPointer(y), SM_DATA_D(J), h->udata);
}

static int ark_ewt_thunk(N_Vector y, N_Vector w, void* user_data) {
  return capi_fill_ewt((sundialr_arkode_handle*) user_data, y, w);
}

static int ark_apply_tol(sundialr_arkode_handle* h) {
  if (h->tol_mode == 2) return ARKodeWFtolerances(h->arkode_mem, ark_ewt_thunk);
  return ARKodeSStolerances(h->arkode_mem, h->rtol_s, h->atol_s);
}

static int ark_apply_steps(sundialr_arkode_handle* h) {
  int flag;
  if (h->has_maxsteps) { flag = ARKodeSetMaxNumSteps(h->arkode_mem, h->maxsteps); if (flag < 0) return flag; }
  if (h->has_hmax)     { flag = ARKodeSetMaxStep(h->arkode_mem, h->hmax);         if (flag < 0) return flag; }
  if (h->has_hmin)     { flag = ARKodeSetMinStep(h->arkode_mem, h->hmin);         if (flag < 0) return flag; }
  return 0;
}

// The one-time ARKStep setup run by the first reinit. On failure the caller
// frees arkode_mem, so a later reinit starts again from a clean slate.
static int ark_first_init(sundialr_arkode_handle* h, double t0) {
  // Implicit only: fe = NULL, fi = the user's RHS.
  h->arkode_mem = ARKStepCreate(NULL, ark_rhs_thunk, t0, h->y, h->sunctx);
  if (!h->arkode_mem) return SUNDIALR_ARK_MEM_FAIL;

  int flag = ARKodeSetUserData(h->arkode_mem, h);
  if (flag < 0) return flag;
  if (h->order > 0) flag = ARKodeSetOrder(h->arkode_mem, h->order);
  else              flag = ARKStepSetTableNum(h->arkode_mem, h->table_id, ARKODE_ERK_NONE);
  if (flag < 0) return flag;
  flag = ARKodeSetLinearSolver(h->arkode_mem, h->LS, h->SM);
  if (flag < 0) return flag;
  if (h->jac) {
    flag = ARKodeSetJacFn(h->arkode_mem, ark_jac_thunk);   // must follow SetLinearSolver
    if (flag < 0) return flag;
  }
  flag = ark_apply_tol(h);
  if (flag < 0) return flag;
  return ark_apply_steps(h);
}

// Shared shape of the setters that may run before or after the integrator
// exists: store the value, and if it already exists apply it now.
#define ARK_SET_OR_APPLY(h, store, apply)                                      \
  CAPI_GUARD(h, SUNDIALR_ARK_MEM_FAIL, {                                       \
    store;                                                                     \
    if ((h)->initialized) {                                                    \
      int flag = (apply);                                                      \
      if (flag < 0) { capi_capture_err(h, flag, "ARKODE"); return flag; }     \
    }                                                                          \
    return SUNDIALR_ARK_SUCCESS;                                               \
  })

extern "C" {

// ---------------------------------------------------------------------------
// Lifecycle
// ---------------------------------------------------------------------------
void* sundialr_arkode_create(int neq, void* udata) {
  if (neq <= 0) return NULL;
  sundialr_arkode_handle* h = NULL;
  try {
    h = new sundialr_arkode_handle();
    h->neq   = neq;
    h->udata = udata;

    if (SUNContext_Create(SUN_COMM_NULL, &h->sunctx) < 0) { sundialr_arkode_free(h); return NULL; }
    SUNContext_PushErrHandler(h->sunctx, sundials_r_err_handler, &h->err);

    h->y = N_VNew_Serial(neq, h->sunctx);
    if (!h->y) { sundialr_arkode_free(h); return NULL; }

    h->SM = SUNDenseMatrix(neq, neq, h->sunctx);
    if (!h->SM) { sundialr_arkode_free(h); return NULL; }

    h->LS = SUNLinSol_Dense(h->y, h->SM, h->sunctx);
    if (!h->LS) { sundialr_arkode_free(h); return NULL; }
  }
  catch (...) {
    if (h) sundialr_arkode_free(h);
    return NULL;
  }
  return h;
}

void sundialr_arkode_free(void* m) {
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) m;
  if (!h) return;
  if (h->y)          N_VDestroy(h->y);
  if (h->LS)         SUNLinSolFree(h->LS);
  if (h->SM)         SUNMatDestroy(h->SM);
  if (h->arkode_mem) ARKodeFree(&h->arkode_mem);
  if (h->sunctx)     SUNContext_Free(&h->sunctx);
  delete h;
}

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------
int sundialr_arkode_set_rhs(void* m, sundialr_rhs f) {
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) m;
  if (!h) return SUNDIALR_ARK_MEM_NULL;
  h->rhs = f;
  return SUNDIALR_ARK_SUCCESS;
}

int sundialr_arkode_set_jac(void* m, sundialr_jac J) {
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) m;
  if (!h) return SUNDIALR_ARK_MEM_NULL;
  ARK_SET_OR_APPLY(h, h->jac = J,
                   ARKodeSetJacFn(h->arkode_mem, J ? ark_jac_thunk : NULL));
}

int sundialr_arkode_set_tol_scalar(void* m, double rtol, double atol) {
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) m;
  if (!h) return SUNDIALR_ARK_MEM_NULL;
  ARK_SET_OR_APPLY(h, (h->tol_mode = 1, h->rtol_s = rtol, h->atol_s = atol),
                   ARKodeSStolerances(h->arkode_mem, rtol, atol));
}

int sundialr_arkode_set_tol_vector(void* m, const double* rtol, const double* atol) {
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) m;
  if (!h) return SUNDIALR_ARK_MEM_NULL;
  if (!rtol || !atol) return SUNDIALR_ARK_ILL_INPUT;
  ARK_SET_OR_APPLY(h,
                   (h->tol_mode = 2,
                    h->rtol_v.assign(rtol, rtol + h->neq),
                    h->atol_v.assign(atol, atol + h->neq)),
                   ARKodeWFtolerances(h->arkode_mem, ark_ewt_thunk));
}

int sundialr_arkode_set_max_steps(void* m, long mxsteps) {
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) m;
  if (!h) return SUNDIALR_ARK_MEM_NULL;
  ARK_SET_OR_APPLY(h, (h->has_maxsteps = true, h->maxsteps = mxsteps),
                   ARKodeSetMaxNumSteps(h->arkode_mem, mxsteps));
}

int sundialr_arkode_set_max_step(void* m, double hmax) {
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) m;
  if (!h) return SUNDIALR_ARK_MEM_NULL;
  ARK_SET_OR_APPLY(h, (h->has_hmax = true, h->hmax = hmax),
                   ARKodeSetMaxStep(h->arkode_mem, hmax));
}

int sundialr_arkode_set_min_step(void* m, double hmin) {
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) m;
  if (!h) return SUNDIALR_ARK_MEM_NULL;
  ARK_SET_OR_APPLY(h, (h->has_hmin = true, h->hmin = hmin),
                   ARKodeSetMinStep(h->arkode_mem, hmin));
}

int sundialr_arkode_set_udata(void* m, void* udata) {
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) m;
  if (!h) return SUNDIALR_ARK_MEM_NULL;
  // The thunks read h->udata at call time (ARKodeSetUserData points at the
  // handle), so no ARKODE call is needed.
  h->udata = udata;
  return SUNDIALR_ARK_SUCCESS;
}

int sundialr_arkode_set_order(void* m, int order) {
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) m;
  if (!h) return SUNDIALR_ARK_MEM_NULL;
  CAPI_GUARD(h, SUNDIALR_ARK_MEM_FAIL, {
    if (h->initialized) {
      h->last_err_msg = "the method must be chosen before the first sundialr_arkode_reinit";
      return SUNDIALR_ARK_ILL_INPUT;
    }
    if (order < 2 || order > 5) {
      h->last_err_msg = "ESDIRK order must be 2, 3, 4 or 5, got " + std::to_string(order);
      return SUNDIALR_ARK_ILL_INPUT;
    }
    h->order = order;
    return SUNDIALR_ARK_SUCCESS;
  })
}

int sundialr_arkode_set_table_name(void* m, const char* name) {
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) m;
  if (!h) return SUNDIALR_ARK_MEM_NULL;
  CAPI_GUARD(h, SUNDIALR_ARK_MEM_FAIL, {
    if (h->initialized) {
      h->last_err_msg = "the method must be chosen before the first sundialr_arkode_reinit";
      return SUNDIALR_ARK_ILL_INPUT;
    }
    if (!name) return SUNDIALR_ARK_ILL_INPUT;
    // Match against the valid IDs here rather than calling
    // ARKStepSetTableName: ARKODE reports an unknown name through its error
    // processor with no memory block, which bypasses the handle's recording
    // handler. The enum is contiguous between its MIN and MAX markers.
    for (int id = ARKODE_MIN_DIRK_NUM; id <= ARKODE_MAX_DIRK_NUM; id++) {
      const char* nm = ARKodeButcherTable_DIRKIDToName((ARKODE_DIRKTableID) id);
      if (nm && std::strcmp(nm, name) == 0) {
        h->table_id = (ARKODE_DIRKTableID) id;
        h->order = 0;
        return SUNDIALR_ARK_SUCCESS;
      }
    }
    h->last_err_msg = std::string("unknown ARKODE DIRK table '") + name + "'";
    return SUNDIALR_ARK_ILL_INPUT;
  })
}

// ---------------------------------------------------------------------------
// Integration
// ---------------------------------------------------------------------------
int sundialr_arkode_reinit(void* m, double t0, const double* y0) {
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) m;
  if (!h) return SUNDIALR_ARK_MEM_NULL;
  CAPI_GUARD(h, SUNDIALR_ARK_MEM_FAIL, {
    capi_clear_err(h);
    double* yp = N_VGetArrayPointer(h->y);
    for (int i = 0; i < h->neq; i++) yp[i] = y0[i];

    if (!h->initialized) {
      if (!h->rhs) {
        h->last_err_msg = "RHS function not set (call sundialr_arkode_set_rhs before reinit)";
        return SUNDIALR_ARK_ILL_INPUT;
      }
      int flag = ark_first_init(h, t0);
      if (flag < 0) {
        capi_capture_err(h, flag, "ARKODE");
        if (h->arkode_mem) ARKodeFree(&h->arkode_mem);
        return flag;
      }
      h->initialized = true;
    } else {
      int flag = ARKodeReset(h->arkode_mem, t0, h->y);
      if (flag < 0) { capi_capture_err(h, flag, "ARKODE"); return flag; }
    }
    return SUNDIALR_ARK_SUCCESS;
  })
}

int sundialr_arkode_solve(void* m, double tout, double* y, double* treached) {
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) m;
  if (!h) return SUNDIALR_ARK_MEM_NULL;
  CAPI_GUARD(h, SUNDIALR_ARK_MEM_FAIL, {
    if (!h->initialized) {
      h->last_err_msg = "solve called before reinit";
      return SUNDIALR_ARK_NO_MALLOC;
    }
    capi_clear_err(h);
    sunrealtype treal = 0.0;
    int flag = ARKodeEvolve(h->arkode_mem, tout, h->y, &treal, ARK_NORMAL);

    double* yp = N_VGetArrayPointer(h->y);
    for (int i = 0; i < h->neq; i++) y[i] = yp[i];
    if (treached) *treached = (double) treal;

    if (flag < 0) capi_capture_err(h, flag, "ARKODE");
    return flag;
  })
}

// ---------------------------------------------------------------------------
// Introspection
// ---------------------------------------------------------------------------
long sundialr_arkode_get_num_steps(void* m) {
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) m;
  if (!h || !h->initialized) return -1;
  long ns = 0;
  if (ARKodeGetNumSteps(h->arkode_mem, &ns) < 0) return -1;
  return ns - h->nst_base;
}

int sundialr_arkode_reset_stats(void* m) {
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) m;
  if (!h) return SUNDIALR_ARK_MEM_NULL;
  if (!h->initialized) { h->nst_base = 0; return SUNDIALR_ARK_SUCCESS; }
  long ns = 0;
  if (ARKodeGetNumSteps(h->arkode_mem, &ns) < 0) return SUNDIALR_ARK_MEM_FAIL;
  h->nst_base = ns;
  return SUNDIALR_ARK_SUCCESS;
}

// White-box hook for the tight-loop test, as sundialr_cvode_state_ptr_internal:
// not in the public header and not registered.
const void* sundialr_arkode_state_ptr_internal(void* m) {
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) m;
  if (!h || !h->y) return NULL;
  return (const void*) N_VGetArrayPointer(h->y);
}

const char* sundialr_arkode_last_err(void* m) {
  sundialr_arkode_handle* h = (sundialr_arkode_handle*) m;
  if (!h) return NULL;
  if (h->last_err_msg.empty()) return NULL;
  return h->last_err_msg.c_str();
}

}  // extern "C"
