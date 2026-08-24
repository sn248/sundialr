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

#ifndef SUNDIALR_CAPI_H
#define SUNDIALR_CAPI_H

/*
 * sundialr C API - a plain-C wrapper around CVODE (BDF, dense) intended to be
 * called from another package's compiled code via
 * R_GetCCallable("sundialr", ...). This is the ONLY header a consumer includes:
 * it deliberately exposes no SUNDIALS types and no Rcpp, only plain double and
 * int arguments and an opaque handle, so the consumer needs neither SUNDIALS
 * headers nor a LinkingTo.
 *
 * Contract (see the notes on each group below):
 *   - Every entry point returns a status code; NONE ever throws or calls
 *     Rf_error(). A C++ exception must not cross this boundary, since it would
 *     unwind into the caller's foreign stack. On failure the code is returned
 *     and a human-readable message is retrievable with sundialr_cvode_last_err().
 *   - Status codes are CVODE's own return codes, mirrored below with a
 *     SUNDIALR_ prefix so a consumer can switch on them without including
 *     cvode.h. The numeric values match CVODE's CV_* macros exactly.
 *   - The handle owns a persistent state N_Vector of length neq. reinit() sets
 *     the state and (re)starts the integrator at t0; solve() advances to tout
 *     and copies the state into the caller's buffer. One handle is meant to be
 *     reused across many segments via reinit(), the supported CVODE pattern.
 *
 * Reuse and allocation guarantee:
 *   After the first reinit()+solve() cycle on a handle, further reinit() and
 *   solve() calls perform no heap allocation. The wrapper's reinit copies y0
 *   into the existing state vector and calls CVodeReInit, which SUNDIALS
 *   documents (and its source confirms) as resetting counters and reusing the
 *   memory allocated by CVodeInit; the linear solver's saved-Jacobian matrix is
 *   cloned once, on the first solve, and reused thereafter. So a host loop
 *   issuing 1e4-1e5 small reinit+solve segments per fit pays the allocation
 *   cost once at create()/first solve, never per segment. The one exception is
 *   a FAILING call, which may allocate to store the message that
 *   sundialr_cvode_last_err() returns. This guarantee is pinned by the tight-
 *   loop test in tests/testthat/test-capi.r.
 *
 * Thread-safety:
 *   One handle per concurrent solve. Every piece of mutable state - the
 *   SUNContext (SUNDIALS 6+ is designed around per-thread contexts), cvode_mem,
 *   state vector, matrix/linear solver, tolerance/step configuration, udata
 *   pointer, error record and message - lives inside the handle; the file has
 *   no global mutable state, the callback thunks are stateless, and no entry
 *   point calls the R API. Distinct handles may therefore be driven from
 *   different threads simultaneously (each thread its own handle), as parallel
 *   host workers do. A single handle must NOT be shared between threads without
 *   external locking. Callbacks run on the thread that called solve(); a
 *   callback that touches R must only ever run on the main R thread.
 *
 * Typical use:
 *   void* m = sundialr_cvode_create(neq, udata);
 *   sundialr_cvode_set_rhs(m, my_rhs);
 *   sundialr_cvode_set_tol_scalar(m, 1e-6, 1e-8);
 *   sundialr_cvode_set_max_steps(m, 5000);
 *   for each subject:                          // one handle serves them all
 *     sundialr_cvode_set_udata(m, subj);       // repoint the callback data
 *     sundialr_cvode_reset_stats(m);           // per-subject step counts
 *     for each segment:
 *       sundialr_cvode_reinit(m, t0, y0);      // only when the state jumps
 *       sundialr_cvode_solve (m, tout, y, &tr); // check return < 0
 *   sundialr_cvode_free(m);
 *
 * History: sundialr_cvode_set_udata and sundialr_cvode_reset_stats were added
 * in sundialr 0.2.1. Both are additive and binary-compatible, so
 * SUNDIALR_ABI_VERSION stays 1 (it moves only on a breaking change); a consumer
 * that needs them should require sundialr >= 0.2.1.
 */

#ifdef __cplusplus
extern "C" {
#endif

/* --- Status codes: numeric values identical to CVODE's CV_* macros --------- */
#define SUNDIALR_CV_SUCCESS            0
#define SUNDIALR_CV_TSTOP_RETURN       1
#define SUNDIALR_CV_ROOT_RETURN        2
#define SUNDIALR_CV_TOO_MUCH_WORK     -1
#define SUNDIALR_CV_TOO_MUCH_ACC      -2
#define SUNDIALR_CV_ERR_FAILURE       -3
#define SUNDIALR_CV_CONV_FAILURE      -4
#define SUNDIALR_CV_LINIT_FAIL        -5
#define SUNDIALR_CV_LSETUP_FAIL       -6
#define SUNDIALR_CV_LSOLVE_FAIL       -7
#define SUNDIALR_CV_RHSFUNC_FAIL      -8
#define SUNDIALR_CV_FIRST_RHSFUNC_ERR -9
#define SUNDIALR_CV_REPTD_RHSFUNC_ERR -10
#define SUNDIALR_CV_UNREC_RHSFUNC_ERR -11
#define SUNDIALR_CV_CONSTR_FAIL       -15
#define SUNDIALR_CV_MEM_FAIL          -20
#define SUNDIALR_CV_MEM_NULL          -21
#define SUNDIALR_CV_ILL_INPUT         -22
#define SUNDIALR_CV_NO_MALLOC         -23
#define SUNDIALR_CV_TOO_CLOSE         -27

/* Incremented on any binary-incompatible change to this header. A consumer
 * should compare its compiled-against value with sundialr_abi_version(). */
#define SUNDIALR_ABI_VERSION 1

/* --- Callback types -------------------------------------------------------- *
 * Both are called from CVODE's integration loop. They must return 0 on success
 * and non-zero to signal a recoverable (>0) or unrecoverable (<0) failure, and
 * must not throw. `udata` is the pointer passed to sundialr_cvode_create().    */

/* Right-hand side: fill ydot[i] = d(y[i])/dt, i = 0..neq-1. */
typedef int (*sundialr_rhs)(double t, const double* y, double* ydot, void* udata);

/* Jacobian: fill the dense neq-by-neq matrix J in COLUMN-MAJOR order, i.e.
 * J[i + j*neq] = d(ydot[i])/d(y[j]). Column-major matches CVODE's dense matrix
 * storage, so the buffer handed in is the matrix's own data - write it in place. */
typedef int (*sundialr_jac)(double t, const double* y, double* J, void* udata);

/* --- Lifecycle ------------------------------------------------------------- */

/* Allocate a handle for a system of neq equations. udata is passed unchanged to
 * every callback. Returns NULL on allocation failure. */
void* sundialr_cvode_create(int neq, void* udata);

/* Release the handle and every SUNDIALS object it owns. NULL-tolerant. */
void  sundialr_cvode_free(void* m);

/* --- Configuration (call after create, before or between solves) ----------- */

/* Required before the first reinit. */
int sundialr_cvode_set_rhs(void* m, sundialr_rhs f);

/* Optional analytic Jacobian; pass NULL to fall back to CVODE's dense
 * finite-difference approximation (the default). */
int sundialr_cvode_set_jac(void* m, sundialr_jac J);

/* Scalar relative and absolute tolerance (CVODE itol == 1). */
int sundialr_cvode_set_tol_scalar(void* m, double rtol, double atol);

/* Per-equation relative AND absolute tolerance (CVODE itol == 4). CVODE has no
 * vector-rtol tolerance setter, so this is realised with a custom error-weight
 * function w[i] = 1/(rtol[i]*|y[i]| + atol[i]); rtol and atol must each point to
 * neq values, which are copied into the handle. */
int sundialr_cvode_set_tol_vector(void* m, const double* rtol, const double* atol);

/* Maximum internal steps between two solve() output points (CVODE default 500). */
int sundialr_cvode_set_max_steps(void* m, long mxsteps);

/* Upper/lower bound on the internal step size. */
int sundialr_cvode_set_max_step(void* m, double hmax);
int sundialr_cvode_set_min_step(void* m, double hmin);

/* Replace the udata pointer handed to every callback. Takes effect from the
 * next callback invocation; no reinit is required (though a host loop will
 * typically follow it with one, since new callback data usually means a new
 * trajectory). Lets one handle serve many parameter sets - e.g. a population
 * of subjects - instead of paying a create()/free() cycle per set. */
int sundialr_cvode_set_udata(void* m, void* udata);

/* --- Integration ----------------------------------------------------------- */

/* Set the state to y0 (length neq) and (re)start the integrator at t0. The first
 * call also performs the one-time CVodeInit and attaches the linear solver,
 * tolerances and Jacobian; later calls are CVodeReInit and keep those settings.
 * Call whenever the state jumps (a dose/event); a plain trajectory needs it once. */
int sundialr_cvode_reinit(void* m, double t0, const double* y0);

/* Advance to tout and copy the state into y (length neq). If treached is not
 * NULL it receives the time actually reached (== tout on success). Returns
 * CVODE's own code: 0 on success, < 0 on failure (state left at the last
 * successful point). */
int sundialr_cvode_solve(void* m, double tout, double* y, double* treached);

/* --- Introspection --------------------------------------------------------- */

/* Total internal steps taken since create() or the last reset_stats(),
 * accumulated ACROSS reinit() calls (CVodeReInit zeroes CVODE's own counter, so
 * the handle carries the running sum over segments). Returns -1 before the
 * first reinit. */
long sundialr_cvode_get_num_steps(void* m);

/* Zero the accumulated step count, so get_num_steps() counts from here - e.g.
 * at a subject boundary, for per-subject step counts from a shared handle.
 * May be called at any point, including mid-segment. */
int sundialr_cvode_reset_stats(void* m);

/* The message recorded for the most recent failure, or NULL if none. The string
 * is owned by the handle and valid until the next entry point or free(). */
const char* sundialr_cvode_last_err(void* m);

/* This build's ABI version (SUNDIALR_ABI_VERSION). */
int sundialr_abi_version(void);

#ifdef __cplusplus
}
#endif

#endif /* SUNDIALR_CAPI_H */
