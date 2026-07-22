//   Copyright (c) 2016-2026, Satyaprakash Nayak
//   Distributed under the BSD-3 licence; see the header of sundialr_capi.h.

// Registration of the sundialr C API as R "CCallable" entry points, so another
// package (e.g. mrgsolve) can reach them with R_GetCCallable("sundialr", name)
// after LinkingTo/Imports sundialr. Rcpp generates R_init_sundialr itself and
// would overwrite anything added there, so per the standard workaround this runs
// from an [[Rcpp::export]] function invoked by the package's .onLoad() (see
// R/zzz.R). The R-side wrapper is named ".register_capi" - the leading dot keeps
// it out of the package namespace (exportPattern only exports alpha-initial
// names), so it is not user-visible and needs no documentation.

#include <Rcpp.h>
#include <R_ext/Rdynload.h>

#include <sundialr_capi.h>

// [[Rcpp::export(".register_capi")]]
void register_capi() {
  R_RegisterCCallable("sundialr", "sundialr_cvode_create",     (DL_FUNC) sundialr_cvode_create);
  R_RegisterCCallable("sundialr", "sundialr_cvode_free",       (DL_FUNC) sundialr_cvode_free);
  R_RegisterCCallable("sundialr", "sundialr_cvode_set_rhs",    (DL_FUNC) sundialr_cvode_set_rhs);
  R_RegisterCCallable("sundialr", "sundialr_cvode_set_jac",    (DL_FUNC) sundialr_cvode_set_jac);
  R_RegisterCCallable("sundialr", "sundialr_cvode_set_tol_scalar", (DL_FUNC) sundialr_cvode_set_tol_scalar);
  R_RegisterCCallable("sundialr", "sundialr_cvode_set_tol_vector", (DL_FUNC) sundialr_cvode_set_tol_vector);
  R_RegisterCCallable("sundialr", "sundialr_cvode_set_max_steps",  (DL_FUNC) sundialr_cvode_set_max_steps);
  R_RegisterCCallable("sundialr", "sundialr_cvode_set_max_step",   (DL_FUNC) sundialr_cvode_set_max_step);
  R_RegisterCCallable("sundialr", "sundialr_cvode_set_min_step",   (DL_FUNC) sundialr_cvode_set_min_step);
  R_RegisterCCallable("sundialr", "sundialr_cvode_reinit",     (DL_FUNC) sundialr_cvode_reinit);
  R_RegisterCCallable("sundialr", "sundialr_cvode_solve",      (DL_FUNC) sundialr_cvode_solve);
  R_RegisterCCallable("sundialr", "sundialr_cvode_get_num_steps", (DL_FUNC) sundialr_cvode_get_num_steps);
  R_RegisterCCallable("sundialr", "sundialr_cvode_last_err",   (DL_FUNC) sundialr_cvode_last_err);
  R_RegisterCCallable("sundialr", "sundialr_abi_version",      (DL_FUNC) sundialr_abi_version);
}
