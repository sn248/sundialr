#ifndef JAC_FUNC_H
#define JAC_FUNC_H

// Prerequisites: Rcpp.h, nvector_serial.h, sunmatrix_dense.h

// Used by cvode, cvodes, cvsolve.
// R function signature: f(t, y, p)  ->  n-by-n matrix of d(ydot_i)/d(y_j)
static inline int jac_eval(sunrealtype t, N_Vector y, SUNMatrix JAC,
                           SEXP jac_eqn, Rcpp::NumericVector params) {
  int n = NV_LENGTH_S(y);
  Rcpp::NumericVector y1(n);
  sunrealtype *y_ptr = N_VGetArrayPointer(y);
  for (int i = 0; i < n; i++) y1[i] = y_ptr[i];

  Rcpp::Function jac_fun(jac_eqn);
  Rcpp::NumericMatrix J = jac_fun(t, y1, params);

  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++)
      SM_ELEMENT_D(JAC, i, j) = J(i, j);
  return 0;
}

// Used by ida only.
// R function signature: f(t, y, ydot, cj, p)  ->  n-by-n matrix of dF/dy + cj * dF/dydot
static inline int jac_eval_ida(sunrealtype t, sunrealtype cj,
                               N_Vector y, N_Vector yp, SUNMatrix JAC,
                               SEXP jac_eqn, Rcpp::NumericVector params) {
  int n = NV_LENGTH_S(y);
  Rcpp::NumericVector y1(n), yp1(n);
  sunrealtype *y_ptr  = N_VGetArrayPointer(y);
  sunrealtype *yp_ptr = N_VGetArrayPointer(yp);
  for (int i = 0; i < n; i++) { y1[i] = y_ptr[i]; yp1[i] = yp_ptr[i]; }

  Rcpp::Function jac_fun(jac_eqn);
  Rcpp::NumericMatrix J = jac_fun(t, y1, yp1, cj, params);

  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++)
      SM_ELEMENT_D(JAC, i, j) = J(i, j);
  return 0;
}

#endif
