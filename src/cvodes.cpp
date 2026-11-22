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

#include <Rcpp.h>
#include <algorithm>                   // to convert SensType to upper case - input cleaning
#include <cmath>                       // std::floor, std::fabs for sensIC
#include <cstring>                     // std::strcmp for sensIC names
#include <vector>

#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sundials/sundials_math.h>    /* definition of ABS */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h>

#include <check_retval.h>
#include <jac_func.h>
#include <sundials_scope_guard.h>
// CRAN fix: replace SUNDIALS' default abort()-based error handler with one that
// records the error for the solver to raise via stop() (see the header)
#include <sundials_err_handler.h>


// #define Ith(v,i)    NV_Ith_S(v,i-1)         /* i-th vector component i=1..NEQ */

using namespace Rcpp;

// struct to use if R or Rcpp function is input as RHS function----------------
struct rhs_func_sens{
  Function rhs_eqn;
  NumericVector params;
  double rtol;
  NumericVector atol;
  SEXP jac_eqn;
  SEXP sens_eqn;              // user sensitivity RHS, or R_NilValue
  sundials_err_record *err;   // collects errors raised inside the callbacks
  SEXP y_names;               // names(IC), set on y, ydot and yS, or R_NilValue
  // The array CVODES perturbs when it approximates the sensitivity equations
  // by difference quotients (the p argument of CVodeSetSensParams). Without
  // IC sensitivities this is params.begin() itself. With them it is a longer
  // array holding the parameters followed by one placeholder per IC
  // sensitivity, and its first params.length() entries are copied into
  // params before every call to R; see sync_params().
  const double *p_sens;
  // Scaling factors, one per sensitivity (the pbar of CVodeSetSensParams).
  // Read by ic_sens_rhs() to size its difference-quotient increment.
  const double *pbar;
};

// Make the parameter vector handed to R reflect any perturbation CVODES has
// applied to p_sens. A no-op unless IC sensitivities are requested, because
// p_sens then aliases params.
static inline void sync_params(struct rhs_func_sens *data) {
  double *p = data->params.begin();
  if (data->p_sens != p) {
    std::copy(data->p_sens, data->p_sens + data->params.length(), p);
  }
}

// function called by CVodeInit if user inputs R function
// Called by SUNDIALS from its own C code, so the body runs under
// sundials_callback_guard; see the note on rhs_function in rhs_func.cpp.
int rhs_function_sens(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data){

  // cast void pointer to pointer to struct and assign rhs to a Function
  struct rhs_func_sens *my_rhs_fun = (struct rhs_func_sens*)user_data;

  // nothing to record against, so just report the failure to SUNDIALS
  if(!my_rhs_fun){ return(-1); }

  return sundials_callback_guard(my_rhs_fun->err, [&]() -> int {

    // convert y (N_Vector) to NumericVector y1
    int y_len = NV_LENGTH_S(y);

    NumericVector y1(y_len);    // filled with zeros
    sunrealtype *y_ptr = N_VGetArrayPointer(y);
    for (int i = 0; i < y_len; i++){
      y1[i] = y_ptr[i];
    }
    set_state_names(y1, my_rhs_fun->y_names);
    sync_params(my_rhs_fun);

    NumericVector ydot1(y_len);    // filled with zeros

    Function rhs_fun = (*my_rhs_fun).rhs_eqn;                // function
    NumericVector p_values = (*my_rhs_fun).params;           // rate parameters

    if (!rhs_fun || (TYPEOF(rhs_fun) != CLOSXP)){
      stop("The RHS function is not an R or Rcpp function, stopping!");
    }

    // use the function to calculate value of RHS ----
    ydot1 = rhs_fun(t, y1, p_values);

    // guards the element-by-element copy below against a short return
    if (ydot1.length() != y_len){
      stop("The RHS function must return a vector of the same length as the state vector: expected %d, got %d",
           y_len, ydot1.length());
    }

    // convert NumericVector ydot1 to N_Vector ydot
    sunrealtype *ydot_ptr = N_VGetArrayPointer(ydot);
    for (int i = 0; i<  y_len; i++){
      ydot_ptr[i] = ydot1[i];
    }

    // If everything went smoothly, return 0
    return(0);
  });
}

// EwtSet Function. Computes the error weights at the current solution.
int ewt(N_Vector y, N_Vector w, void *user_data)
{
  sunrealtype yy, ww, rtol;

  struct rhs_func_sens *my_rhs_fun = (struct rhs_func_sens*)user_data;
  if(!my_rhs_fun){ return(-1); }

  rtol = my_rhs_fun->rtol;

  // Read the tolerances through a pointer rather than copying the vector.
  // CVODES calls this on every internal step, far more often than the RHS,
  // so a copy here is on the hottest path in the solve.
  const double *atol = my_rhs_fun->atol.begin();

  sunindextype y_len = N_VGetLength_Serial(y);

  for (int i=1; i<=y_len; i++) {

    yy = NV_Ith_S(y,i-1);
    ww = rtol * SUNRabs(yy) + atol[i-1];
    if (ww <= 0.0) return (-1);
    NV_Ith_S(w,i-1) = 1.0/ww;

    // Rcpp::Rcout << NV_Ith_S(w,i-1) << "\n";
  }

  return(0);
}

// Call the jacobian function if manual jacobian is provided
static int jac_cvodes(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix JAC,
                      void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    struct rhs_func_sens *data = (struct rhs_func_sens*)user_data;
    if (!data) { return -1; }
    return sundials_callback_guard(data->err, [&]() -> int {
      sync_params(data);
      return jac_eval(t, y, JAC, data->jac_eqn, data->params, data->y_names);
    });
}

// Calls the user's RHS at (t, y) and checks the length of what comes back.
// Used by ic_sens_rhs(); the caller has already run sync_params().
static NumericVector call_rhs(struct rhs_func_sens *data, sunrealtype t,
                              const NumericVector &y) {
  NumericVector f = data->rhs_eqn(t, y, data->params);
  if (f.length() != y.length()) {
    stop("The RHS function must return a vector of the same length as the state vector: expected %d, got %d",
         (int)y.length(), (int)f.length());
  }
  return f;
}

// d(yS)/dt for a sensitivity with respect to an initial condition, which is
// J %*% yS: an initial value enters no right-hand side, so there is no df/dp
// term. sundialr computes this itself rather than asking the user's
// `sensitivity` function, so that function only ever sees parameter indices
// and one written for parameters alone keeps working when sensIC is added.
//
// With a user Jacobian the product is exact. Without one it is the centred
// directional difference (f(y + d yS) - f(y - d yS)) / (2 d), the scheme
// CVODES' own difference quotients use by default (CV_CENTERED), with the
// increment chosen as CVODES chooses it in cvSensRhs1InternalDQ:
//   d = 1 / max(||yS||_WRMS, 1 / (sqrt(max(rtol, uround)) * pbar))
// The weights are those of ewt() below, inlined because it is defined later.
static void ic_sens_rhs(struct rhs_func_sens *data, sunrealtype t,
                        const NumericVector &y1, const sunrealtype *yS_ptr,
                        int iS, sunrealtype *ySdot_ptr) {
  int n = y1.length();

  if (data->jac_eqn != R_NilValue) {
    Function jac_fun(data->jac_eqn);
    NumericMatrix J = jac_fun(t, y1, data->params);
    if (J.nrow() != n || J.ncol() != n) {
      stop("The Jacobian function must return a %d-by-%d matrix; got %d-by-%d",
           n, n, J.nrow(), J.ncol());
    }
    for (int i = 0; i < n; i++) {
      double s = 0.0;
      for (int j = 0; j < n; j++) s += J(i, j) * yS_ptr[j];
      ySdot_ptr[i] = s;
    }
    return;
  }

  const double *atol = data->atol.begin();
  double wrms = 0.0;
  for (int i = 0; i < n; i++) {
    double w = 1.0 / (data->rtol * std::fabs(y1[i]) + atol[i]);
    wrms += (yS_ptr[i] * w) * (yS_ptr[i] * w);
  }
  wrms = std::sqrt(wrms / n);
  double Delta  = std::sqrt(std::max(data->rtol, (double)SUN_UNIT_ROUNDOFF));
  double rdelta = 1.0 / (Delta * data->pbar[iS]);
  double d      = 1.0 / std::max(wrms, rdelta);

  NumericVector yp(n), ym(n);
  for (int i = 0; i < n; i++) {
    yp[i] = y1[i] + d * yS_ptr[i];
    ym[i] = y1[i] - d * yS_ptr[i];
  }
  set_state_names(yp, data->y_names);
  set_state_names(ym, data->y_names);
  NumericVector fp = call_rhs(data, t, yp);
  NumericVector fm = call_rhs(data, t, ym);
  for (int i = 0; i < n; i++) ySdot_ptr[i] = (fp[i] - fm[i]) / (2.0 * d);
}

// Sensitivity right-hand side, used when the caller supplies `sensitivity`.
// CVODES calls this once per sensitivity (a CVSensRhs1Fn). For a parameter it
// delegates to the R function stored in data->sens_eqn, whose signature is
//   sens_rhs(t, y, ydot, iS, yS, p)  ->  numeric vector of length(y)
// returning d(yS_iS)/dt = J %*% yS_iS + df/dp_iS. iS is passed 1-based to R.
// The sensitivities past the parameters are the sensIC ones, which never reach
// the R function: ic_sens_rhs() computes them.
// Called by SUNDIALS from its own C code, so the body runs under
// sundials_callback_guard; see the note on rhs_function in rhs_func.cpp.
static int sens_rhs1_cvodes(int Ns, sunrealtype t, N_Vector y, N_Vector ydot,
                            int iS, N_Vector yS, N_Vector ySdot,
                            void *user_data, N_Vector tmp1, N_Vector tmp2) {
    struct rhs_func_sens *data = (struct rhs_func_sens*)user_data;
    if (!data) { return -1; }
    return sundials_callback_guard(data->err, [&]() -> int {

      int n = NV_LENGTH_S(y);

      if (iS >= data->params.length()) {
        NumericVector y1(n);
        sunrealtype *y_ptr = N_VGetArrayPointer(y);
        for (int i = 0; i < n; i++) y1[i] = y_ptr[i];
        set_state_names(y1, data->y_names);
        sync_params(data);
        ic_sens_rhs(data, t, y1, N_VGetArrayPointer(yS), iS,
                    N_VGetArrayPointer(ySdot));
        return 0;
      }

      NumericVector y1(n), ydot1(n), yS1(n);
      sunrealtype *y_ptr    = N_VGetArrayPointer(y);
      sunrealtype *ydot_ptr = N_VGetArrayPointer(ydot);
      sunrealtype *yS_ptr   = N_VGetArrayPointer(yS);
      for (int i = 0; i < n; i++) { y1[i] = y_ptr[i]; ydot1[i] = ydot_ptr[i]; yS1[i] = yS_ptr[i]; }
      set_state_names(y1, data->y_names);
      set_state_names(ydot1, data->y_names);
      set_state_names(yS1, data->y_names);
      sync_params(data);

      Function sens_fun(data->sens_eqn);
      // iS arrives 0-based from CVODES; hand R the 1-based parameter index
      NumericVector ySdot1 = sens_fun(t, y1, ydot1, iS + 1, yS1, data->params);

      // guards the element-by-element copy below against a short return
      if (ySdot1.length() != n) {
        stop("The sensitivity function must return a vector of the same length as the state vector: expected %d, got %d",
             n, ySdot1.length());
      }

      sunrealtype *ySdot_ptr = N_VGetArrayPointer(ySdot);
      for (int i = 0; i < n; i++) ySdot_ptr[i] = ySdot1[i];

      return 0;
    });
}
//------------------------------------------------------------------------------
//' cvodes
//'
//' CVODES solver to solve ODEs and calculate sensitivities
//'@param time_vector time vector
//'@param IC Initial Conditions. If named, the names are passed on to every function the solver calls: \code{y}, \code{ydot} and \code{yS} arrive named, so they can be indexed by name (\code{y["A"]})
//'@param input_function Right Hand Side function of ODEs
//'@param Parameters Parameters input to ODEs. Their names, if any, are kept
//'@param reltolerance Relative Tolerance (a scalar, default value  = 1e-04)
//'@param abstolerance Absolute Tolerance (a scalar or vector with length equal to ydot, default = 1e-04)
//'@param SensType Sensitivity Type - allowed values are "STG" (for Staggered, default) or "SIM" (for Simultaneous)
//'@param ErrCon Error Control - allowed values are TRUE or FALSE (default). When TRUE the sensitivities are included in the local error test, which makes them more accurate at the cost of more steps. With the difference-quotient sensitivities (no \code{sensitivity} function) very tight tolerances, near \code{reltolerance = 1e-10}, may then exhaust the step limit; loosen them or supply \code{sensitivity}
//'@param jacobian (Optional) Jacobian of the RHS with signature \code{function(t, y, p)}. Default is NULL
//'@param sensitivity (Optional) Sensitivity right-hand side with signature \code{function(t, y, ydot, iS, yS, p)} returning the derivative \code{d(yS_iS)/dt = J \%*\% yS_iS + df/dp_iS} as a numeric vector of \code{length(y)}, where \code{iS} is the 1-based parameter index. The function is only ever called for parameters, so \code{iS} never exceeds \code{length(Parameters)}, with or without \code{sensIC}: the sensitivities with respect to initial conditions satisfy \code{d(yS)/dt = J \%*\% yS}, which the solver computes itself, exactly from \code{jacobian} when one is supplied and otherwise by a centred difference of \code{input_function} along \code{yS}. Supply \code{jacobian} as well when the initial-condition sensitivities need to be exact. Default is NULL, in which case the sensitivity equations are approximated by finite differences of the RHS
//'@param sensIC (Optional) States whose initial values are also sensitivity parameters, as 1-based indices into \code{IC} or as names of \code{IC}. Each adds one sensitivity, started at the unit vector for its state, so its block holds \code{d(y)/d(IC[i])}. Default is NULL: sensitivities with respect to \code{Parameters} only
//'@param Events (Optional) Discontinuities in the solution, such as bolus doses, as a data frame with exactly three columns (names ignored): the 1-based index of the state, the time of the event, and the value to add to that state at that time. This is the format of \code{Events} in \code{\link{cvsolve}}, and the value is always added, including at the initial time, where it adds to \code{IC}. The solver integrates up to each event, applies it and restarts, carrying the sensitivities across unchanged: the event amounts and times are fixed, not sensitivity parameters, so an event moves the states but not their derivatives. Unlike \code{cvsolve}, the result still has one row per element of \code{time_vector}; an output time that coincides with an event reports the state just after the event. Default is NULL
//'@returns A matrix with \code{1 + length(IC) * (length(Parameters) + length(sensIC))} columns and one row per element of \code{time_vector}. The first column is the time vector. The remaining columns hold the sensitivities in blocks of \code{length(IC)}, one block per parameter: the sensitivities of all states with respect to the first parameter, then of all states with respect to the second parameter, and so on. The sensitivity of state \code{j} with respect to parameter \code{k} is in column \code{1 + (k - 1) * length(IC) + j}. The blocks for \code{sensIC} follow the parameter blocks, in the order given in \code{sensIC}, so they are numbered as if they were parameters \code{length(Parameters) + 1}, \code{length(Parameters) + 2}, and so on.
//'
//' The state solution is attached as the attribute \code{"solution"}: a matrix with the same layout as the output of \code{\link{cvode}}, time in the first column followed by one column per state. Retrieve it with \code{attr(result, "solution")}. As with any attribute, subsetting the returned matrix drops it, so extract it first.
//'@example /inst/examples/cvs_Roberts_dns.r
// [[Rcpp::export]]
NumericMatrix cvodes(NumericVector time_vector, NumericVector IC,
                      SEXP input_function,
                      NumericVector Parameters,
                      double reltolerance = 0.0001,
                      NumericVector abstolerance = 0.0001,
                      std::string SensType = "STG",
                      bool ErrCon = false,
                      Nullable<Function> jacobian = R_NilValue,
                      Nullable<Function> sensitivity = R_NilValue,
                      SEXP sensIC = R_NilValue,
                      Nullable<DataFrame> Events = R_NilValue){

  int flag;

  int time_vec_len = time_vector.length();
  double time;
  int NOUT = time_vec_len;
  sunrealtype T0 = SUN_RCONST(time_vector[0]);     // Initial Time

  // Initial Conditions
  int y_len = IC.length();

  // Relative Tolerance
  sunrealtype reltol = reltolerance;

  // absolute tolerance is either length == 1 or equal to length of IC
  // If abstol is not equal to 1 or abstol is not equal to IC, then stop
  int abstol_len = abstolerance.length();
  if(abstol_len != 1 && abstol_len != y_len){
     stop("Absolute tolerance must be a scalar or a vector of same length as IC \n");
  }
  // Create an abstolerance vector equal to IC length
  Rcpp::NumericVector abstol(y_len);
  if(abstol_len == 1){
    for (int i = 0; i<y_len; i++){
      abstol[i] = abstolerance[0];
    }
  }
  if (abstol_len == y_len){
    for (int i = 0; i<y_len; i++){
      abstol[i] = abstolerance[i];
    }
  }

  // Number of model parameters
  int NP = Parameters.length();

  // States whose initial values are sensitivity parameters (0-based), from
  // sensIC. Validated here, before any SUNDIALS object exists.
  std::vector<int> ic_idx;
  if (sensIC != R_NilValue && Rf_length(sensIC) > 0) {
    int n_ic = Rf_length(sensIC);
    if (TYPEOF(sensIC) == STRSXP) {
      SEXP ic_names = Rf_getAttrib(IC, R_NamesSymbol);
      if (ic_names == R_NilValue) {
        stop("sensIC is given by name, but IC has no names");
      }
      for (int m = 0; m < n_ic; m++) {
        SEXP s = STRING_ELT(sensIC, m);
        if (s == NA_STRING) stop("sensIC must not contain NA");
        int found = -1;
        for (int j = 0; j < y_len; j++) {
          SEXP nj = STRING_ELT(ic_names, j);
          if (nj != NA_STRING && std::strcmp(CHAR(nj), CHAR(s)) == 0) { found = j; break; }
        }
        if (found < 0) stop("sensIC names a state that is not in names(IC): %s", CHAR(s));
        ic_idx.push_back(found);
      }
    } else if (TYPEOF(sensIC) == INTSXP || TYPEOF(sensIC) == REALSXP) {
      NumericVector idx(sensIC);
      for (int m = 0; m < n_ic; m++) {
        double v = idx[m];
        if (ISNAN(v) || v != std::floor(v) || v < 1 || v > y_len) {
          stop("sensIC must hold whole numbers between 1 and length(IC) = %d", y_len);
        }
        ic_idx.push_back((int)v - 1);
      }
    } else {
      stop("sensIC must be NULL, a vector of state indices or a vector of state names");
    }
    for (size_t a = 0; a < ic_idx.size(); a++)
      for (size_t b = a + 1; b < ic_idx.size(); b++)
        if (ic_idx[a] == ic_idx[b]) stop("sensIC must not list a state more than once");
  }

  // Events, validated here, before any SUNDIALS object exists, by the same
  // rules as cvsolve(). The state index is used to subscript the state vector
  // directly, so anything but a whole number in 1..y_len would read or write
  // outside it. Held as records sorted by time; a stable sort keeps events at
  // one time in the order given, though they commute since all of them add.
  struct sens_event { double t; int state; double value; };
  std::vector<sens_event> events;
  if (Events.isNotNull()) {
    DataFrame ev(Events);
    if (ev.size() != 3) {
      stop("Events must be a data frame with exactly three columns: state index, time and value");
    }
    NumericVector ev_state = ev[0];
    NumericVector ev_time  = ev[1];
    NumericVector ev_value = ev[2];
    double t_max = Rcpp::max(time_vector);
    for (int i = 0; i < ev_state.length(); i++) {
      double s = ev_state[i], t = ev_time[i], v = ev_value[i];
      if (ISNAN(s) || s < 1 || s > y_len || s != std::floor(s)) {
        stop("The state index in the first column of the Events dataframe must be a whole number between 1 and the number of states");
      }
      if (ISNAN(t) || t < time_vector[0] || t > t_max) {
        stop("The event time in the second column of the Events dataframe must lie between the first and last time points");
      }
      if (!R_FINITE(v)) {
        stop("The event value in the third column of the Events dataframe must be a finite number");
      }
      events.push_back({t, (int)s - 1, v});
    }
    std::stable_sort(events.begin(), events.end(),
                     [](const sens_event &a, const sens_event &b) { return a.t < b.t; });
  }

  // Total number of sensitivities: one per parameter, then one per sensIC
  // state, in that order. Needed by the guard to size the yS free.
  int NIC = (int)ic_idx.size();
  int NS  = NP + NIC;

  // The parameter array CVODES perturbs for its difference quotients, and the
  // scaling factors pbar that go with it. Without IC sensitivities this is
  // Parameters itself, exactly as before. With them, CVODES needs one entry
  // per sensitivity, so the parameters are followed by one placeholder per
  // sensIC state. The RHS never reads a placeholder, so perturbing it changes
  // nothing: the difference quotient then yields J %*% yS with no df/dp term,
  // which is precisely the IC sensitivity equation. Because CVODES perturbs
  // this array rather than Parameters, the R function is handed a separate
  // length-NP vector that sync_params() refreshes from it before each call.
  NumericVector p_cvodes, pbar, p_for_R;
  if (NIC == 0) {
    p_cvodes = Parameters;
    pbar     = Parameters;
    p_for_R  = Parameters;
  } else {
    p_cvodes = NumericVector(NS);
    pbar     = NumericVector(NS);
    p_for_R  = NumericVector(NP);
    for (int k = 0; k < NP; k++) {
      p_cvodes[k] = Parameters[k];
      pbar[k]     = Parameters[k];
      p_for_R[k]  = Parameters[k];
    }
    if (NP > 0) p_for_R.attr("names") = Rf_getAttrib(Parameters, R_NamesSymbol);
    for (int m = 0; m < NIC; m++) {
      double ic = IC[ic_idx[m]];
      p_cvodes[NP + m] = ic;
      // pbar scales the sensitivity's tolerances (atol / pbar under
      // CVodeSensEEtolerances) and the difference-quotient increment. An IC
      // sensitivity d(y)/d(y0_i) has the size of y / y0_i, which is the same
      // convention as pbar = p for a parameter; an IC of zero falls back to 1.
      pbar[NP + m] = (ic != 0.0) ? std::fabs(ic) : 1.0;
    }
  }

  // Receives SUNDIALS errors. Declared before the guard so that it is
  // destroyed after it - the SUNContext freed there holds a pointer to it.
  sundials_err_record sun_err;

  // SUNDIALS objects, released by the guard below on every exit path
  SUNContext sunctx      = NULL;
  void *cvode_mem        = NULL;
  N_Vector y0            = NULL;
  N_Vector *yS           = NULL;
  SUNMatrix SM           = NULL;
  SUNLinearSolver LS     = NULL;

  // Free in the same order the trailing free calls used to; runs on the
  // normal return and on the exceptions thrown by stop() below. yS must be
  // freed with the same count it was cloned with (NS).
  auto sundials_cleanup = make_scope_guard([&]{
    if (y0)        N_VDestroy(y0);
    if (yS)        N_VDestroyVectorArray(yS, NS);
    if (cvode_mem) CVodeFree(&cvode_mem);
    if (LS)        SUNLinSolFree(LS);
    if (SM)        SUNMatDestroy(SM);
    if (sunctx)    SUNContext_Free(&sunctx);
  });

  // Set context for sundials
  SUNContext_Create(SUN_COMM_NULL, &sunctx);
  // CRAN fix: redirect SUNDIALS fatal errors to R instead of calling abort()
  SUNContext_PushErrHandler(sunctx, sundials_r_err_handler, &sun_err);
  sundials_check(sun_err);   // context creation is not otherwise checked

  // Set the initial conditions-------------------------------------------------
  y0 = N_VNew_Serial(y_len, sunctx);
  sundials_check(sun_err);   // vector allocations are not otherwise checked
  sunrealtype *y0_ptr = N_VGetArrayPointer(y0);
  for (int i = 0; i<y_len; i++){
    y0_ptr[i] = IC[i];
  }
  // Events at the initial time add to the initial conditions, as in cvsolve().
  // ev_next is the first event still to be applied in the solve loop. The
  // sensitivities start from yS(0) as usual: the added amount is a constant,
  // so d(y0)/dp is unchanged and d(y0)/d(IC[i]) is still the unit vector.
  size_t ev_next = 0;
  while (ev_next < events.size() && events[ev_next].t == T0) {
    y0_ptr[events[ev_next].state] += events[ev_next].value;
    ev_next++;
  }

  // Check SensType input ---------------------------------------------------------
  // convert the SensType argument to upper case first
  std::transform(SensType.begin(), SensType.end(), SensType.begin(), ::toupper);
  if(SensType.compare("STG") != 0 && SensType.compare("SIM") != 0){
    stop("SensType argument can be STG (for Staggered) or SIM (for simulated) \nsee cvode help for example usage");
  }

  //-------------------------------------------------------------------------------

  // Check Error Control input ----------------------------------------------------
  if(ErrCon != TRUE && ErrCon != FALSE){
    stop("ErrCon can only be either TRUE or FALSE");
  }

  // //----------------------------------------------------------------------------
  /* Set sensitivity initial conditions */
  yS = N_VCloneVectorArray(NS, y0);   // vector of sensitivities (yS)

  // Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval(cvode_mem, "CVodeCreate")) { sundials_stop(sun_err, "CVodeCreate", "Stopping cvodes, cannot allocate memory for CVODES!"); }

  //-- assign user input to the struct based on SEXP type of input_function
  if (!input_function){
    stop("Something is wrong with input function, stopping!");
  }

  // Rcout << TYPEOF(input_function) << "\n";

  if (!input_function){ stop("There is no input function, stopping!"); }

  if(TYPEOF(input_function) != CLOSXP) { stop("Incorrect input function type - input function can be an R or Rcpp function"); }

  // realtype *params = Parameters.begin();
  // Initialize the struct for user data
  // Jacobian and sensitivity RHS have initial value of NULL
  SEXP jac_sexp = R_NilValue;
  if (jacobian.isNotNull()) jac_sexp = as<SEXP>(jacobian);
  SEXP sens_sexp = R_NilValue;
  if (sensitivity.isNotNull()) sens_sexp = as<SEXP>(sensitivity);
  struct rhs_func_sens my_rhs_function = {input_function,
                                          p_for_R,
                                          reltol,
                                          abstol,
                                          jac_sexp,
                                          sens_sexp,
                                          &sun_err,
                                          Rf_getAttrib(IC, R_NamesSymbol),
                                          p_cvodes.begin(),
                                          pbar.begin()};

  // setting the user_data in rhs function
  flag = CVodeSetUserData(cvode_mem, (void*)&my_rhs_function);
  if (check_retval(flag, "CVodeSetUserData")) { sundials_stop(sun_err, "CVodeSetUserData", "Stopping cvodes, something went wrong in setting user data!"); }

  /* Allocate space for CVODES */
  flag = CVodeInit(cvode_mem, rhs_function_sens, T0, y0);
  if (check_retval(flag, "CVodeInit")) { sundials_stop(sun_err, "CVodeInit", "Stopping cvodes, something went wrong in allocating space for CVODES!"); }

  /* Use private function to compute error weights */
  flag = CVodeWFtolerances(cvode_mem, ewt);
  if (check_retval(flag, "CVodeWFtolerances")) { sundials_stop(sun_err, "CVodeWFtolerances", "Stopping cvodes, something went wrong in computing error weights!"); }

  /* Create dense SUNMatrix */
  sunindextype y_len_M = y_len;
  SM = SUNDenseMatrix(y_len_M, y_len_M, sunctx);
  if (check_retval(SM, "SUNDenseMatrix")) { sundials_stop(sun_err, "SUNDenseMatrix", "Stopping cvodes, something went wrong in setting SUNDenseMatrix!"); }

  /* Create dense SUNLinearSolver */
  LS = SUNLinSol_Dense(y0, SM, sunctx);
  if (check_retval(LS, "SUNLinSol_Dense")) { sundials_stop(sun_err, "SUNLinSol_Dense", "Stopping cvodes, something went wrong in setting Linear Solver!"); }

  /* Attach the matrix and linear solver */
  flag = CVodeSetLinearSolver(cvode_mem, LS, SM);
  if (check_retval(flag, "CVodeSetLinearSolver")) { sundials_stop(sun_err, "CVodeSetLinearSolver", "Stopping cvodes, something went wrong in attaching SUNDenseMatrix and Linear Solver!"); }

  /* If manual Jacobian is provided, use the jacobian */
  if (jacobian.isNotNull()){
    flag = CVodeSetJacFn(cvode_mem, jac_cvodes);
    if(check_retval(flag, "CVodeSetJacFn")) { sundials_stop(sun_err, "CVodeSetJacFn", "Stopping cvodes, something went wrong in setting the Jacobian function!"); }
  }

  if (check_retval(yS, "N_VCloneVectorArray")) { sundials_stop(sun_err, "N_VCloneVectorArray", "Stopping cvodes, something went wrong in setting Sensitivity Array!"); }
  for (int is=0;is<NS;is++) N_VConst(SUN_RCONST(0.0), yS[is]);
  // An IC sensitivity starts at the unit vector of its state: d(y0)/d(y0_i) = e_i
  for (int m=0;m<NIC;m++) NV_Ith_S(yS[NP + m], ic_idx[m]) = SUN_RCONST(1.0);

  /* Call CVodeSensInit1 to activate forward sensitivity computations
   and allocate internal memory for CVODES related to sensitivity
   calculations. Computes the right-hand sides of the sensitivity
   ODE, one at a time */
  //
  // The fourth argument is the sensitivity right-hand side. When the caller
  // supplies `sensitivity`, sens_rhs1_cvodes is passed and CVODES uses that R
  // function to evaluate the sensitivity equations analytically. Passing NULL
  // (the default) instead tells CVODES to approximate it by finite differences
  // of the state right-hand side, which needs no extra function from the caller
  // but evaluates the state RHS an extra NS times per internal step - each an R
  // callback, the most expensive thing in the solve - and is only as accurate
  // as the difference quotients.
  //
  // (int) flag to set sensitivity solution method - see CVodeSensInit1
  int ism = CV_STAGGERED;
  if (SensType.compare("SIM") == 0) ism = CV_SIMULTANEOUS;
  CVSensRhs1Fn fS1 = sensitivity.isNotNull() ? sens_rhs1_cvodes : NULL;
  flag = CVodeSensInit1(cvode_mem, NS, ism, fS1, yS);
  if(check_retval(flag, "CVodeSensInit1")) { sundials_stop(sun_err, "CVodeSensInit1", "Stopping cvodes, something went wrong in calculating Sensitivities!"); }

  /* Call CVodeSensEEtolerances to estimate tolerances for sensitivity
   variables based on the rolerances supplied for states variables and
   the scaling factor pbar */
  flag = CVodeSensEEtolerances(cvode_mem);
  if(check_retval(flag, "CVodeSensEEtolerances")) { sundials_stop(sun_err, "CVodeSensEEtolerances", "Stopping cvodes, something went wrong in estimating tolerances for sensitivities!"); }

  /* Include the sensitivities in the local error test if ErrCon is TRUE.
   Before 0.2.2 ErrCon was validated but never passed on, so this was always
   CVODES' default of SUNFALSE whatever the caller asked for. The setting
   survives CVodeSensReInit, so it holds across Events too. */
  flag = CVodeSetSensErrCon(cvode_mem, ErrCon ? SUNTRUE : SUNFALSE);
  if(check_retval(flag, "CVodeSetSensErrCon")) { sundials_stop(sun_err, "CVodeSetSensErrCon", "Stopping cvodes, something went wrong in setting sensitivity error control!"); }

  /* Call CVodeSetSensParams to specify problem parameter information for
   sensitivity calculations */
  // struct rhs_func_sens *ptr = &my_rhs_function;        // struct UserData storing my data
  // p = (my_rhs_function.params).begin();
  // Rcout << (my_rhs_function.params).begin() << "\n";
  // Rcout << &(ptr->params[0]);
  flag = CVodeSetSensParams(cvode_mem, p_cvodes.begin(), pbar.begin(), NULL);
  if (check_retval(flag, "CVodeSetSensParams")) { sundials_stop(sun_err, "CVodeSetSensParams", "Stopping cvodes, something went wrong in setting Sensitivity Parameters!"); }

  // First row for initial conditions, First column is for time
  int y_len_1 = y_len + 1;
  NumericMatrix soln(Dimension(time_vec_len,y_len_1));

  // fill the first row of soln matrix with the initial state (IC plus any
  // events at the initial time)
  soln(0,0) = time_vector[0];   // get the first time value
  for(int i = 0; i<y_len; i++){
    soln(0,i+1) = y0_ptr[i];
  }

  NumericMatrix sens(Dimension(time_vec_len, (y_len * NS) + 1));
  // Store the Sensitivity Results. The first row is yS(0): zero for every
  // parameter, the unit vector of its state for every IC sensitivity.
  sens(0,0) = time_vector[0];           // get the first time value
  for (int m = 0; m < NIC; m++){
    sens(0, y_len*(NP + m) + ic_idx[m] + 1) = 1.0;
  }

  // Integrates from the current time to t, leaving the state in y0 and the
  // sensitivities in yS.
  auto advance_to = [&](sunrealtype t) {
    flag = CVode(cvode_mem, t, y0, &time, CV_NORMAL);
    if (check_retval(flag, "CVode")) { sundials_stop(sun_err, "CVode", "Stopping cvodes, something went wrong in solving the system using CVODE!"); } // Something went wrong in solving it!
    flag = CVodeGetSens(cvode_mem, &time, yS);
    if (check_retval(flag, "CVodeGetSens")) { sundials_stop(sun_err, "CVodeGetSens", "Stopping cvodes, something went wrong in calculating Sensitivities!"); }
  };

  // Writes the current state and sensitivities into output row `row`.
  auto store_row = [&](int row, sunrealtype t) {
    soln(row, 0) = t;                     // first column is for time
    for (int i = 0; i<y_len; i++){
      soln(row, i+1) = y0_ptr[i];
    }
    sens(row, 0) = t;                     // first column is for time
    for (int i = 0; i < NS; i++){
      sunrealtype *yS_ptr = N_VGetArrayPointer(yS[i]);   // sensitivities w.r.t. param i
      for(int j = 0; j < y_len; j++){
        sens(row, y_len*i+j+1) = yS_ptr[j]; // sensitivities of y1 ... yJ w.r.t param i
      }
    }
  };

  // t_start is the time CVODES was last (re)initialised at, and at_start says
  // no step has been taken since, so y0 and yS still hold the state there.
  // An output time equal to it is then filled from y0 and yS directly:
  // CVode() refuses an output time at its initial time.
  sunrealtype t_start = T0;
  bool at_start = true;

  sunrealtype tout;  // For output times
  for(int iout = 0; iout < NOUT-1; iout++) {

    // output times start from the index after initial time
    tout = time_vector[iout+1];

    // Apply every event up to and including tout. For each event time:
    // integrate to it, add all the events at that time to the state, and
    // restart both the state and the sensitivity integration from there.
    // The sensitivities are carried across unchanged - an event adds a fixed
    // amount, which depends on no parameter and no initial condition - but
    // CVODES must still be restarted, because its step-size and order history
    // describe a solution that no longer exists past the jump.
    while (ev_next < events.size() && events[ev_next].t <= tout) {
      // te is past the solver's current time: events at T0 were folded into
      // the IC, and every earlier event or output time has been consumed.
      sunrealtype te = events[ev_next].t;
      advance_to(te);
      while (ev_next < events.size() && events[ev_next].t == te) {
        y0_ptr[events[ev_next].state] += events[ev_next].value;
        ev_next++;
      }
      flag = CVodeReInit(cvode_mem, te, y0);
      if (check_retval(flag, "CVodeReInit")) { sundials_stop(sun_err, "CVodeReInit", "Stopping cvodes, something went wrong in reinitializing the ODE system!"); }
      flag = CVodeSensReInit(cvode_mem, ism, yS);
      if (check_retval(flag, "CVodeSensReInit")) { sundials_stop(sun_err, "CVodeSensReInit", "Stopping cvodes, something went wrong in reinitializing the sensitivities!"); }
      t_start  = te;
      at_start = true;
    }

    if (at_start && tout == t_start) {
      // At the (re)initialisation time: nothing to integrate. With an event
      // at tout this reports the state just after it.
      store_row(iout+1, tout);
    } else {
      advance_to(tout);
      at_start = false;
      store_row(iout+1, time);
    }

  }

  /* SUNDIALS objects are released by sundials_cleanup on scope exit */

  sens.attr("solution") = soln;
  return sens;



}

