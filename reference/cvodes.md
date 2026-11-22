# cvodes

CVODES solver to solve ODEs and calculate sensitivities

## Usage

``` r
cvodes(
  time_vector,
  IC,
  input_function,
  Parameters,
  reltolerance = 1e-04,
  abstolerance = 1e-04,
  SensType = "STG",
  ErrCon = FALSE,
  jacobian = NULL,
  sensitivity = NULL,
  sensIC = NULL,
  Events = NULL
)
```

## Arguments

- time_vector:

  time vector

- IC:

  Initial Conditions. If named, the names are passed on to every
  function the solver calls: `y`, `ydot` and `yS` arrive named, so they
  can be indexed by name (`y["A"]`)

- input_function:

  Right Hand Side function of ODEs

- Parameters:

  Parameters input to ODEs. Their names, if any, are kept

- reltolerance:

  Relative Tolerance (a scalar, default value = 1e-04)

- abstolerance:

  Absolute Tolerance (a scalar or vector with length equal to ydot,
  default = 1e-04)

- SensType:

  Sensitivity Type - allowed values are "STG" (for Staggered, default)
  or "SIM" (for Simultaneous)

- ErrCon:

  Error Control - allowed values are TRUE or FALSE (default). When TRUE
  the sensitivities are included in the local error test, which makes
  them more accurate at the cost of more steps. With the
  difference-quotient sensitivities (no `sensitivity` function) very
  tight tolerances, near `reltolerance = 1e-10`, may then exhaust the
  step limit; loosen them or supply `sensitivity`

- jacobian:

  (Optional) Jacobian of the RHS with signature `function(t, y, p)`.
  Default is NULL

- sensitivity:

  (Optional) Sensitivity right-hand side with signature
  `function(t, y, ydot, iS, yS, p)` returning the derivative
  `d(yS_iS)/dt = J %*% yS_iS + df/dp_iS` as a numeric vector of
  `length(y)`, where `iS` is the 1-based parameter index. The function
  is only ever called for parameters, so `iS` never exceeds
  `length(Parameters)`, with or without `sensIC`: the sensitivities with
  respect to initial conditions satisfy `d(yS)/dt = J %*% yS`, which the
  solver computes itself, exactly from `jacobian` when one is supplied
  and otherwise by a centred difference of `input_function` along `yS`.
  Supply `jacobian` as well when the initial-condition sensitivities
  need to be exact. Default is NULL, in which case the sensitivity
  equations are approximated by finite differences of the RHS

- sensIC:

  (Optional) States whose initial values are also sensitivity
  parameters, as 1-based indices into `IC` or as names of `IC`. Each
  adds one sensitivity, started at the unit vector for its state, so its
  block holds `d(y)/d(IC[i])`. Default is NULL: sensitivities with
  respect to `Parameters` only

- Events:

  (Optional) Discontinuities in the solution, such as bolus doses, as a
  data frame with exactly three columns (names ignored): the 1-based
  index of the state, the time of the event, and the value to add to
  that state at that time. This is the format of `Events` in
  [`cvsolve`](http://sn248.github.io/sundialr/reference/cvsolve.md), and
  the value is always added, including at the initial time, where it
  adds to `IC`. The solver integrates up to each event, applies it and
  restarts, carrying the sensitivities across unchanged: the event
  amounts and times are fixed, not sensitivity parameters, so an event
  moves the states but not their derivatives. Unlike `cvsolve`, the
  result still has one row per element of `time_vector`; an output time
  that coincides with an event reports the state just after the event.
  Default is NULL

## Value

A matrix with `1 + length(IC) * (length(Parameters) + length(sensIC))`
columns and one row per element of `time_vector`. The first column is
the time vector. The remaining columns hold the sensitivities in blocks
of `length(IC)`, one block per parameter: the sensitivities of all
states with respect to the first parameter, then of all states with
respect to the second parameter, and so on. The sensitivity of state `j`
with respect to parameter `k` is in column
`1 + (k - 1) * length(IC) + j`. The blocks for `sensIC` follow the
parameter blocks, in the order given in `sensIC`, so they are numbered
as if they were parameters `length(Parameters) + 1`,
`length(Parameters) + 2`, and so on.

The state solution is attached as the attribute `"solution"`: a matrix
with the same layout as the output of
[`cvode`](http://sn248.github.io/sundialr/reference/cvode.md), time in
the first column followed by one column per state. Retrieve it with
`attr(result, "solution")`. As with any attribute, subsetting the
returned matrix drops it, so extract it first.

## Examples

``` r
# Example of solving a set sensitivity equations for ODEs with cvodes function
# ODEs described by an R function
ODE_R <- function(t, y, p){

  # vector containing the right hand side gradients
  ydot = vector(mode = "numeric", length = length(y))

  # R indices start from 1
  ydot[1] = -p[1]*y[1] + p[2]*y[2]*y[3]
  ydot[2] = p[1]*y[1] - p[2]*y[2]*y[3] - p[3]*y[2]*y[2]
  ydot[3] = p[3]*y[2]*y[2]

  # ydot[1] = -0.04 * y[1] + 10000 * y[2] * y[3]
  # ydot[3] = 30000000 * y[2] * y[2]
  # ydot[2] = -ydot[1] - ydot[3]

  ydot

}

# R code to genrate time vector, IC and solve the equations
time_vec <- c(0.0, 0.4, 4.0, 40.0, 4E2, 4E3, 4E4, 4E5, 4E6, 4E7, 4E8, 4E9, 4E10)
IC <- c(1,0,0)
params <- c(0.04, 10000, 30000000)
reltol <- 1e-04
abstol <- c(1e-8,1e-14,1e-6)

## Solving the ODEs and Sensitivities using cvodes function
df1 <- cvodes(time_vec, IC, ODE_R , params, reltol, abstol,"STG",FALSE)           ## using R

## Solving with a manual Jacobian  J[i,j] = d(ydot_i)/d(y_j)
JAC_R <- function(t, y, p) {
  matrix(c(
    -p[1],         p[1],                        0,
     p[2]*y[3],   -p[2]*y[3] - 2*p[3]*y[2],   2*p[3]*y[2],
     p[2]*y[2],   -p[2]*y[2],                  0
  ), nrow = 3, ncol = 3)
}
df3 <- cvodes(time_vec, IC, ODE_R, params, reltol, abstol, "STG", FALSE, jacobian = JAC_R)

## Solving with a manual sensitivity right-hand side. CVODES calls it once per
## parameter with the 1-based index iS, and it returns
##   d(yS_iS)/dt = J %*% yS_iS + df/dp_iS  (a vector of length(y))
## This avoids the finite-difference approximation used when sensitivity = NULL.
SENS_R <- function(t, y, ydot, iS, yS, p) {
  J <- matrix(c(
    -p[1],         p[1],                        0,
     p[2]*y[3],   -p[2]*y[3] - 2*p[3]*y[2],   2*p[3]*y[2],
     p[2]*y[2],   -p[2]*y[2],                  0
  ), nrow = 3, ncol = 3)
  dfdp <- switch(iS,
                 c(-y[1],      y[1],       0),   # d f / d p1
                 c( y[2]*y[3], -y[2]*y[3], 0),   # d f / d p2
                 c( 0,        -y[2]^2,     y[2]^2))  # d f / d p3
  as.numeric(J %*% yS + dfdp)
}
df4 <- cvodes(time_vec, IC, ODE_R, params, reltol, abstol, "STG", FALSE, sensitivity = SENS_R)

## ODEs can also be described using Rcpp. Compiling the C++ version invokes the
## toolchain, which takes far longer than the solve itself, so it is not run here.
if (FALSE) { # \dontrun{
Rcpp::sourceCpp(code = '

                #include <Rcpp.h>
                using namespace Rcpp;

                // ODE functions defined using Rcpp
                // [[Rcpp::export]]
                NumericVector ODE_Rcpp (double t, NumericVector y, NumericVector p){

                // Initialize ydot filled with zeros
                NumericVector ydot(y.length());

                ydot[0] = -p[0]*y[0] + p[1]*y[1]*y[2];
                ydot[1] = p[0]*y[0] - p[1]*y[1]*y[2] - p[2]*y[1]*y[1];
                ydot[2] = p[2]*y[1]*y[1];

                return ydot;

                }')

df2 <- cvodes(time_vec, IC, ODE_Rcpp , params, reltol, abstol,"STG",FALSE)        ## using Rcpp

## Check that both solutions are identical
# identical(df1, df2)
} # }
```
