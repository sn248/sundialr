# cvsolve

CVSOLVE solver to solve stiff ODEs with discontinuties

## Usage

``` r
cvsolve(
  time_vector,
  IC,
  input_function,
  Parameters,
  Events = NULL,
  reltolerance = 1e-04,
  abstolerance = 1e-04,
  jacobian = NULL
)
```

## Arguments

- time_vector:

  time vector

- IC:

  Initial Conditions

- input_function:

  Right Hand Side function of ODEs

- Parameters:

  Parameters input to ODEs

- Events:

  Discontinuities in the solution (a DataFrame, default value is NULL).
  Three columns, names ignored: the 1-based index of the state, the time
  of the discontinuity, and the value to add to that state at that time.
  The value is always added to the current value of the state, including
  at the initial time, so the initial conditions in `IC` are the
  starting point and an event at t = 0 adds to them.

- reltolerance:

  Relative Tolerance (a scalar, default value = 1e-04)

- abstolerance:

  Absolute Tolerance (a scalar or vector with length equal to ydot,
  default = 1e-04)

- jacobian:

  (Optional) Jacobian of the RHS with signature `function(t, y, p)`
  returning an n-by-n matrix where entry \[i,j\] is d(ydot_i)/d(y_j).
  Default is NULL and SUNDIALS uses internal finite-difference
  approximation.

## Value

A Matrix. First column is the time-vector, the other columns are values
of y in order they are provided.

## Examples

``` r
# Example of solving a set of ODEs with multiple discontinuities using cvsolve
# A simple One dimensional equation, y = -0.1 * y
# ODEs described by an R function
ODE_R <- function(t, y, p){

  # vector containing the right hand side gradients
  ydot = vector(mode = "numeric", length = length(y))

  # R indices start from 1
  ydot[1] = -p[1]*y[1]

  ydot

}

# R code to generate time vector, IC and solve the equations
TSAMP <- seq(from = 0, to = 100, by = 0.1)      # sampling time points
IC <- c(1)
params <- c(0.1)

# A dataset describing the dosing at times at which additions to y[1] are to be done
# Names of the columns don't matter, but they MUST be in the order of state index,
# times and Values at discontinuity.
TDOSE <- data.frame(ID = 1, TIMES = c(0, 10, 20, 30, 40, 50), VAL = 100)
df1 <- cvsolve(TSAMP, c(1), ODE_R, params)         # solving without any discontinuity
df2 <- cvsolve(TSAMP, c(1), ODE_R, params, TDOSE, 0.001, 0.001, NULL)  # solving with discontinuity

## Solving with a manual Jacobian  J[1,1] = d(ydot[1])/d(y[1]) = -p[1]
JAC_R <- function(t, y, p) matrix(-p[1], nrow = 1, ncol = 1)
df3 <- cvsolve(TSAMP, IC, ODE_R, params, jacobian = JAC_R)
df4 <- cvsolve(TSAMP, IC, ODE_R, params, TDOSE, jacobian = JAC_R)
```
