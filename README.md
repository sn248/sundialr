# sundialr <img src="man/figures/sundialr_hex_image.png" align="right" alt="" width="150" />

<!-- badges: start -->

[![CRAN](https://www.r-pkg.org/badges/version/sundialr)](https://cran.r-project.org/package=sundialr) 
[![R-CMD-check](https://github.com/sn248/sundialr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sn248/sundialr/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

## sundialr

`sundialr` is a wrapper around a few of the solvers in the `SUNDIALS` ODE solving C library produced by the Lawrence Livermore National Laboratory and Southern Methodist University. More
information about `SUNDIALS` can be found [here](https://computing.llnl.gov/projects/sundials).
`SUNDIALS` is one of the most popular and well-respected ODE solving libraries available and 
`sundialr` provides a way to interface some of the `SUNDIALS` solvers in `R`. 

Currently `sundialr` provides an interface to the `serial` versions of `cvode` (for solving ODES), `cvodes` (for solving ODE with sensitivity equations) and `ida` (for solving differential-algebraic equations) using the Linear Solver (dense version).

A convenience function `cvsolve` is provided which allows solving a system of equations with
multiple discontinuities in solution. An application of such a system of equations would be 
to simulate the effect of multiple bolus doses of a drug in clinical pharmacokinetics. See the 
vignette for more details.

## Installation

```r
install.packages("sundialr")
```

The development version can be installed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("sn248/sundialr")
```

`SUNDIALS` itself is bundled with the package and is compiled from source during
installation, so no system installation of `SUNDIALS` is needed. Building from
source does require `cmake` to be available.

## Quick start

Every solver takes the system of equations as an `R` function. For `cvode`, that
function is given the time, the state vector and the parameter vector, and
returns the derivatives:

```r
library(sundialr)

# The Roberts problem, a classic stiff system of three equations
ODE_R <- function(t, y, p) {
  ydot <- vector(mode = "numeric", length = length(y))

  ydot[1] <- -p[1] * y[1] + p[2] * y[2] * y[3]
  ydot[2] <-  p[1] * y[1] - p[2] * y[2] * y[3] - p[3] * y[2] * y[2]
  ydot[3] <-  p[3] * y[2] * y[2]

  ydot
}

time_vec <- c(0.0, 0.4, 4.0, 40.0, 4e2, 4e3, 4e4)
IC       <- c(1, 0, 0)
params   <- c(0.04, 1e4, 3e7)
reltol   <- 1e-6
abstol   <- c(1e-8, 1e-14, 1e-6)

# The result is a matrix: the first column is time, the remaining columns are
# the states, in the order they were given in IC
soln <- cvode(time_vec, IC, ODE_R, params, reltol, abstol)
```

`reltolerance` is a single value. `abstolerance` is either a single value,
applied to every state, or one value per state as above, which is useful when
the states differ in magnitude.

### Supplying a Jacobian

By default `SUNDIALS` approximates the Jacobian by finite differences. For a
stiff system it is usually faster, and more accurate, to supply it. The function
takes the same arguments as the system itself and returns an n-by-n matrix whose
`[i, j]` entry is `d(ydot[i])/d(y[j])`:

```r
JAC_R <- function(t, y, p) {
  matrix(
    c(-p[1],       p[1],                                  0,
       p[2]*y[3], -p[2]*y[3] - 2*p[3]*y[2],   2*p[3]*y[2],
       p[2]*y[2], -p[2]*y[2],                             0),
    nrow = 3, ncol = 3
  )
}

soln <- cvode(time_vec, IC, ODE_R, params, reltol, abstol, jacobian = JAC_R)
```

The system and the Jacobian can also be written in `C++` with `Rcpp` instead of
in `R`, which avoids the cost of calling back into `R` at every step. The
examples in `?cvode` show both forms.

## The other solvers

**`cvodes`** solves the same systems as `cvode` and additionally computes
forward sensitivities, the partial derivative of every state with respect to
every parameter, alongside the solution. These are what parameter estimation and
identifiability work need. The sensitivity equations can be solved either
staggered or simultaneously with the state equations, and can be included in or
excluded from the error control. See `?cvodes`.

**`ida`** solves differential-algebraic equations, systems that mix differential
equations with algebraic constraints and so cannot be written as `y' = f(t, y)`.
The system is given as a residual `F(t, y, y', p) = 0`, and an initial value of
`y'` consistent with the constraints at the initial time is required alongside
the initial values of `y`. See `?ida`.

**`cvsolve`** is a convenience wrapper around `cvode` for systems whose solution
has discontinuities. Discontinuities are described by a data frame giving the
state, the time, and the value to add to that state at that time; the solver
integrates up to each one, applies it, and restarts. The motivating case is
simulating repeated bolus doses of a drug in clinical pharmacokinetics, where
each dose is an instantaneous addition to a compartment. See `?cvsolve` and the
vignette.

## What's new?

### Version 0.2.0
+ **New**: all four solvers accept an optional `jacobian` argument for supplying the Jacobian analytically instead of letting `SUNDIALS` approximate it by finite differences. It defaults to `NULL`, so existing code is unaffected. See the quick start above for `cvode`; `ida` takes a different form, documented in `?ida`
+ Updated the underlying `SUNDIALS` library to v7.8.0 (Jun 2026)
+ **Behaviour change**: in `cvsolve`, an event at time 0 now adds its value to the initial conditions rather than replacing them, so the initial conditions are always the starting state and time 0 behaves like every other event time. Results change for any call that applies an event at time 0
+ `cvsolve` now takes the initial time from the first element of the time vector rather than assuming it is 0, and validates the state indices and times given in `Events`
+ Fixed several memory faults, including an out-of-bounds access when `abstolerance` was given as a single value to `ida` or `cvsolve`, and memory that was not released when a solve ended in an error
+ Errors raised inside a user-supplied `R` function are now reported properly, and error messages name the `SUNDIALS` routine that failed
+ Added test coverage for `ida` and `cvsolve`

For the full list, and for earlier versions, see [NEWS.md](NEWS.md).
