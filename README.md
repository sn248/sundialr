# sundialr <img src="man/figures/sundialr_hex_image.png" align="right" alt="" width="150" />

<!-- badges: start -->

[![Build status](https://ci.appveyor.com/api/projects/status/3mp1p26lpqp16t3d?svg=true)](https://ci.appveyor.com/project/sn248/sundialr)    [![Travis-CI Build Status](https://travis-ci.org/sn248/sundialr.svg?branch=master)](https://travis-ci.org/sn248/sundialr)   [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/sundialr)](https://cran.r-project.org/package=sundialr) [![](https://cranlogs.r-pkg.org/badges/sundialr)]( https://CRAN.R-project.org/package=sundialr)

<!-- badges: end -->

## sundialr

`sundialr` is a wrapper around a few of the solvers in the `SUNDIALS` ODE solving C library produced by the Lawrence Livermore National Laboratory and Southern Methodist University. More
information about `SUNDIALS` can be found [here](https://computing.llnl.gov/projects/sundials).
`SUNDIALS` is one of the most popular and well-respected ODE solving libraries available and 
`sundialr` provides a way to interface some of the `SUNDIALS` solvers in `R`. 

Currently `sundialr` provides an interface to the `serial` versions of `cvode` (for solving ODES), `cvodes` (for solving ODE with sensitivity equations) and `ida` (for solving differential-algebraic equations) using the Linear Solver (dense version).

A convenience function `cvsolve` is provided which allows solving a system of equations with
multiple discontinutities in solution. An application of such a system of equations would be 
to simulate the effect of multiple bolus doses of a drug in clinical pharmacokinetics. See the 
vignette for more details.

## What's new?

### Release 0.1.4
+ This version has version 5.2.0 of `SUNDIALS` (released in March 2020) at the back end.
+ A new function `cvsolve` is added. It allows solving ODEs with multiple discontinuities in the solution. See a complete use case in the vignette.
+ A `pkgdown` create site of the package is added.
+ A hex sticker for the package is released!

### Release 0.1.3 
+ This version has version 4.0.1 of `SUNDIALS` (released in Dec 2018) at the back end.
+ An interface to `CVODES` is added. It calculates forward sensitivities w.r.t all parameters of the ODE system.
+ Parameters can now be defined as an input parameter to the ODE function. This will allow performing parameter optimization via numerical optimizers.
