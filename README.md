# sundialr <img src="man/figures/sundialr_hex_image.png" align="right" alt="" width="150" />

<!-- badges: start -->

[![Build status](https://ci.appveyor.com/api/projects/status/3mp1p26lpqp16t3d?svg=true)](https://ci.appveyor.com/project/sn248/sundialr)    [![Travis-CI Build Status](https://travis-ci.org/sn248/sundialr.svg?branch=master)](https://travis-ci.org/sn248/sundialr)   [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/sundialr)](https://cran.r-project.org/package=sundialr) 
[![](https://cranlogs.r-pkg.org/badges/grand-total/sundialr)]( https://CRAN.R-project.org/package=sundialr)
[![DOI](https://zenodo.org/badge/60889710.svg)](https://zenodo.org/badge/latestdoi/60889710)
[![R-CMD-check](https://github.com/sn248/sundialr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sn248/sundialr/actions/workflows/R-CMD-check.yaml)
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
### Release 0.1.4.1
+ Fixed the linking bug due to multiple defined symbols. No other change.

### Release 0.1.4
+ This version has version 5.2.0 of `SUNDIALS` (released in March 2020) at the back end.
+ A new function `cvsolve` is added. It allows solving ODEs with multiple discontinuities in the solution. See a complete use case in the vignette.
+ A `pkgdown` create site of the package is added.
+ A hex sticker for the package is released!

### Release 0.1.3 
+ This version has version 4.0.1 of `SUNDIALS` (released in Dec 2018) at the back end.
+ An interface to `CVODES` is added. It calculates forward sensitivities w.r.t all parameters of the ODE system.
+ Parameters can now be defined as an input parameter to the ODE function. This will allow performing parameter optimization via numerical optimizers.

## Citations help! 

I will be much-obliged if you cite this package if you use it in your work. You can cite this package like this " we used the _sundialr_ R package (Satyaprakash Nayak, 2021)". Here is the full bibliographic reference to include in your reference list (don't forget to update the 'last accessed' date):

> Satyaprakash Nayak. (2021). sundialr: An Interface to 'SUNDIALS' Ordinary Differential Equation (ODE) Solvers (Version v1.4.1). Zenodo. https://doi.org/10.5281/zenodo.5501631. Last accessed 12 September 2021

The full citation entry is provided below:

    citEntry(
      entry    = "Misc",
      title    = "The sundialr package: An Interface to 'SUNDIALS' Ordinary Differential Equation (ODE) Solvers",
      author   = personList(as.person("Satyaprakash Nayak")),
      year     = "2021",
      doi      = "10.5281/zenodo.5501631",
      url      = "https://github.com/sn248/sundialr",
      textVersion =
        paste("Satyaprakash Nayak(2019).",
              "The sundialr package: An Interface to 'SUNDIALS' Ordinary Differential Equation (ODE) Solvers",
              "https://doi.org/10.5281/zenodo.5501631"
        )
    )

A BibTeX entry for LaTeX users can be generated using `citation("sundialr")` or can be copied from below:

    @Misc{,
    title = {The sundialr package: An Interface to 'SUNDIALS' Ordinary Differential Equation (ODE) Solvers},
    author = {Satyaprakash Nayak},
    year = {2021},
    doi = {10.5281/zenodo.5501631},
    url = {https://github.com/sn248/sundialr},
    }
