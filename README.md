[![Build status](https://ci.appveyor.com/api/projects/status/3mp1p26lpqp16t3d?svg=true)](https://ci.appveyor.com/project/sn248/sundialr)    [![Travis-CI Build Status](https://travis-ci.org/sn248/sundialr.svg?branch=master)](https://travis-ci.org/sn248/sundialr)   [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/sundialr)](https://cran.r-project.org/package=sundialr)

# sundialr

A wrapper around the SUNDIALS ODE solving C library.

# What's new?

## Release 0.1.3 
+ The branch works with the version 4.0.1 of SUNDIALS (released in Dec 2018)
+ An interface to `CVODES` is added. It calculates forward sensitivities w.r.t all parameters of the ODE system.
+ Parameters can now be defined as an input parameter to the ODE function. This will allow performing parameter optimization via numerical optimizers.
