## Comments for version 0.1.3

## Test environments
* local macOS install, R 3.6.0
* ubuntu 14.04 (on travis-ci), R 3.6.0
* win-builder (R-release)

## R win-builder check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  
Possibly mis-spelled words in DESCRIPTION:
CVODE (8:167)
CVODES (8:219)
Livermore (8:108)

Explanation: This is my first submission to CRAN. CVODE, CVODES and names of solvers in SUNDIALS C library, developed at Lawrence Livermore National Laboratory.

## Downstream dependencies
There are currently one downstream dependency for this package (CausalKinetiX) which breaks
as there due to the way differential equations are now defined. I have informed
the maintainer of the upcoming changes about a month back, they have agreed to change their
code once the newer version of `sundialr` is released.



  
  

 
