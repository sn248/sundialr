## Comments for version 0.1.6
+ Updated the upstream `SUNDIALS` to version 7.2.0
+ Complete overall of the build system to use `cmake` based installation

## Comments for 0.1.5
+ Updated the upstream `SUNDIALS` to version 7.1.1. 
+ Fixed the `pkgdown` website
+ There was a bug in assigning absolute tolerance in equations. Fixed now.

## Comments for version 0.1.4.1.2
* Fixed the issue with forwarded email address.

## Comments for version 0.1.4.1
* Fixed the problem with linking due to multiple definitions. The downstream 
dependency should be unaffected by this change.

## Comments for version 0.1.4

## Test environments
* local macOS install, R 4.0.0
* ubuntu 16.04.6 LTS (on travis-ci), R 4.0.0
* win-builder (R-release)
* `rhub::check_on_cran()` 

## R win-builder check results
There were no ERRORs or WARNINGs or NOTEs.

## Downstream dependencies
There are currently one downstream dependency for this package (CausalKinetiX) 
which should be unaffected.
