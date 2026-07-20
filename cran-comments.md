## Comments for version 0.1.8
+ Updated the upstream `SUNDIALS` to version 7.8.0
+ CRAN-compliance patches (removal of `abort()`, `stdout`/`stderr` writes and `sprintf` from the bundled library) moved to a dedicated, self-verifying script
+ Fixed an out-of-bounds read/write in `ida()` and `cvsolve()`: a scalar `abstolerance` (the default) allocated a length-1 tolerance vector but one entry per state was written to and read from it. With more than one state both solvers failed at the initial time. The vector `abstolerance` form was unaffected, as were `cvode()` and `cvodes()`
+ `cvode()`, `cvodes()`, `ida()` and `cvsolve()` now release their `SUNDIALS` objects (context, vectors, matrix, linear/nonlinear solver, integrator memory) from a scope guard rather than from trailing free calls, so they are freed when an error unwinds the stack as well as on normal return. Previously any error raised between allocation and those free calls leaked every object allocated up to that point; these are C allocations, so `R`'s garbage collector did not reclaim them. The `SUNDIALS` error handler no longer raises the error itself either: it recorded the message and returns, and the solver reports it with `stop()`, because raising from the handler unwound with `longjmp` and so skipped the same cleanup. In `cvsolve()` the constraints vector was additionally never freed on any path
+ `src/Makevars.in` no longer sets `CXX`. It hardcoded `clang++`, which is not the package's choice to make and would fail on a machine without that compiler installed; the compiler is now left to `R`'s configuration

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
