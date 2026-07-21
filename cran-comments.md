## Comments for version 0.1.8
+ Updated the upstream `SUNDIALS` to version 7.8.0
+ CRAN-compliance patches (removal of `abort()`, `stdout`/`stderr` writes and `sprintf` from the bundled library) moved to a dedicated, self-verifying script
+ Fixed `cvsolve()` assuming the initial time is 0 rather than reading it from `time_vector`. With a time vector starting elsewhere, an event at the initial time was silently dropped and its value lost. Event times outside the range of `time_vector`, and `NA` times, are now rejected rather than producing a row of zeros in the result, and rows that do not advance the solution are filled with the current state instead of being left zeroed
+ Behaviour change in `cvsolve()`: an `Events` record at time 0 now adds to the initial condition rather than replacing it, matching what events at every later time already did and what the argument is documented to mean. Calls that dose at time 0 return different values, so the first row of the `cvsolve` example changes from 100 to 101. There are no reverse dependencies affected; the only downstream package, `rxode2`, does not call `cvsolve()`
+ Fixed an out-of-bounds read/write in `ida()` and `cvsolve()`: a scalar `abstolerance` (the default) allocated a length-1 tolerance vector but one entry per state was written to and read from it. With more than one state both solvers failed at the initial time. The vector `abstolerance` form was unaffected, as were `cvode()` and `cvodes()`
+ `cvode()`, `cvodes()`, `ida()` and `cvsolve()` now release their `SUNDIALS` objects (context, vectors, matrix, linear/nonlinear solver, integrator memory) from a scope guard rather than from trailing free calls, so they are freed when an error unwinds the stack as well as on normal return. Previously any error raised between allocation and those free calls leaked every object allocated up to that point; these are C allocations, so `R`'s garbage collector did not reclaim them. The `SUNDIALS` error handler no longer raises the error itself either: it recorded the message and returns, and the solver reports it with `stop()`, because raising from the handler unwound with `longjmp` and so skipped the same cleanup. In `cvsolve()` the constraints vector was additionally never freed on any path
+ Fixed an out-of-bounds write in `cvsolve()`: the state index in the first column of the `Events` dataframe was only checked against an upper bound, but it is converted to a position in the state vector and used to subscript it directly, so an index of `0` or below wrote outside that vector and could abort the `R` session with a heap corruption error. The index is now required to be a whole number within range, which also rejects fractional indices that were previously truncated to a different state
+ Exceptions raised in the `R` callbacks are no longer allowed to propagate through SUNDIALS' C frames. The callbacks now catch them, record the message and return SUNDIALS' unrecoverable-error code, so the failing routine reports the error through its return value and SUNDIALS unwinds its own state. `R` errors and interrupts are deliberately re-thrown, since only `Rcpp`'s handler at the `.Call` boundary can resume the `R` unwind, and swallowing one would discard the condition
+ Fixed an out-of-bounds read when a user-supplied `R` function returns the wrong number of values: a right-hand side or residual function returning fewer values than there are states, or a Jacobian returning a matrix of the wrong shape, was copied out element by element without a length check. The lengths are now verified and a wrong one is reported
+ Added `tests/testthat/test-ida.r`, `tests/testthat/test-cvsolve.r` and `tests/testthat/test-callbacks.r`; `ida()` and `cvsolve()` previously had no automated tests, which is why the absolute tolerance fault above went unnoticed
+ The generated `inst/include/sundials/sundials_config.h` is now reproducible. `cmake` stamps it with the compiler, its full flag list and a build timestamp, so the committed copy carried one machine's build paths and changed on every build; those strings are provenance metadata and are now blanked after the `cmake` install step. The regenerated header also drops the `SUNDIALS_ARKODE` and `SUNDIALS_KINSOL` defines, which it had continued to declare even though the build disables both modules and ships neither library
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
