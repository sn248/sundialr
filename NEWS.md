sundialr v0.1.8
===============
* Updated the underlying SUNDIALS library to v7.8.0 (Jun 2026)
* CRAN-compliance patches moved to a dedicated script (`src/scripts/cran_patches.sh`) which now verifies that no flagged call survives patching
* Added `tools/strip_sundials_tarball.sh` to reproducibly strip the bundled SUNDIALS tarball of examples, docs, tests and benchmarks
* **Behaviour change**: an `Events` record at time 0 in `cvsolve()` now adds its value to the initial condition instead of replacing it, so the `IC` supplied by the caller is always the starting state and t = 0 behaves like every other event time. Several events at t = 0 on one state now accumulate rather than the last one winning. Results change for any call that doses at time 0: with `IC = c(1)` and a value of 100 at t = 0 the state now starts at 101 rather than 100, and `inst/examples/cvsolve_1D.r` starts at 101 accordingly. To reproduce the old numbers, subtract the initial condition from the value of any event at time 0
* Fixed `cvsolve()` treating the initial time as 0 rather than as the first element of `time_vector`. An event at the initial time is now applied to the initial conditions as intended; previously, if the time vector began anywhere other than 0, such an event was silently dropped and its value lost with no error or warning
* `cvsolve()` now rejects an event time outside the range of `time_vector` instead of accepting it. An event before the first output time cannot be applied, since the solve starts there; it previously produced a row of zeros in the returned matrix and discarded the event. Event times that are `NA` are rejected for the same reason
* `cvsolve()` no longer returns rows of zeros. A record that does not advance the solution — one at the initial time, or one sharing a time with the record before it — was skipped without storing anything, leaving that row filled with zeros instead of the state
* Fixed a bug in `ida()` and `cvsolve()` where a scalar `abstolerance` (the default) allocated a length-1 tolerance vector instead of one entry per state. On systems with more than one state this read and wrote past the end of that vector, and the solvers failed immediately with a corrector convergence error at the initial time. Passing `abstolerance` as a vector of the same length as `IC` was unaffected, as were `cvode()` and `cvodes()`
* `cvode()`, `cvodes()`, `ida()` and `cvsolve()` no longer leak memory when a call ends in an error. The SUNDIALS objects behind each solve are now released by a scope guard instead of by free calls at the end of the function, so they are also freed when an error unwinds the stack. Previously a failed call leaked everything allocated up to the point of the error, which accumulated when solvers were called repeatedly (for example from an optimiser or a fitting loop). This includes errors raised by SUNDIALS itself, such as a failure to converge: these previously aborted the call in a way that skipped all cleanup, and were the most common source of the leak, since a solver that fails to converge is exactly the case a user is likely to retry in a loop. Separately, `cvsolve()` never freed its constraints vector on any path, including successful ones
* `cvsolve()` now rejects an invalid state index in the first column of `Events`. The index is converted to a position in the state vector and used to subscript it directly, so a value outside `1:length(IC)` wrote outside that vector: an index of `0`, the natural mistake when counting from zero, corrupted the heap and could abort the `R` session. Fractional and `NA` indices are rejected for the same reason, a fractional index having previously been truncated to a different state without warning
* An error raised inside a user-supplied `R` function is now returned to the solver as a failure code rather than being thrown out through SUNDIALS' own C code, which lets SUNDIALS unwind its internal state instead of being abandoned part-way through a step. The `R` error itself is unchanged: the original condition still reaches the caller, so `tryCatch()` on a specific condition class continues to work
* A user-supplied `R` function returning the wrong number of values is now reported instead of being read past the end. Previously a right-hand side or residual function returning fewer values than there are states, or a Jacobian function returning a matrix of the wrong shape, was copied out element by element regardless and the solve continued on uninitialised memory
* Added test coverage for `ida()` and `cvsolve()`, which previously had none, including a check that a scalar and an equivalent vector `abstolerance` give identical results, and for the callback error handling above
* Error messages now name the SUNDIALS routine that actually failed. Where SUNDIALS reports a reason of its own, that message is passed through in full, giving the routine, source location and cause instead of a generic description. Five call sites also reported the wrong name, which is now corrected: `CVodeSetLinearSolver` was reported as `CVDlsSetLinearSolver` (a name removed in SUNDIALS 4) in `cvode()` and `cvsolve()`, `CVodeWFtolerances` as `CVodeSetEwtFn` and `CVodeSensInit1` as `CVodeSensInit` in `cvodes()`, and `IDASetUserData` and `IDASolve` as their `cvode()` equivalents in `ida()`

sundialr v0.1.7
===============
* Changes in `/tools/cmake_call.sh` file to remove warnings with respect to stdout/stderr/fprintf etc in upstream C library
* Updated the underlying SUNDIALS library to v7.6.0 (Jan 2026)

sundialr v0.1.6
================
Updated build system to use cmake

sundialr v0.1.5
================
Bug fixes:
* Fixed the bug in assigning absolute tolerances

sundialr v0.1.4.1
================
Bug fixes:
* Fixed the linking problem due to multiple defined symbols

sundialr v0.1.4
================

New features:
* added the `cvsolve` function to solve ODEs with multiple discontinuities

Bug fixes:
* None

API change:
* NONE is existing functions
* A new function `cvsolve` is added

Internal changes:
* Now uses SUNDIALS version 5.2.0 (released March 2020) at the back end
