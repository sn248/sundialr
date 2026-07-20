# Changelog

## sundialr v0.1.8

- Updated the underlying SUNDIALS library to v7.8.0 (Jun 2026)
- CRAN-compliance patches moved to a dedicated script
  (`src/scripts/cran_patches.sh`) which now verifies that no flagged
  call survives patching
- Added `tools/strip_sundials_tarball.sh` to reproducibly strip the
  bundled SUNDIALS tarball of examples, docs, tests and benchmarks
- Fixed a bug in
  [`ida()`](http://sn248.github.io/sundialr/reference/ida.md) and
  [`cvsolve()`](http://sn248.github.io/sundialr/reference/cvsolve.md)
  where a scalar `abstolerance` (the default) allocated a length-1
  tolerance vector instead of one entry per state. On systems with more
  than one state this read and wrote past the end of that vector, and
  the solvers failed immediately with a corrector convergence error at
  the initial time. Passing `abstolerance` as a vector of the same
  length as `IC` was unaffected, as were
  [`cvode()`](http://sn248.github.io/sundialr/reference/cvode.md) and
  [`cvodes()`](http://sn248.github.io/sundialr/reference/cvodes.md)
- [`cvode()`](http://sn248.github.io/sundialr/reference/cvode.md),
  [`cvodes()`](http://sn248.github.io/sundialr/reference/cvodes.md),
  [`ida()`](http://sn248.github.io/sundialr/reference/ida.md) and
  [`cvsolve()`](http://sn248.github.io/sundialr/reference/cvsolve.md) no
  longer leak memory when a call ends in an error. The SUNDIALS objects
  behind each solve are now released by a scope guard instead of by free
  calls at the end of the function, so they are also freed when an error
  unwinds the stack. Previously a failed call leaked everything
  allocated up to the point of the error, which accumulated when solvers
  were called repeatedly (for example from an optimiser or a fitting
  loop). This includes errors raised by SUNDIALS itself, such as a
  failure to converge: these previously aborted the call in a way that
  skipped all cleanup, and were the most common source of the leak,
  since a solver that fails to converge is exactly the case a user is
  likely to retry in a loop. Separately,
  [`cvsolve()`](http://sn248.github.io/sundialr/reference/cvsolve.md)
  never freed its constraints vector on any path, including successful
  ones
- Error messages now name the SUNDIALS routine that actually failed.
  Where SUNDIALS reports a reason of its own, that message is passed
  through in full, giving the routine, source location and cause instead
  of a generic description. Five call sites also reported the wrong
  name, which is now corrected: `CVodeSetLinearSolver` was reported as
  `CVDlsSetLinearSolver` (a name removed in SUNDIALS 4) in
  [`cvode()`](http://sn248.github.io/sundialr/reference/cvode.md) and
  [`cvsolve()`](http://sn248.github.io/sundialr/reference/cvsolve.md),
  `CVodeWFtolerances` as `CVodeSetEwtFn` and `CVodeSensInit1` as
  `CVodeSensInit` in
  [`cvodes()`](http://sn248.github.io/sundialr/reference/cvodes.md), and
  `IDASetUserData` and `IDASolve` as their
  [`cvode()`](http://sn248.github.io/sundialr/reference/cvode.md)
  equivalents in
  [`ida()`](http://sn248.github.io/sundialr/reference/ida.md)

## sundialr v0.1.7

CRAN release: 2026-04-21

- Changes in `/tools/cmake_call.sh` file to remove warnings with respect
  to stdout/stderr/fprintf etc in upstream C library
- Updated the underlying SUNDIALS library to v7.6.0 (Jan 2026)

## sundialr v0.1.6

CRAN release: 2024-12-15

Updated build system to use cmake

## sundialr v0.1.5

CRAN release: 2024-11-11

Bug fixes: \* Fixed the bug in assigning absolute tolerances

## sundialr v0.1.4.1

CRAN release: 2021-05-16

Bug fixes: \* Fixed the linking problem due to multiple defined symbols

## sundialr v0.1.4

CRAN release: 2020-05-31

New features: \* added the `cvsolve` function to solve ODEs with
multiple discontinuities

Bug fixes: \* None

API change: \* NONE is existing functions \* A new function `cvsolve` is
added

Internal changes: \* Now uses SUNDIALS version 5.2.0 (released March
2020) at the back end
