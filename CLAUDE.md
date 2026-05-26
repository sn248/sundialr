# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package is

`sundialr` is an R package (v0.1.7) wrapping the SUNDIALS C ODE-solving library (v7.6.0). It exposes four solvers to R via Rcpp:

- `cvode` — stiff ODE solver (CVODE, BDF method, dense linear solver)
- `cvodes` — ODE solver with forward sensitivity equations (CVODES)
- `ida` — differential-algebraic equation (DAE) solver (IDA)
- `cvsolve` — convenience wrapper around `cvode` for ODEs with multiple discontinuities (e.g. pharmacokinetic dosing simulations); uses RcppArmadillo

The package bundles SUNDIALS and builds it from source at install time via cmake — no system SUNDIALS installation required.

## Build commands

```r
# Standard R package build/check
R CMD build .
R CMD check sundialr_*.tar.gz --as-cran

# Install from source (triggers the cmake build of SUNDIALS)
R CMD INSTALL .

# Run all tests
Rscript -e "devtools::test()"

# Run a single test file
Rscript -e "testthat::test_file('tests/testthat/test-cvode.R')"

# Regenerate Rcpp exports and documentation
Rscript -e "Rcpp::compileAttributes(); devtools::document()"
```

If `configure.ac` is edited, regenerate `configure` with `autoconf` before committing.

## Build system architecture

SUNDIALS is compiled from source on every `R CMD INSTALL`:

1. `configure` (generated from `configure.ac`) runs `tools/cmake_call.sh`
2. `tools/cmake_call.sh` orchestrates:
   - Extracts `src/sundials-mod-7.6.0.tar.gz` → `src/sundials-src/`
   - Applies source patches (perl/python3 in-place edits) to make SUNDIALS CRAN-compliant: removes `abort()`, `fprintf(stderr,...)`, `printf`, `stdout`/`stderr` references, and replaces `sprintf` with `strcpy` in `cvodes`/`idas`
   - Runs cmake to build static SUNDIALS libs → installs headers to `inst/include/`, libs to `inst/lib/`
   - Moves compiled SUNDIALS C sources to `src/sundials/`; removes `src/sundials-src/` and `src/sundials-build/`
3. `src/Makevars.in` → `src/Makevars` (substituted by configure): links the R package against the freshly built static libs in `inst/lib/`

**The `inst/include/` headers are the installed SUNDIALS headers (committed). The actual SUNDIALS source tarball lives in `src/sundials-mod-7.6.0.tar.gz`. The `src/sundials/` directory contains the compiled SUNDIALS C sources moved there post-build (also committed).**

## Key CRAN compliance constraints

Several non-obvious patches are applied to the SUNDIALS source at build time to satisfy CRAN's "Writing R Extensions §1.6.4" rules:

- All `abort()` calls removed; SUNDIALS errors are redirected to `Rf_error()` via a registered `SUNContext_PushErrHandler` (see `sundials_err_handler.h` in `inst/include/`, included in each solver `.cpp`)
- All `fprintf(stderr,...)`, `printf`, `stdout`/`stderr` references in SUNDIALS source are removed or set to NULL
- `sprintf(name, "LITERAL")` replaced with `strcpy(name, "LITERAL")` in cvodes/idas to eliminate the `___sprintf_chk` symbol on macOS
- `SUNDIALS_LOGGING_LEVEL=0` cmake flag disables all logging output
- No explicit `CXX_STD = CXX11` in `Makevars.in` (CRAN flags it as unnecessary for modern R)

Any upgrade of the SUNDIALS version must re-apply all these patches in `tools/cmake_call.sh`.

## Source file layout

| File | Purpose |
|------|---------|
| `src/cvode.cpp` | `cvode()` R-callable function |
| `src/cvodes.cpp` | `cvodes()` R-callable function |
| `src/ida.cpp` | `ida()` R-callable function |
| `src/cvsolve.cpp` | `cvsolve()` R-callable function; uses RcppArmadillo, includes `src/utils/sortTimes.cpp` |
| `src/rhs_func.cpp` | Shared `rhs_function()` callback; converts R `NumericVector` ↔ SUNDIALS `N_Vector` |
| `src/check_retval.cpp` | `check_retval()` helper used by all solvers |
| `src/scripts/` | cmake config, R toolchain detection, SUNDIALS extraction, cleanup scripts |
| `inst/include/` | Installed SUNDIALS headers plus three custom package headers (see below) |
| `inst/examples/` | Runnable R examples referenced in Roxygen `@example` tags |

Custom package headers in `inst/include/` (not part of SUNDIALS):
- `rhs_func.h` — declares `struct rhs_func` and `rhs_function()` used by cvode/cvsolve
- `check_retval.h` — declares `check_retval()` used by all solvers
- `sundials_err_handler.h` — the static `sundials_r_err_handler` that routes SUNDIALS fatal errors to `Rf_error()`; must be called via `SUNContext_PushErrHandler(sunctx, sundials_r_err_handler, NULL)` after `SUNContext_Create()`

## Solver function interfaces

### cvode / cvsolve
The user-supplied RHS function must have the signature `f(t, y, p)` where `t` is time (scalar), `y` is the state vector (`NumericVector`), and `p` is the parameter vector (`NumericVector`). It must return `ydot` as a `NumericVector` of the same length as `y`. This contract is enforced in `rhs_func.cpp` via `TYPEOF(rhs_fun) == CLOSXP`.

### cvodes
Same RHS signature as `cvode`. Additional parameters:
- `SensType`: `"STG"` (staggered, default) or `"SIM"` (simultaneous) — controls how sensitivity equations are solved relative to the state equations
- `ErrCon`: `TRUE`/`FALSE` — whether sensitivities are included in error control

The return matrix has `1 + length(IC) * length(params)` columns: time followed by sensitivity of each state with respect to each parameter.

### ida
The residual function has signature `f(t, y, ydot, p)` where `ydot` is the current derivative vector — the function returns residuals `F(t, y, y', p) = 0` as a `NumericVector`. Call signature: `ida(time_vector, IC, IRes, input_function, Parameters, reltolerance, abstolerance)` where `IRes` is the initial value of `ydot` (must be consistent with the DAE at `t0`).

### cvsolve Events
The `Events` argument is a `data.frame` with exactly three columns (names ignored): state index (1-based integer), event time, and value to add to that state at that time. Example: `data.frame(state=1, time=c(10,20), value=100)` adds 100 to `y[1]` at t=10 and t=20.

## Tests

Tests live in `tests/testthat/`. Only `cvode` and `cvodes` have test files (`test-cvode.r`, `test-cvodes.r`). There are no automated tests for `ida` or `cvsolve`; the `inst/examples/` files serve as manual verification.
