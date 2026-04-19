# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## What this package is

`sundialr` is an R package (v0.1.7) wrapping the SUNDIALS C ODE-solving
library (v7.6.0). It exposes three solvers to R via Rcpp:

- `cvode` — stiff ODE solver (CVODE, BDF method, dense linear solver)
- `cvodes` — ODE solver with forward sensitivity equations (CVODES)
- `ida` — differential-algebraic equation solver (IDA)
- `cvsolve` — convenience wrapper around `cvode` for ODEs with multiple
  discontinuities (e.g. pharmacokinetic dosing simulations); uses
  RcppArmadillo

The user supplies an ODE RHS as an R or Rcpp function; the package does
not require SUNDIALS to be installed on the host system — it bundles and
builds SUNDIALS from source at install time via cmake.

## Build commands

``` r
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

## Build system architecture

SUNDIALS is not a pre-installed dependency — it is compiled from source
on every `R CMD INSTALL`:

1.  `configure` (generated from `configure.ac`) runs
    `tools/cmake_call.sh`
2.  `tools/cmake_call.sh` orchestrates:
    - Extracts `src/sundials-mod-7.6.0.tar.gz` → `src/sundials-src/`
    - Applies source patches (perl/python3 in-place edits) to make
      SUNDIALS CRAN-compliant: removes `abort()`, `fprintf(stderr,...)`,
      `printf`, `stdout`/`stderr` references, and replaces `sprintf`
      with `strcpy` in `cvodes`/`idas`
    - Runs cmake to build static SUNDIALS libs → installs headers to
      `inst/include/`, libs to `inst/lib/`
    - Moves compiled SUNDIALS C sources to `src/sundials/`; removes
      `src/sundials-src/` and `src/sundials-build/`
3.  `src/Makevars.in` → `src/Makevars` (substituted by configure): links
    the R package against the freshly built static libs in `inst/lib/`

**The `inst/include/` headers are the installed SUNDIALS headers
(committed). The actual SUNDIALS source tarball lives in
`src/sundials-mod-7.6.0.tar.gz`.**

## Key CRAN compliance constraints

Several non-obvious patches are applied to the SUNDIALS source at build
time to satisfy CRAN’s “Writing R Extensions §1.6.4” rules:

- All `abort()` calls removed; SUNDIALS errors are redirected to
  `Rf_error()` via a registered `SUNContext_PushErrHandler` (see
  `sundials_err_handler.h` included in each solver `.cpp`)
- All `fprintf(stderr,...)`, `printf`, `stdout`/`stderr` references in
  SUNDIALS source are removed or set to NULL
- `sprintf(name, "LITERAL")` replaced with `strcpy(name, "LITERAL")` in
  cvodes/idas to eliminate the `___sprintf_chk` symbol on macOS
- `SUNDIALS_LOGGING_LEVEL=0` cmake flag disables all logging output
- No explicit `CXX_STD = CXX11` in `Makevars.in` (CRAN flags it as
  unnecessary for modern R)

Any upgrade of the SUNDIALS version must re-apply all these patches in
`tools/cmake_call.sh`.

## Source file layout

| File                   | Purpose                                                                                                                                         |
|------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------|
| `src/cvode.cpp`        | [`cvode()`](http://sn248.github.io/sundialr/reference/cvode.md) R-callable function                                                             |
| `src/cvodes.cpp`       | [`cvodes()`](http://sn248.github.io/sundialr/reference/cvodes.md) R-callable function                                                           |
| `src/ida.cpp`          | [`ida()`](http://sn248.github.io/sundialr/reference/ida.md) R-callable function                                                                 |
| `src/cvsolve.cpp`      | [`cvsolve()`](http://sn248.github.io/sundialr/reference/cvsolve.md) R-callable function; uses RcppArmadillo, includes `src/utils/sortTimes.cpp` |
| `src/rhs_func.cpp`     | Shared `rhs_function()` callback; converts R `NumericVector` ↔︎ SUNDIALS `N_Vector`                                                              |
| `src/check_retval.cpp` | `check_retval()` helper used by all solvers                                                                                                     |
| `src/scripts/`         | cmake config, R toolchain detection, SUNDIALS extraction, cleanup scripts                                                                       |
| `inst/examples/`       | Runnable R examples referenced in Roxygen `@example` tags                                                                                       |

## ODE function interface

The user-supplied RHS function must have the signature `f(t, y, p)`
where `t` is time (scalar), `y` is the state vector (`NumericVector`),
and `p` is the parameter vector (`NumericVector`). It must return `ydot`
as a `NumericVector` of the same length as `y`. This contract is
enforced in `rhs_func.cpp` via `TYPEOF(rhs_fun) == CLOSXP`.
