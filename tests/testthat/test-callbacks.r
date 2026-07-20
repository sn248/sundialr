context("Checking callback error handling")

## SUNDIALS calls the R functions from its own C code. An error there must come
## back as an ordinary, catchable R error rather than escaping through those C
## frames or being copied out of a wrongly sized return value.

time_vec <- c(0, 1, 2)
IC       <- c(1, 0, 0)
IRes     <- c(-0.5, 0.5, 0)
params   <- c(0.5)

ODE_R <- function(t, y, p) c(-p[1] * y[1], p[1] * y[1], 0)
DAE_R <- function(t, y, ydot, p) {
  c(-p[1] * y[1] - ydot[1], p[1] * y[1] - ydot[2], y[3] - y[1] - y[2])
}

## each solver called with a right-hand side / residual of the caller's choosing
with_rhs <- list(
  cvode   = function(f) cvode(time_vec, IC, f, params),
  cvodes  = function(f) cvodes(time_vec, IC, f, params, 1e-4, 1e-4, "STG", FALSE),
  cvsolve = function(f) cvsolve(time_vec, IC, f, params),
  ida     = function(f) ida(time_vec, IC, IRes, f, params)
)

test_that("An R error raised inside a callback surfaces as an R error", {

  for (solver in names(with_rhs)) {
    ## the residual takes an extra argument, so give each solver a signature
    ## it can call, and have both raise the same condition
    f <- if (solver == "ida") {
      function(t, y, ydot, p) stop("failure inside the user function")
    } else {
      function(t, y, p) stop("failure inside the user function")
    }

    err <- tryCatch(with_rhs[[solver]](f), error = function(e) e)

    expect_s3_class(err, "error")
    ## the original R message must survive, not be replaced by a generic one
    expect_match(conditionMessage(err), "failure inside the user function",
                 info = solver)
  }

})

test_that("A callback returning the wrong number of values is rejected", {

  for (solver in names(with_rhs)) {
    short <- if (solver == "ida") {
      function(t, y, ydot, p) c(1, 2)
    } else {
      function(t, y, p) c(1, 2)
    }
    long <- if (solver == "ida") {
      function(t, y, ydot, p) c(1, 2, 3, 4)
    } else {
      function(t, y, p) c(1, 2, 3, 4)
    }

    ## too few would otherwise be copied out of the end of the returned vector
    expect_error(with_rhs[[solver]](short), "same length as the state vector",
                 info = solver)
    expect_error(with_rhs[[solver]](long), "same length as the state vector",
                 info = solver)
  }

})

test_that("A Jacobian of the wrong shape is rejected", {

  ## cvode and cvsolve share the f(t, y, p) Jacobian signature
  for (solver in c("cvode", "cvsolve")) {
    call_with_jac <- function(j) {
      if (solver == "cvode") cvode(time_vec, IC, ODE_R, params, jacobian = j)
      else                   cvsolve(time_vec, IC, ODE_R, params, jacobian = j)
    }
    expect_error(call_with_jac(function(t, y, p) matrix(0, 2, 2)),
                 "must return a 3-by-3 matrix", info = solver)
    expect_error(call_with_jac(function(t, y, p) matrix(0, 3, 2)),
                 "must return a 3-by-3 matrix", info = solver)
  }

  ## ida's Jacobian takes cj and ydot as well
  expect_error(
    ida(time_vec, IC, IRes, DAE_R, params,
        jacobian = function(t, y, ydot, cj, p) matrix(0, 2, 2)),
    "must return a 3-by-3 matrix"
  )

})

test_that("A solver still works normally after callbacks have failed", {

  ## the failures above must leave nothing behind that breaks a later solve
  df <- cvode(time_vec, IC, ODE_R, params, 1e-8, 1e-10)

  expect_equal(nrow(df), length(time_vec))
  expect_lt(abs(df[3, 2] - exp(-params[1] * 2)), 1e-6)

})
