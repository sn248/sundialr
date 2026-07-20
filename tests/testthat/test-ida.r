context("Checking ida Solution")

## A well-conditioned index-1 DAE with a closed-form solution, so the tests
## assert against exact values rather than pinned regression numbers:
##   y1' = -k1*y1              ->  y1 = exp(-k1*t)
##   y2' =  k1*y1 - k2*y2      ->  y2 = k1/(k2-k1) * (exp(-k1*t) - exp(-k2*t))
##   0   =  y3 - (y1 + y2)     ->  y3 = y1 + y2      (algebraic)
DAE_R <- function(t, y, ydot, p) {
  c(
    -p[1] * y[1] - ydot[1],
     p[1] * y[1] - p[2] * y[2] - ydot[2],
     y[3] - y[1] - y[2]
  )
}

## J[i,j] = dF_i/dy_j + cj * dF_i/dydot_j; matrix() fills column-major
JAC_IDA <- function(t, y, ydot, cj, p) {
  matrix(c(
    -p[1] - cj,   p[1],         -1,
     0,          -p[2] - cj,    -1,
     0,           0,             1
  ), nrow = 3, ncol = 3)
}

analytic <- function(t, p) {
  y1 <- exp(-p[1] * t)
  y2 <- p[1] / (p[2] - p[1]) * (exp(-p[1] * t) - exp(-p[2] * t))
  cbind(y1, y2, y1 + y2)
}

time_vec <- seq(0, 20, by = 2)
IC       <- c(1, 0, 1)
IRes     <- c(-0.5, 0.5, 0)   # consistent: ydot = (-k1*y1, k1*y1-k2*y2, ydot1+ydot2)
params   <- c(0.5, 0.2)
reltol   <- 1e-8
abstol   <- rep(1e-10, 3)

test_that("Size of solution is as expected", {

  df1 <- ida(time_vec, IC, IRes, DAE_R, params, reltol, abstol)

  expect_equal(length(time_vec), nrow(df1))
  expect_equal(length(IC) + 1, ncol(df1))
  expect_equal(time_vec, df1[, 1])

})

test_that("Solution matches the closed-form solution", {

  df1 <- ida(time_vec, IC, IRes, DAE_R, params, reltol, abstol)

  expect_lt(max(abs(df1[, 2:4] - analytic(time_vec, params))), 1e-6)

})

test_that("Algebraic constraint is satisfied at every output time", {

  df1 <- ida(time_vec, IC, IRes, DAE_R, params, reltol, abstol)

  ## y3 - (y1 + y2) = 0 is enforced by the DAE, so it holds to the
  ## nonlinear-solve tolerance - far tighter than the solution accuracy,
  ## but not exactly zero
  expect_lt(max(abs(df1[, 4] - df1[, 2] - df1[, 3])), 1e-8)

})

test_that("Manual Jacobian gives same solution as finite-difference approximation", {

  df_no_jac <- ida(time_vec, IC, IRes, DAE_R, params, reltol, abstol)
  df_jac    <- ida(time_vec, IC, IRes, DAE_R, params, reltol, abstol,
                   jacobian = JAC_IDA)

  expect_equal(dim(df_no_jac), dim(df_jac))
  expect_lt(max(abs(df_no_jac - df_jac)), 1e-6)

  ## and the manual-Jacobian solution is itself correct
  expect_lt(max(abs(df_jac[, 2:4] - analytic(time_vec, params))), 1e-6)

})

test_that("Scalar and vector absolute tolerances agree", {

  ## Regression test: a scalar abstolerance is expanded to one entry per
  ## state internally, so it must give exactly the same answer as passing
  ## that value repeated. Sizing that expansion by the length of the input
  ## rather than by the number of states previously overran the buffer and
  ## made every multi-state solve with a scalar tolerance fail at t0.
  df_scalar <- ida(time_vec, IC, IRes, DAE_R, params, reltol, 1e-10)
  df_vector <- ida(time_vec, IC, IRes, DAE_R, params, reltol, rep(1e-10, 3))

  expect_equal(df_scalar, df_vector)
  expect_lt(max(abs(df_scalar[, 2:4] - analytic(time_vec, params))), 1e-6)

})

test_that("Invalid input is rejected", {

  ## IC and IRes must be the same length
  expect_error(ida(time_vec, IC, c(-0.5, 0.5), DAE_R, params, reltol, abstol))

  ## abstolerance must be a scalar or one value per state
  expect_error(ida(time_vec, IC, IRes, DAE_R, params, reltol, rep(1e-10, 2)))

})
