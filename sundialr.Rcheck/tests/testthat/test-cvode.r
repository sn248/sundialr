context("Checking cvode Solution")

ODE_R <- function(t, y, p){
  ydot = vector(mode = "numeric", length = length(y))
  ydot[1] = -p[1]*y[1] + p[2]*y[2]*y[3]
  ydot[2] = p[1]*y[1] - p[2]*y[2]*y[3] - p[3]*y[2]*y[2]
  ydot[3] = p[3]*y[2]*y[2]
  ydot
}

JAC_R <- function(t, y, p){
  matrix(c(
    -p[1],                              p[1],                   0,
     p[2]*y[3],   -p[2]*y[3] - 2*p[3]*y[2],   2*p[3]*y[2],
     p[2]*y[2],               -p[2]*y[2],                   0
  ), nrow = 3, ncol = 3)
}

time_vec <- c(0.0, 0.4, 4.0, 40.0, 4E2, 4E3, 4E4, 4E5, 4E6, 4E7, 4E8, 4E9, 4E10)
IC       <- c(1, 0, 0)
params   <- c(0.04, 10000, 30000000)
reltol   <- 1e-04
abstol   <- c(1e-8, 1e-14, 1e-6)

test_that("Size of solution is as expected", {

  df1 <- cvode(time_vec, IC, ODE_R, params, reltol, abstol)

  expect_equal(length(time_vec), nrow(df1))
  expect_equal(length(IC) + 1, ncol(df1))

  expect_lt(abs(6.934511e-08 - df1[nrow(df1),2]), 1e-6)
  expect_lt(abs(2.773804e-13 - df1[nrow(df1),3]), 1e-6)
  expect_lt(abs(9.999999e-01 - df1[nrow(df1),4]), 1e-6)

})

test_that("Manual Jacobian gives same solution as finite-difference approximation", {

  df_no_jac <- cvode(time_vec, IC, ODE_R, params, reltol, abstol)
  df_jac    <- cvode(time_vec, IC, ODE_R, params, reltol, abstol, jacobian = JAC_R)

  ## dimensions unchanged
  expect_equal(dim(df_no_jac), dim(df_jac))

  ## solutions agree to within some tolerance
  expect_lt(max(abs(df_no_jac - df_jac)), 1e-4)

})

test_that("Manual Jacobian solution meets known reference values", {

  df_jac <- cvode(time_vec, IC, ODE_R, params, reltol, abstol, jacobian = JAC_R)

  expect_lt(abs(6.934511e-08 - df_jac[nrow(df_jac),2]), 1e-6)
  expect_lt(abs(2.773804e-13 - df_jac[nrow(df_jac),3]), 1e-6)
  expect_lt(abs(9.999999e-01 - df_jac[nrow(df_jac),4]), 1e-6)

})
