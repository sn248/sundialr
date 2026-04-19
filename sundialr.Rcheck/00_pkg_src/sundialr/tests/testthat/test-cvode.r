context("Checking cvode Solution")

test_that("Size of solution is as expected", {

  ODE_R <- function(t, y, p){

    # vector containing the right hand side gradients
    ydot = vector(mode = "numeric", length = length(y))

    # R indices start from 1
    ydot[1] = -p[1]*y[1] + p[2]*y[2]*y[3]
    ydot[2] = p[1]*y[1] - p[2]*y[2]*y[3] - p[3]*y[2]*y[2]
    ydot[3] = p[3]*y[2]*y[2]

    # ydot[1] = -0.04 * y[1] + 10000 * y[2] * y[3]
    # ydot[3] = 30000000 * y[2] * y[2]
    # ydot[2] = -ydot[1] - ydot[3]

    ydot

  }

  # R code to genrate time vector, IC and solve the equations
  time_vec <- c(0.0, 0.4, 4.0, 40.0, 4E2, 4E3, 4E4, 4E5, 4E6, 4E7, 4E8, 4E9, 4E10)
  IC <- c(1,0,0)
  params <- c(0.04, 10000, 30000000)
  reltol <- 1e-04
  abstol <- c(1e-8,1e-14,1e-6)

  ## Solving the ODEs using cvode function
  df1 <- cvode(time_vec, IC, ODE_R , params, reltol, abstol)

  ## Expect solution to have same rows as time vector
  expect_equal(length(time_vec), nrow(df1))

  ## Expect solution to have same columns as IC + 1 (first column is time)
  expect_equal(length(IC) + 1, ncol(df1))

  ## checking for accuracy of the actual solution
  ## values in sundials example at last time point are considered actual solution
  ## the slight difference I see is could be due to the fact that sundials
  ## provides a Jacobian manually, numbers are very small for the first two tests
  expect_lt(abs(6.934511e-08 - df1[nrow(df1),2]), 1e-6)
  expect_lt(abs(2.773804e-13 - df1[nrow(df1),3]), 1e-6)
  expect_lt(abs(9.999999e-01 - df1[nrow(df1),4]), 1e-6)

})
