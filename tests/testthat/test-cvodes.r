context("Checking cvodes Solution Size")

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

  ## Solving the ODEs and Sensitivities using cvodes function
  df1 <- cvodes(time_vec, IC, ODE_R , params, reltol, abstol,"STG",FALSE)

  ## Expect solution to have same rows as time vector
  expect_equal(length(time_vec), nrow(df1))

  ## Expect solution to have same columns as length(IC) * length(params) + 1 (first column is time)
  expect_equal(length(IC) * length(params) + 1, ncol(df1))

})

test_that("analytic sensitivity RHS matches the closed-form sensitivity", {

  # Scalar decay: y' = -k y, y(0) = y0  =>  y(t) = y0 exp(-k t)
  # Sensitivity  s = dy/dk = -t y0 exp(-k t) = -t y(t), with s(0) = 0.
  ODE  <- function(t, y, p) -p[1] * y
  # sens RHS = J %*% yS + df/dk = -k * yS - y   (J = -k, df/dk = -y)
  SENS <- function(t, y, ydot, iS, yS, p) -p[1] * yS - y

  time_vec <- seq(0, 5, by = 0.5)
  y0 <- 2
  k  <- 0.7

  out <- cvodes(time_vec, y0, ODE, k, 1e-10, 1e-12, "STG", FALSE,
                sensitivity = SENS)

  # columns: time, sensitivity of y1 w.r.t k
  s_num   <- out[, 2]
  s_exact <- -time_vec * y0 * exp(-k * time_vec)

  expect_equal(s_num, s_exact, tolerance = 1e-5)
})

test_that("analytic sensitivity RHS agrees with finite-difference sensitivities", {

  # Roberts problem; the analytic sensitivity RHS should reproduce the
  # finite-difference sensitivities CVODES computes when `sensitivity` is NULL.
  ODE_R <- function(t, y, p){
    ydot    = vector(mode = "numeric", length = length(y))
    ydot[1] = -p[1]*y[1] + p[2]*y[2]*y[3]
    ydot[2] =  p[1]*y[1] - p[2]*y[2]*y[3] - p[3]*y[2]*y[2]
    ydot[3] =  p[3]*y[2]*y[2]
    ydot
  }

  # J[i,j] = d(ydot_i)/d(y_j); dfdp[[iS]] = d(ydot)/d(p_iS)
  SENS_R <- function(t, y, ydot, iS, yS, p) {
    J <- matrix(c(
      -p[1],        p[1],                      0,
       p[2]*y[3],  -p[2]*y[3] - 2*p[3]*y[2],   2*p[3]*y[2],
       p[2]*y[2],  -p[2]*y[2],                 0
    ), nrow = 3, ncol = 3)
    dfdp <- switch(iS,
                   c(-y[1],      y[1],       0),
                   c( y[2]*y[3], -y[2]*y[3], 0),
                   c( 0,        -y[2]^2,     y[2]^2))
    as.numeric(J %*% yS + dfdp)
  }

  time_vec <- c(0.0, 0.4, 4.0, 40.0, 4E2, 4E3, 4E4)
  IC     <- c(1, 0, 0)
  params <- c(0.04, 10000, 30000000)
  reltol <- 1e-06
  abstol <- c(1e-10, 1e-14, 1e-8)

  fd  <- cvodes(time_vec, IC, ODE_R, params, reltol, abstol, "STG", FALSE)
  ana <- cvodes(time_vec, IC, ODE_R, params, reltol, abstol, "STG", FALSE,
                sensitivity = SENS_R)

  expect_equal(dim(fd), dim(ana))
  # relative agreement; sensitivities span many orders of magnitude
  expect_equal(ana, fd, tolerance = 1e-3)
})
