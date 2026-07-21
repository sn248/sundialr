context("Checking cvsolve Solution")

## First-order decay, y' = -k*y, so y = y0 * exp(-k*t) between discontinuities.
ODE_R <- function(t, y, p) {
  ydot <- vector(mode = "numeric", length = length(y))
  ydot[1] <- -p[1] * y[1]
  ydot
}

JAC_R <- function(t, y, p) matrix(-p[1], nrow = 1, ncol = 1)

## Two-compartment transit chain, used where more than one state is needed.
ODE_3 <- function(t, y, p) {
  c(-p[1] * y[1], p[1] * y[1] - p[2] * y[2], p[2] * y[2])
}

TSAMP  <- seq(0, 20, by = 2)
IC     <- c(1)
params <- c(0.1)
reltol <- 1e-8
abstol <- 1e-10

test_that("Size of solution is as expected without events", {

  df1 <- cvsolve(TSAMP, IC, ODE_R, params)

  expect_equal(length(TSAMP), nrow(df1))
  expect_equal(length(IC) + 1, ncol(df1))
  expect_equal(TSAMP, df1[, 1])

})

test_that("Solution without events matches the closed-form solution", {

  df1 <- cvsolve(TSAMP, IC, ODE_R, params, NULL, reltol, abstol)

  expect_lt(max(abs(df1[, 2] - IC[1] * exp(-params[1] * TSAMP))), 1e-6)

})

test_that("An event adds its value to the state at the event time", {

  ## Dose deliberately placed off the sampling grid so it contributes its
  ## own output row and the jump is unambiguous.
  t_dose <- 9
  dose   <- 10
  TDOSE  <- data.frame(state = 1, time = t_dose, value = dose)

  df1 <- cvsolve(TSAMP, IC, ODE_R, params, TDOSE, reltol, abstol)

  ## the event time appears as an extra output row
  expect_equal(nrow(df1), length(TSAMP) + 1)
  expect_true(t_dose %in% df1[, 1])

  ## closed form: decay to the dose, add the dose, then decay onward.
  ## The value reported AT the dose time already includes the dose.
  expected <- ifelse(
    df1[, 1] < t_dose,
    IC[1] * exp(-params[1] * df1[, 1]),
    (IC[1] * exp(-params[1] * t_dose) + dose) *
      exp(-params[1] * (df1[, 1] - t_dose))
  )

  expect_lt(max(abs(df1[, 2] - expected)), 1e-6)

  ## the jump at the dose is the dose amount
  i <- which(df1[, 1] == t_dose)
  expect_lt(abs((df1[i, 2] - IC[1] * exp(-params[1] * t_dose)) - dose), 1e-6)

})

test_that("An event at t = 0 adds to the initial condition", {

  ## t = 0 is not a special case: the event adds to the IC the user supplied,
  ## exactly as events at later times add to the current state. Before 0.1.8
  ## an event here replaced the IC, discarding it.
  dose  <- 10
  TDOSE <- data.frame(state = 1, time = 0, value = dose)
  df1   <- cvsolve(TSAMP, IC, ODE_R, params, TDOSE, reltol, abstol)

  expect_equal(df1[1, 2], IC[1] + dose)

  ## and the whole trajectory decays from that starting value
  expect_lt(max(abs(df1[, 2] - (IC[1] + dose) * exp(-params[1] * df1[, 1]))), 1e-6)

  ## several events at t = 0 on the same state accumulate
  df2 <- cvsolve(TSAMP, IC, ODE_R, params,
                 data.frame(state = 1, time = c(0, 0), value = c(10, 5)),
                 reltol, abstol)
  expect_equal(df2[1, 2], IC[1] + 15)

})

test_that("The initial time is taken from the time vector, not assumed to be 0", {

  ## The solve starts at time_vector[0]. An event there must be applied to the
  ## initial condition exactly as one at t = 0 is; before 0.1.8 the initial
  ## time was hard-coded as 0, so with a time vector starting elsewhere the
  ## event was silently dropped and the dose vanished without an error.
  TS    <- seq(5, 20, by = 5)
  dose  <- 100
  df1   <- cvsolve(TS, IC, ODE_R, params,
                   data.frame(state = 1, time = 5, value = dose), reltol, abstol)

  expect_equal(df1[1, 2], IC[1] + dose)
  ## relative, not absolute: the state starts at 101 here, so an absolute
  ## bound would really be testing the magnitude rather than the accuracy
  expected <- (IC[1] + dose) * exp(-params[1] * (TS - 5))
  expect_lt(max(abs(df1[, 2] - expected) / expected), 1e-6)

  ## an event later in the same window is unaffected
  df2 <- cvsolve(TS, IC, ODE_R, params,
                 data.frame(state = 1, time = 10, value = dose), reltol, abstol)
  expect_lt(abs(df2[2, 2] - (IC[1] * exp(-params[1] * 5) + dose)), 1e-6)

})

test_that("Events outside the output window are rejected", {

  TS <- seq(0, 20, by = 5)

  ## before the first output time: the solve starts there, so it cannot be
  ## applied. This used to emit a row of zeros and drop the event.
  expect_error(cvsolve(TS, IC, ODE_R, params,
                       data.frame(state = 1, time = -5, value = 100),
                       reltol, abstol),
               "between the first and last time points")

  ## after the last output time
  expect_error(cvsolve(TS, IC, ODE_R, params,
                       data.frame(state = 1, time = 99, value = 100),
                       reltol, abstol),
               "between the first and last time points")

  expect_error(cvsolve(TS, IC, ODE_R, params,
                       data.frame(state = 1, time = NA, value = 100),
                       reltol, abstol),
               "between the first and last time points")

  ## both ends of the window are valid
  expect_silent(cvsolve(TS, IC, ODE_R, params,
                        data.frame(state = 1, time = 0, value = 100), reltol, abstol))
  expect_silent(cvsolve(TS, IC, ODE_R, params,
                        data.frame(state = 1, time = 20, value = 100), reltol, abstol))

})

test_that("No output row is left unfilled when events share a time", {

  ## Records that do not advance the solution - one at the initial time, or one
  ## at the same time as the previous record - used to be skipped without
  ## storing anything, leaving the row zero-filled.
  df <- cvsolve(seq(0, 10, by = 2), IC, ODE_R, params,
                data.frame(state = 1, time = c(0, 0), value = c(10, 5)),
                reltol, abstol)

  expect_false(any(df[, 2] == 0))
  expect_equal(df[1, 2], IC[1] + 15)
  ## the duplicate row at the same time carries the same state
  expect_equal(df[2, 1], 0)
  expect_equal(df[2, 2], df[1, 2])

})

test_that("Repeated events accumulate", {

  TDOSE <- data.frame(state = 1, time = c(5, 15), value = 10)
  df1   <- cvsolve(TSAMP, IC, ODE_R, params, TDOSE, reltol, abstol)

  y_at <- function(tt) df1[which(df1[, 1] == tt), 2]

  ## after the first dose
  expect_lt(abs(y_at(5) - (IC[1] * exp(-params[1] * 5) + 10)), 1e-6)

  ## after the second, carrying the decayed remainder of the first
  before2 <- (IC[1] * exp(-params[1] * 5) + 10) * exp(-params[1] * (15 - 5))
  expect_lt(abs(y_at(15) - (before2 + 10)), 1e-6)

})

test_that("Manual Jacobian gives same solution as finite-difference approximation", {

  TDOSE <- data.frame(state = 1, time = 9, value = 10)

  df_no_jac <- cvsolve(TSAMP, IC, ODE_R, params, TDOSE, reltol, abstol)
  df_jac    <- cvsolve(TSAMP, IC, ODE_R, params, TDOSE, reltol, abstol,
                       jacobian = JAC_R)

  expect_equal(dim(df_no_jac), dim(df_jac))
  expect_lt(max(abs(df_no_jac - df_jac)), 1e-6)

})

test_that("Scalar and vector absolute tolerances agree", {

  ## Regression test: see the matching test in test-ida.r. Needs more than
  ## one state, which is exactly the case the old sizing bug broke.
  IC3 <- c(1, 0, 0)
  p3  <- c(0.5, 0.2)

  df_scalar <- cvsolve(TSAMP, IC3, ODE_3, p3, NULL, reltol, 1e-10)
  df_vector <- cvsolve(TSAMP, IC3, ODE_3, p3, NULL, reltol, rep(1e-10, 3))

  expect_equal(df_scalar, df_vector)

  ## y1 is an independent decay with a known solution, and the three
  ## states conserve mass
  expect_lt(max(abs(df_scalar[, 2] - exp(-p3[1] * TSAMP))), 1e-6)
  expect_lt(max(abs(rowSums(df_scalar[, 2:4]) - 1)), 1e-6)

})

test_that("Invalid input is rejected", {

  ## abstolerance must be a scalar or one value per state
  expect_error(cvsolve(TSAMP, c(1, 0, 0), ODE_3, c(0.5, 0.2), NULL,
                       reltol, rep(1e-10, 2)))

  ## an event cannot refer to a state that does not exist
  expect_error(cvsolve(TSAMP, IC, ODE_R, params,
                       data.frame(state = 5, time = 9, value = 10),
                       reltol, abstol))

  ## The state index is turned into a 0-based subscript into the state
  ## vector, so anything that is not a whole number in 1..n would index
  ## outside it. Each of these must be rejected rather than corrupting memory.
  IC3 <- c(1, 0, 0)
  p3  <- c(0.5, 0.2)
  bad_state <- function(s) {
    cvsolve(TSAMP, IC3, ODE_3, p3, data.frame(state = s, time = 9, value = 10),
            reltol, 1e-10)
  }
  expect_error(bad_state(0))    # 0-based typo -> subscript -1
  expect_error(bad_state(-1))   # negative    -> subscript -2
  expect_error(bad_state(1.5))  # fractional  -> silently truncated to state 1
  expect_error(bad_state(NA))   # NA          -> undefined subscript
  expect_error(bad_state(4))    # above n     -> subscript past the end

  ## the valid boundaries are still accepted
  expect_silent(bad_state(1))
  expect_silent(bad_state(3))

  ## an event cannot happen after the last output time
  expect_error(cvsolve(TSAMP, IC, ODE_R, params,
                       data.frame(state = 1, time = 999, value = 10),
                       reltol, abstol))

})
