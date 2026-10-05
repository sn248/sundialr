context("sundialr C API (ARKODE, ESDIRK) - internal test drivers")

# The drivers shared with test-capi.r take a trailing `solver` argument:
# 0 = CVODE, 1 = ARKODE with its default ESDIRK. Every contract pinned there for
# CVODE is re-checked here for ARKODE.
ARK <- 1L

decay    <- sundialr:::.capi_test_decay
twocmt   <- sundialr:::.capi_test_twocmt
method   <- sundialr:::.capi_test_arkode_method
roberts  <- sundialr:::.capi_test_robertson

test_that("ESDIRK decay matches the closed form, with and without reinits", {
  times <- seq(0, 8, by = 0.5)
  y0 <- 3; k <- 0.6
  exact  <- y0 * exp(-k * times)
  march  <- decay(times, y0, k, tol_mode = 1L, reinit_each = FALSE, solver = ARK)
  perseg <- decay(times, y0, k, tol_mode = 1L, reinit_each = TRUE,  solver = ARK)
  expect_equal(march,  exact, tolerance = 1e-6)
  expect_equal(perseg, exact, tolerance = 1e-6)
})

test_that("ESDIRK uniform vector tolerances equal scalar tolerances", {
  times <- seq(0, 8, by = 0.5)
  sc <- decay(times, 3, 0.6, tol_mode = 1L, reinit_each = FALSE, solver = ARK)
  ve <- decay(times, 3, 0.6, tol_mode = 2L, reinit_each = FALSE, solver = ARK)
  expect_equal(ve, sc, tolerance = 1e-8)
})

test_that("ESDIRK two-compartment solve matches the closed form, Jacobian or not", {
  times <- seq(0, 10, by = 0.5)
  A <- 10; k1 <- 1.0; k2 <- 0.3
  y1 <- A * exp(-k1 * times)
  y2 <- A * k1 / (k2 - k1) * (exp(-k1 * times) - exp(-k2 * times))
  with_jac <- twocmt(times, A, k1, k2, use_jac = TRUE,  solver = ARK)
  fd       <- twocmt(times, A, k1, k2, use_jac = FALSE, solver = ARK)
  expect_equal(with_jac[, 1], y1, tolerance = 1e-6)
  expect_equal(with_jac[, 2], y2, tolerance = 1e-6)
  expect_equal(with_jac, fd, tolerance = 1e-6)
})

test_that("ESDIRK and BDF agree on the stiff Robertson problem", {
  # Rate constants spanning nine orders of magnitude; an explicit method would
  # need millions of steps here. Mass is conserved exactly by the equations.
  times <- c(0.4, 4, 40, 400, 4e3, 4e4)
  ark   <- roberts(times, use_jac = TRUE,  solver = ARK)
  ark_fd <- roberts(times, use_jac = FALSE, solver = ARK)
  bdf   <- roberts(times, use_jac = TRUE,  solver = 0L)
  expect_equal(rowSums(ark), rep(1, length(times)), tolerance = 1e-6)
  expect_equal(ark, bdf, tolerance = 1e-4)
  expect_equal(ark_fd, ark, tolerance = 1e-4)
  # Reference value at t = 0.4 (Hairer & Wanner; the CVODE example output).
  expect_equal(ark[1, 1], 0.9851712, tolerance = 1e-5)
})

test_that("every ESDIRK order and a named table reach the closed form", {
  times <- seq(0, 8, by = 0.5)
  y0 <- 3; k <- 0.6
  exact <- y0 * exp(-k * times)
  steps <- numeric(0)
  for (q in 2:5) {
    res <- method(times, y0, k, order = q, table = "")
    # The order-2 table's error estimate comes from a first-order embedding, so
    # its global error at rtol = 1e-8 accumulates to ~5e-5; orders 3-5 stay
    # below 2e-7.
    expect_equal(res$y, exact, tolerance = if (q == 2) 1e-4 else 1e-6)
    steps[q - 1] <- res$steps
  }
  # At a fixed tolerance a higher-order method needs far fewer steps.
  expect_true(steps[4] < steps[1])
  named <- method(times, y0, k, order = 0L, table = "ARKODE_ESDIRK547L2SA2_7_4_5")
  expect_equal(named$y, exact, tolerance = 1e-6)
  # ARKODE's default order-5 table is that same ESDIRK.
  expect_equal(named$steps, steps[4])
})

test_that("bad method choices return ILL_INPUT with a message and leave the handle usable", {
  res <- sundialr:::.capi_test_arkode_method_errors()
  expect_equal(res$bad_order, -22)
  expect_match(res$bad_order_msg, "order")
  expect_equal(res$bad_name, -22)
  expect_match(res$bad_name_msg, "ARKODE_NO_SUCH_TABLE", fixed = TRUE)
  expect_equal(res$reinit, 0)
  expect_equal(res$late, -22)          # the method is fixed once integrating
  expect_equal(res$solve, 0)
  expect_equal(res$y1, exp(-1), tolerance = 1e-5)
})

test_that("a failed ESDIRK solve returns ARKODE's negative code and a message", {
  res <- sundialr:::.capi_test_force_error(solver = ARK)
  expect_equal(res$flag, -1)           # ARK_TOO_MUCH_WORK
  expect_true(nchar(res$err) > 0)
  expect_true(sundialr:::.capi_test_clean_err(solver = ARK))
  expect_true(sundialr:::.capi_test_num_steps(8.0, 3.0, 0.6, solver = ARK) > 0)
})

test_that("ESDIRK set_udata lets one handle serve several parameter sets", {
  times <- seq(0, 8, by = 0.5)
  num <- sundialr:::.capi_test_set_udata(times, 3, 0.6, 1.7, solver = ARK)
  expect_equal(num[, 1], 3 * exp(-0.6 * times), tolerance = 1e-6)
  expect_equal(num[, 2], 3 * exp(-1.7 * times), tolerance = 1e-6)
})

test_that("an ESDIRK tight reinit+solve loop reuses the handle's memory and stays exact", {
  res <- sundialr:::.capi_test_tight_loop(10000L, 0.001, 3, 0.6, solver = ARK)
  expect_true(res$ptr_stable)
  expect_equal(res$y_end, 3 * exp(-0.6 * res$t_end), tolerance = 1e-6)
  expect_true(res$steps_mid > 100)
  expect_true(res$steps_end > res$steps_mid)
})

test_that("ESDIRK reset_stats zeroes the step count and counting resumes", {
  res <- sundialr:::.capi_test_reset_stats(solver = ARK)
  expect_true(res$before > 0)
  expect_equal(res$at_reset, 0)
  expect_true(res$after_solve > 0)
  expect_true(res$after_reinit > res$after_solve)
})

test_that("independent ESDIRK handles can be driven from concurrent threads", {
  times <- seq(0, 8, by = 0.25)
  res <- sundialr:::.capi_test_concurrent(times, 3, 0.6, 1.7, solver = ARK)
  expect_equal(res$y1, 3 * exp(-0.6 * times), tolerance = 1e-6)
  expect_equal(res$y2, 3 * exp(-1.7 * times), tolerance = 1e-6)
})
