context("sundialr C API (CVODE) - internal test drivers")

# The drivers are unexported (dot-prefixed), so reach them through the namespace.
decay      <- sundialr:::.capi_test_decay
twocmt     <- sundialr:::.capi_test_twocmt
forceerr   <- sundialr:::.capi_test_force_error
numsteps   <- sundialr:::.capi_test_num_steps
cleanerr   <- sundialr:::.capi_test_clean_err
abi        <- sundialr:::.capi_test_abi

test_that("scalar-tolerance decay matches the closed form", {
  # y' = -k y  ->  y(t) = y0 exp(-k t)
  times <- seq(0, 8, by = 0.5)
  y0 <- 3; k <- 0.6
  num <- decay(times, y0, k, tol_mode = 1L, reinit_each = FALSE)
  exact <- y0 * exp(-k * times)
  expect_equal(num, exact, tolerance = 1e-6)
})

test_that("reinit-per-segment gives the same trajectory as a single reinit", {
  # State continuity across CVodeReInit: reinitialising at every output point
  # with the current state must reproduce the uninterrupted solve (and the
  # closed form). This is the pattern a host loop drives per dose/event.
  times <- seq(0, 8, by = 0.5)
  y0 <- 3; k <- 0.6
  march  <- decay(times, y0, k, tol_mode = 1L, reinit_each = FALSE)
  perseg <- decay(times, y0, k, tol_mode = 1L, reinit_each = TRUE)
  exact  <- y0 * exp(-k * times)
  expect_equal(perseg, exact, tolerance = 1e-6)
  expect_equal(perseg, march, tolerance = 1e-6)
})

test_that("uniform vector tolerances equal scalar tolerances", {
  # itol = 4 path (CVodeWFtolerances) with uniform rtol/atol must agree with the
  # scalar (itol = 1) path - the same invariant the ida/cvsolve tests check.
  times <- seq(0, 8, by = 0.5)
  y0 <- 3; k <- 0.6
  sc <- decay(times, y0, k, tol_mode = 1L, reinit_each = FALSE)
  ve <- decay(times, y0, k, tol_mode = 2L, reinit_each = FALSE)
  expect_equal(ve, sc, tolerance = 1e-8)
})

test_that("two-compartment solve with analytic Jacobian matches the closed form", {
  # y1 = A e^{-k1 t};  y2 = A k1/(k2-k1) (e^{-k1 t} - e^{-k2 t})
  times <- seq(0, 10, by = 0.5)
  A <- 10; k1 <- 1.0; k2 <- 0.3
  num <- twocmt(times, A, k1, k2, use_jac = TRUE)
  y1 <- A * exp(-k1 * times)
  y2 <- A * k1 / (k2 - k1) * (exp(-k1 * times) - exp(-k2 * times))
  expect_equal(num[, 1], y1, tolerance = 1e-6)
  expect_equal(num[, 2], y2, tolerance = 1e-6)
})

test_that("analytic Jacobian agrees with the finite-difference default", {
  times <- seq(0, 10, by = 0.5)
  A <- 10; k1 <- 1.0; k2 <- 0.3
  with_jac <- twocmt(times, A, k1, k2, use_jac = TRUE)
  fd       <- twocmt(times, A, k1, k2, use_jac = FALSE)
  expect_equal(with_jac, fd, tolerance = 1e-6)
})

test_that("a failed solve returns a negative CVODE code and records a message", {
  res <- forceerr()
  # maxsteps = 1 over a long interval -> CV_TOO_MUCH_WORK (-1)
  expect_true(res$flag < 0)
  expect_equal(res$flag, -1)
  expect_true(nchar(res$err) > 0)
})

test_that("last_err is NULL after a successful solve", {
  expect_true(cleanerr())
})

test_that("introspection: ABI version and step count", {
  expect_equal(abi(), 1L)
  expect_true(numsteps(8.0, 3.0, 0.6) > 0)
})

test_that("set_udata lets one handle serve several parameter sets", {
  # Subject 1 solves with k1; set_udata repoints the callback data at k2 and
  # subject 2 solves on the SAME handle. Each column must match its own closed
  # form - if the old pointer were still read, column 2 would reproduce column 1.
  setud <- sundialr:::.capi_test_set_udata
  times <- seq(0, 8, by = 0.5)
  y0 <- 3; k1 <- 0.6; k2 <- 1.7
  num <- setud(times, y0, k1, k2)
  expect_equal(num[, 1], y0 * exp(-k1 * times), tolerance = 1e-6)
  expect_equal(num[, 2], y0 * exp(-k2 * times), tolerance = 1e-6)
})

test_that("a tight reinit+solve loop reuses the handle's memory and stays exact", {
  # The handle-reuse guarantee at host-loop volume: 1e4 reinit+solve segments
  # (a fit's order of magnitude) on one handle. ptr_stable pins that the state
  # buffer is never reallocated; y_end pins correctness after 1e4 restarts.
  tight <- sundialr:::.capi_test_tight_loop
  y0 <- 3; k <- 0.6
  res <- tight(10000L, 0.001, y0, k)
  expect_true(res$ptr_stable)
  expect_equal(res$y_end, y0 * exp(-k * res$t_end), tolerance = 1e-6)
  # get_num_steps accumulates across reinits (CVodeReInit zeroes CVODE's own
  # counter, so without the handle's running sum steps_end would be ~0 here).
  expect_true(res$steps_mid > 100)
  expect_true(res$steps_end > res$steps_mid)
})

test_that("reset_stats zeroes the step count and counting resumes", {
  rstats <- sundialr:::.capi_test_reset_stats
  res <- rstats()
  expect_true(res$before > 0)             # steps taken before the reset
  expect_equal(res$at_reset, 0)           # reads zero immediately, mid-segment
  expect_true(res$after_solve > 0)        # resumes counting
  expect_true(res$after_reinit > res$after_solve)  # and accumulates across reinit
})

test_that("independent handles can be driven from concurrent threads", {
  # One handle per thread, each with its own SUNContext and no shared state -
  # the documented thread-safety contract. Both trajectories must match their
  # closed forms after running simultaneously.
  conc <- sundialr:::.capi_test_concurrent
  times <- seq(0, 8, by = 0.25)
  y0 <- 3; k1 <- 0.6; k2 <- 1.7
  res <- conc(times, y0, k1, k2)
  expect_equal(res$y1, y0 * exp(-k1 * times), tolerance = 1e-6)
  expect_equal(res$y2, y0 * exp(-k2 * times), tolerance = 1e-6)
})
