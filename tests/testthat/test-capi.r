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
