context("Names on IC and Parameters reach the user's functions")

# Every solver copies the state into a fresh vector before calling R, which
# used to drop names(IC). A right-hand side written by name must now work
# without an adapter, in every callback.

test_that("cvode passes names to the RHS and the Jacobian", {
  seen <- new.env()
  rhs <- function(t, y, p) { seen$rhs <- names(y); seen$p <- names(p); -p[["k"]] * y[["A"]] }
  jac <- function(t, y, p) { seen$jac <- names(y); matrix(-p[["k"]]) }
  out <- cvode(c(0, 1), c(A = 2), rhs, c(k = 0.5), 1e-10, 1e-12, jacobian = jac)
  expect_equal(seen$rhs, "A")
  expect_equal(seen$jac, "A")
  expect_equal(seen$p, "k")
  expect_equal(out[2, 2], 2 * exp(-0.5), tolerance = 1e-7)
})

test_that("cvsolve passes names to the RHS", {
  seen <- new.env()
  rhs <- function(t, y, p) { seen$rhs <- names(y); -p[["k"]] * y[["A"]] }
  out <- cvsolve(c(0, 1), c(A = 2), rhs, c(k = 0.5), reltolerance = 1e-10,
                 abstolerance = 1e-12)
  expect_equal(seen$rhs, "A")
  expect_equal(out[2, 2], 2 * exp(-0.5), tolerance = 1e-7)
})

test_that("cvodes passes names to the RHS, Jacobian and sensitivity RHS", {
  seen <- new.env()
  rhs  <- function(t, y, p) { seen$rhs <- names(y); seen$p <- names(p); -p[["k"]] * y[["A"]] }
  jac  <- function(t, y, p) { seen$jac <- names(y); matrix(-p[["k"]]) }
  sens <- function(t, y, ydot, iS, yS, p) {
    seen$sens <- c(names(y), names(ydot), names(yS)); seen$sp <- names(p)
    -p[["k"]] * yS - if (iS == 1) y[["A"]] else 0
  }
  y0 <- c(A = 2); k <- c(k = 0.5); tt <- c(0, 1)
  # with and without sensIC: with it, the R function gets a separate copy of
  # the parameters, which must keep their names too
  for (s in list(NULL, "A")) {
    seen$sens <- seen$sp <- NULL
    out <- cvodes(tt, y0, rhs, k, 1e-10, 1e-12, jacobian = jac,
                  sensitivity = sens, sensIC = s)
    expect_equal(seen$rhs, "A")
    expect_equal(seen$jac, "A")
    expect_equal(seen$p, "k")
    expect_equal(seen$sens, c("A", "A", "A"))
    expect_equal(seen$sp, "k")
    expect_equal(out[2, 2], -1 * 2 * exp(-0.5), tolerance = 1e-6)
  }
})

test_that("ida passes names to the residual and the Jacobian", {
  seen <- new.env()
  res <- function(t, y, yp, p) { seen$res <- c(names(y), names(yp)); yp[["A"]] + p[["k"]] * y[["A"]] }
  jac <- function(t, y, yp, cj, p) { seen$jac <- c(names(y), names(yp)); matrix(p[["k"]] + cj) }
  out <- ida(c(0, 1), c(A = 2), c(A = -1), res, c(k = 0.5), 1e-10, 1e-12,
             jacobian = jac)
  expect_equal(seen$res, c("A", "A"))
  expect_equal(seen$jac, c("A", "A"))
  expect_equal(out[2, 2], 2 * exp(-0.5), tolerance = 1e-6)
})

test_that("unnamed IC still gives unnamed vectors", {
  seen <- new.env()
  rhs <- function(t, y, p) { seen$n <- names(y); -p * y }
  cvode(c(0, 1), 2, rhs, 0.5)
  expect_null(seen$n)
})
