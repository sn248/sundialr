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

# Two-state chain y1 -> y2 with two parameters, used by the tests below.
# y1 = exp(-a t), y2 = a/(b - a) (exp(-a t) - exp(-b t)), so every state and
# every sensitivity has a closed form.
chain_rhs <- function(t, y, p) c(-p[1] * y[1], p[1] * y[1] - p[2] * y[2])
chain_t   <- c(0, 0.5, 1, 2, 4)
chain_y0  <- c(1, 0)
chain_p   <- c(0.5, 0.2)

chain_exact <- function(t, a, b) {
  e_a <- exp(-a * t); e_b <- exp(-b * t)
  y2  <- a / (b - a) * (e_a - e_b)
  list(
    y1 = e_a, y2 = y2,
    dy1_da = -t * e_a,
    dy1_db = 0 * t,
    dy2_da = b / (b - a)^2 * (e_a - e_b) - a * t * e_a / (b - a),
    dy2_db = -a / (b - a)^2 * (e_a - e_b) + a * t * e_b / (b - a)
  )
}

test_that("cvodes runs with its default arguments", {
  # the default ErrCon used to be the character 'F', which Rcpp refused to
  # convert to a logical, so a call relying on the defaults always failed
  expect_true(isFALSE(formals(cvodes)$ErrCon))
  out <- cvodes(chain_t, chain_y0, chain_rhs, chain_p)
  expect_equal(dim(out), c(length(chain_t), 1 + 2 * 2))
})

test_that("sensitivity columns are grouped by parameter, as documented", {
  # rtol 1e-8: with ErrCon = TRUE, difference-quotient sensitivities are too
  # noisy to pass the error test at 1e-10, and CVODES runs out of steps
  out <- cvodes(chain_t, chain_y0, chain_rhs, chain_p, 1e-8, 1e-10, ErrCon = TRUE)
  ex  <- chain_exact(chain_t, chain_p[1], chain_p[2])

  # no sensitivity RHS is supplied, so CVODES uses difference quotients whose
  # increment scales with sqrt(rtol): accurate to ~1e-6 here, far tighter than
  # needed to tell the four columns apart
  expect_equal(out[, 1], chain_t)
  # column 1 + (k - 1) * length(IC) + j holds d(y_j)/d(p_k)
  expect_equal(out[, 2], ex$dy1_da, tolerance = 1e-4)
  expect_equal(out[, 3], ex$dy2_da, tolerance = 1e-4)
  expect_equal(out[, 4], ex$dy1_db, tolerance = 1e-4)
  expect_equal(out[, 5], ex$dy2_db, tolerance = 1e-4)
})

test_that("the state solution is returned as the 'solution' attribute", {
  out  <- cvodes(chain_t, chain_y0, chain_rhs, chain_p, 1e-10, 1e-12)
  soln <- attr(out, "solution")
  ex   <- chain_exact(chain_t, chain_p[1], chain_p[2])

  expect_true(is.matrix(soln))
  expect_equal(dim(soln), c(length(chain_t), 1 + length(chain_y0)))
  expect_equal(soln[, 1], chain_t)
  expect_equal(soln[1, -1], chain_y0)
  expect_equal(soln[, 2], ex$y1, tolerance = 1e-6)
  expect_equal(soln[, 3], ex$y2, tolerance = 1e-6)

  # the same states cvode() computes
  ref <- cvode(chain_t, chain_y0, chain_rhs, chain_p, 1e-10, 1e-12)
  expect_equal(unname(soln), unname(ref), tolerance = 1e-6)
})

# Sensitivities with respect to initial conditions (sensIC). For the chain
# y1 -> y2 the states are linear in y(0), so
#   d(y1)/d(y1_0) = exp(-a t)        d(y1)/d(y2_0) = 0
#   d(y2)/d(y1_0) = a/(b-a) (exp(-a t) - exp(-b t))
#   d(y2)/d(y2_0) = exp(-b t)
chain_ic_exact <- function(t, a, b) {
  e_a <- exp(-a * t); e_b <- exp(-b * t)
  cbind(e_a, a / (b - a) * (e_a - e_b), 0 * t, e_b)
}
chain_y0_ic <- c(1, 0.3)
# Analytic sensitivity RHS written for the parameters alone, in the usual
# J %*% yS + df/dp[, iS] form. It indexes past the end of dfdp for any iS
# beyond length(p), so the stopifnot pins that sundialr never passes it an
# IC sensitivity: those are computed by the solver itself.
chain_jac  <- function(t, y, p) matrix(c(-p[1], p[1], 0, -p[2]), 2)
chain_sens <- function(t, y, ydot, iS, yS, p) {
  stopifnot(iS >= 1, iS <= length(p))
  dfdp <- cbind(c(-y[1], y[1]), c(0, -y[2]))
  as.vector(chain_jac(t, y, p) %*% yS) + dfdp[, iS]
}

test_that("sensIC appends IC sensitivity blocks after the parameter blocks", {
  # with a Jacobian, the IC sensitivities are J %*% yS exactly
  out <- cvodes(chain_t, chain_y0_ic, chain_rhs, chain_p, 1e-10, 1e-12,
                ErrCon = TRUE, jacobian = chain_jac, sensitivity = chain_sens,
                sensIC = 1:2)
  expect_equal(ncol(out), 1 + 2 * (2 + 2))

  # the parameter blocks are unaffected by the extra sensitivities
  ref <- cvodes(chain_t, chain_y0_ic, chain_rhs, chain_p, 1e-10, 1e-12,
                ErrCon = TRUE, jacobian = chain_jac, sensitivity = chain_sens)
  expect_equal(out[, 1:5], ref[, 1:5], tolerance = 1e-8)

  # yS(0) is the unit vector of each listed state
  expect_equal(out[1, 6:9], c(1, 0, 0, 1))
  expect_equal(unname(out[, 6:9]),
               unname(chain_ic_exact(chain_t, chain_p[1], chain_p[2])),
               tolerance = 1e-7)
})

test_that("a parameters-only sensitivity function never sees an IC index", {
  seen_iS <- integer(0)
  sens <- function(t, y, ydot, iS, yS, p) {
    seen_iS <<- union(seen_iS, iS)
    chain_sens(t, y, ydot, iS, yS, p)
  }
  # without a Jacobian the IC sensitivities come from a centred difference of
  # the RHS along yS, so they are only as accurate as that difference
  out <- cvodes(chain_t, chain_y0_ic, chain_rhs, chain_p, 1e-8, 1e-10,
                ErrCon = TRUE, sensitivity = sens, sensIC = 1:2)
  expect_setequal(seen_iS, 1:2)
  expect_equal(unname(out[, 6:9]),
               unname(chain_ic_exact(chain_t, chain_p[1], chain_p[2])),
               tolerance = 1e-4)
  # the parameter blocks still come from the user's function, exactly
  ex <- chain_exact(chain_t, chain_p[1], chain_p[2])
  expect_equal(out[, 2], ex$dy1_da, tolerance = 1e-7)
})

test_that("sensIC works with the difference-quotient sensitivity RHS", {
  # no sensitivity RHS: CVODES perturbs a placeholder the RHS never reads, so
  # the difference quotient reduces to J %*% yS, the IC sensitivity equation
  out <- cvodes(chain_t, chain_y0_ic, chain_rhs, chain_p, 1e-8, 1e-10,
                ErrCon = TRUE, sensIC = c(2, 1))
  ex  <- chain_ic_exact(chain_t, chain_p[1], chain_p[2])
  # blocks follow the order given in sensIC: state 2 first, then state 1
  expect_equal(unname(out[, 6:7]), unname(ex[, 3:4]), tolerance = 1e-4)
  expect_equal(unname(out[, 8:9]), unname(ex[, 1:2]), tolerance = 1e-4)
  # d(y1)/d(y2_0) is exactly zero: nothing couples y2 back into y1
  expect_equal(out[, 6], rep(0, length(chain_t)))
})

test_that("sensIC accepts state names", {
  y0 <- c(A = 1, B = 0.3)
  by_name  <- cvodes(chain_t, y0, chain_rhs, chain_p, 1e-10, 1e-12, sensIC = "B")
  by_index <- cvodes(chain_t, y0, chain_rhs, chain_p, 1e-10, 1e-12, sensIC = 2)
  expect_identical(by_name, by_index)
  expect_equal(ncol(by_name), 1 + 2 * 3)
})

test_that("sensIC is validated", {
  y0 <- c(A = 1, B = 0.3)
  run <- function(s) cvodes(chain_t, y0, chain_rhs, chain_p, sensIC = s)
  expect_error(run(0),       "whole numbers between 1 and length\\(IC\\)")
  expect_error(run(3),       "whole numbers between 1 and length\\(IC\\)")
  expect_error(run(1.5),     "whole numbers between 1 and length\\(IC\\)")
  expect_error(run(NA_real_), "whole numbers between 1 and length\\(IC\\)")
  expect_error(run(c(1, 1)), "more than once")
  expect_error(run("Z"),     "not in names\\(IC\\): Z")
  expect_error(run(TRUE),    "state indices or a vector of state names")
  expect_error(cvodes(chain_t, unname(y0), chain_rhs, chain_p, sensIC = "A"),
               "IC has no names")
  # NULL and an empty vector both mean no IC sensitivities
  expect_identical(run(NULL), run(integer(0)))
})

test_that("sensIC works when there are no parameters", {
  # y' = -0.5 y  =>  d(y)/d(y0) = exp(-0.5 t)
  out <- cvodes(c(0, 1, 2), c(A = 2), function(t, y, p) -0.5 * y, numeric(0),
                1e-8, 1e-10, ErrCon = TRUE, sensIC = "A")
  expect_equal(ncol(out), 2)
  expect_equal(out[, 2], exp(-0.5 * c(0, 1, 2)), tolerance = 1e-5)
})

# Events (doses). The chain is linear and a dose into y1 is a fixed amount,
# so by superposition every state and every parameter sensitivity is the
# undosed closed form plus, for each dose D at td, D times the same closed
# form started at td. The IC sensitivities do not change at all.
dosed_exact <- function(t, a, b, doses) {
  ex <- chain_exact(t, a, b)
  for (d in seq_len(nrow(doses))) {
    td <- doses$time[d]; D <- doses$value[d]
    s  <- pmax(t - td, 0); on <- t >= td
    exd <- chain_exact(s, a, b)
    for (nm in names(ex)) ex[[nm]] <- ex[[nm]] + on * D * exd[[nm]]
  }
  ex
}
dose_t <- seq(0, 12, by = 1)
# one dose off the output grid, one on it, and one at the final time
doses  <- data.frame(state = 1, time = c(2.5, 6, 12), value = c(2, 1.5, 3))

test_that("Events add doses and carry the sensitivities across them", {
  out  <- cvodes(dose_t, chain_y0, chain_rhs, chain_p, 1e-10, 1e-12,
                 Events = doses)
  soln <- attr(out, "solution")
  ex   <- dosed_exact(dose_t, chain_p[1], chain_p[2], doses)

  # one row per output time, whatever the events
  expect_equal(dim(out), c(length(dose_t), 1 + 2 * 2))
  expect_equal(out[, 1], dose_t)
  expect_equal(soln[, 1], dose_t)

  # an output time at a dose reports the state just after it
  expect_equal(soln[, 2], ex$y1, tolerance = 1e-6)
  expect_equal(soln[, 3], ex$y2, tolerance = 1e-6)
  expect_equal(out[, 2], ex$dy1_da, tolerance = 1e-5)
  expect_equal(out[, 3], ex$dy2_da, tolerance = 1e-5)
  expect_equal(out[, 4], ex$dy1_db, tolerance = 1e-5)
  expect_equal(out[, 5], ex$dy2_db, tolerance = 1e-5)
})

test_that("Events agree between the analytic and difference-quotient paths", {
  ana <- cvodes(dose_t, chain_y0_ic, chain_rhs, chain_p, 1e-10, 1e-12,
                jacobian = chain_jac, sensitivity = chain_sens,
                sensIC = 1:2, Events = doses)
  dq  <- cvodes(dose_t, chain_y0_ic, chain_rhs, chain_p, 1e-10, 1e-12,
                sensIC = 1:2, Events = doses)
  expect_equal(ana, dq, tolerance = 1e-5)

  # the IC sensitivities pass through the doses unchanged
  expect_equal(unname(ana[, 6:9]),
               unname(chain_ic_exact(dose_t, chain_p[1], chain_p[2])),
               tolerance = 1e-6)
})

test_that("an event at the initial time adds to IC and leaves yS(0) alone", {
  ev  <- data.frame(state = 2, time = 0, value = 0.5)
  out <- cvodes(chain_t, chain_y0, chain_rhs, chain_p, 1e-10, 1e-12,
                sensIC = 2, Events = ev)
  ref <- cvodes(chain_t, chain_y0 + c(0, 0.5), chain_rhs, chain_p, 1e-10, 1e-12,
                sensIC = 2)
  expect_equal(attr(out, "solution")[1, -1], c(1, 0.5))
  expect_equal(out[1, -1], c(0, 0, 0, 0, 0, 1))
  expect_equal(out, ref, tolerance = 1e-8)
})

test_that("Events at one time accumulate, and match cvsolve's states", {
  ev   <- data.frame(state = 1, time = c(3, 3, 7), value = c(1, 2, 1))
  out  <- cvodes(dose_t, chain_y0, chain_rhs, chain_p, 1e-10, 1e-12, Events = ev)
  ref  <- cvsolve(dose_t, chain_y0, chain_rhs, chain_p, ev, 1e-10, 1e-12)
  # cvsolve gives each event record its own row, so the two events at t = 3
  # make two rows there; its last row at each time is the post-event state
  ref <- ref[!duplicated(ref[, 1], fromLast = TRUE), ]
  expect_equal(unname(attr(out, "solution")), unname(ref), tolerance = 1e-6)
  ex <- dosed_exact(dose_t, chain_p[1], chain_p[2],
                    data.frame(time = c(3, 7), value = c(3, 1)))
  expect_equal(out[, 2], ex$dy1_da, tolerance = 1e-5)
})

test_that("No Events, or an empty Events, gives the undosed result", {
  base <- cvodes(chain_t, chain_y0, chain_rhs, chain_p, 1e-10, 1e-12)
  expect_identical(
    cvodes(chain_t, chain_y0, chain_rhs, chain_p, 1e-10, 1e-12, Events = NULL),
    base)
  expect_identical(
    cvodes(chain_t, chain_y0, chain_rhs, chain_p, 1e-10, 1e-12,
           Events = data.frame(state = numeric(0), time = numeric(0),
                               value = numeric(0))),
    base)
})

test_that("Events are validated", {
  run <- function(ev) cvodes(chain_t, chain_y0, chain_rhs, chain_p, Events = ev)
  expect_error(run(data.frame(state = 1, time = 1)), "exactly three columns")
  expect_error(run(data.frame(state = 3, time = 1, value = 1)), "whole number between 1")
  expect_error(run(data.frame(state = 0, time = 1, value = 1)), "whole number between 1")
  expect_error(run(data.frame(state = 1.5, time = 1, value = 1)), "whole number between 1")
  expect_error(run(data.frame(state = 1, time = -1, value = 1)), "between the first and last")
  expect_error(run(data.frame(state = 1, time = 5, value = 1)), "between the first and last")
  expect_error(run(data.frame(state = 1, time = NA, value = 1)), "between the first and last")
  expect_error(run(data.frame(state = 1, time = 1, value = NA)), "finite number")
})

test_that("ErrCon = TRUE puts the sensitivities under error control", {
  # Before 0.2.2 ErrCon was never passed to CVODES, so TRUE and FALSE gave
  # identical results. With the sensitivities in the error test they are
  # markedly more accurate at the same tolerance (about 30x here).
  t  <- seq(0, 20, by = 1)
  ex <- chain_exact(t, chain_p[1], chain_p[2])
  ex <- cbind(ex$dy1_da, ex$dy2_da, ex$dy1_db, ex$dy2_db)
  err <- function(ec, ev = NULL) {
    out <- cvodes(t, chain_y0, chain_rhs, chain_p, 1e-5, 1e-8, ErrCon = ec,
                  Events = ev)
    max(abs(out[, -1] - ex))
  }
  expect_lt(err(TRUE), err(FALSE) / 5)

  # the setting survives the restarts at events (an event at the final time
  # leaves every earlier row unchanged, so the closed form still applies)
  ev <- data.frame(state = 1, time = 20, value = 0)
  expect_lt(err(TRUE, ev), err(FALSE, ev) / 5)
})
