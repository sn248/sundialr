## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
DAE_R <- function(t, y, ydot, p){

  # vector containing the residuals
  res = vector(mode = "numeric", length = length(y))

  # R indices start from 1
  res[1] <- -0.04 * y[1] + 10000 * y[2] * y[3] - ydot[1]
  res[2] <- -res[1] - 30000000 * y[2] * y[2] - ydot[2]
  res[3] <- y[1] + y[2] + y[3] - 1.0

  res
}

# R code to genrate time vector, IC and solve the equations
time_vec <- c(0.0, 0.4, 4.0, 40.0, 4E2, 4E3, 4E4, 4E5, 4E6, 4E7, 4E8, 4E9, 4E10)
IC <- c(1,0,0)                  # Initial value of y
IRes <- c(-0.4, 0.4, 0)         # Initial value of ydot
params <- c(0.04, 10000, 30000000)
reltol <- 1e-04
abstol <- c(1e-8,1e-14,1e-6)

## Solving the DAEs using the ida function
df1 <- sundialr::ida(time_vec, IC, IRes, DAE_R , params, reltol, abstol) 

## -----------------------------------------------------------------------------
ODErepeated_R <- function(t, y, p){

  # vector containing the right hand side gradients
  ydot = vector(mode = "numeric", length = length(y))

  # R indices start from 1
  ydot[1] = -p[1]*y[1]

  ydot

}

## -----------------------------------------------------------------------------
TDOSE <- data.frame(ID = 1, TIMES = c(0, 10, 20, 30, 40, 50), VAL = 100)
TDOSE

## ---- fig.height=5, fig.width=7, fig.align="C"--------------------------------
# Example of solving a set of ODEs with multiple discontinuities using cvsolve
# A simple One dimensional equation, y = -0.1 * y
# ODEs described by an R function
ODE_R <- function(t, y, p){

  # vector containing the right hand side gradients
  ydot = vector(mode = "numeric", length = length(y))

  # R indices start from 1
  ydot[1] = -p[1]*y[1]

  ydot

}

# R code to generate time vector, IC and solve the equations
TSAMP <- seq(from = 0, to = 100, by = 0.1)      # sampling time points
IC <- c(1)
params <- c(0.1)

# A dataset describing the dosing at times at which additions to y[1] are to be done
# Names of the columns don't matter, but they MUST be in the order of state index,
# times and Values at discontinuity.
TDOSE <- data.frame(ID = 1, TIMES = c(0, 10, 20, 30, 40, 50), VAL = 100)

df1 <- sundialr::cvsolve(TSAMP, c(1), ODE_R, params)         # solving without any discontinuity
df2 <- sundialr::cvsolve(TSAMP, c(1), ODE_R, params, TDOSE)  # solving with discontinuity

## Plot the solution with discontinuities
## first column is time, second column is the state
time <- df2[,1]
y1 <- df2[,2]
plot(time, y1, type = "l", lty = 1, main = "An ODE system with discontinuties", frame.plot = F)   

