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
df1 <- cvsolve(TSAMP, c(1), ODE_R, params)         # solving without any discontinuity
df2 <- cvsolve(TSAMP, c(1), ODE_R, params, TDOSE)  # solving with discontinuity
