# Example of solving a set of Differential Algebraic Equations (DAEs)
# with IDA function
# DAEs (residuals) described by an R function
DAE_R <- function(t, y, ydot, p){

  # vector containing the residuals
  res = vector(mode = "numeric", length = length(y))

  # R indices start from 1
  res[1] <- -0.04 * y[1] + 10000 * y[2] * y[3] - ydot[1]
  res[2] <- -res[1] - 30000000 * y[2] * y[2] - ydot[2]
  res[3] <- y[1] + y[2] + y[3] - 1.0

  res
}

# DAEs can also be described using Rcpp
Rcpp::sourceCpp(code = '

                #include <Rcpp.h>
                using namespace Rcpp;

                // ODE functions defined using Rcpp
                // [[Rcpp::export]]
                NumericVector DAE_Rcpp (double t, NumericVector y, NumericVector ydot, NumericVector p){

                // Initialize ydot filled with zeros
                NumericVector res(y.length());

                res[0] = -0.04 * y[0] + 10000 * y[1] * y[2];
                res[1] = -res[0] - 30000000 * y[1] * y[1] - ydot[1];
                res[0] = res[0] - ydot[0];
                res[2] = y[0] + y[1] + y[2] - 1.0;

                return res;

                }')

# R code to genrate time vector, IC and solve the equations
time_vec <- c(0.0, 0.4, 4.0, 40.0, 4E2, 4E3, 4E4, 4E5, 4E6, 4E7, 4E8, 4E9, 4E10)
IC <- c(1,0,0)
IRes <- c(-0.4, 0.4, 0)
params <- c(0.04, 10000, 30000000)
reltol <- 1e-04
abstol <- c(1e-8,1e-14,1e-6)

## Solving the DAEs using the ida function
df1 <- ida(time_vec, IC, IRes, DAE_R , params, reltol, abstol)           ## using R
df2 <- ida(time_vec, IC, IRes, DAE_Rcpp , params, reltol, abstol)        ## using Rcpp
