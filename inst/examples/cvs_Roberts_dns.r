# Example of solving a set sensitivity equations for ODEs with cvodes function
# ODEs described by an R function
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

# ODEs can also be described using Rcpp
Rcpp::sourceCpp(code = '

                #include <Rcpp.h>
                using namespace Rcpp;

                // ODE functions defined using Rcpp
                // [[Rcpp::export]]
                NumericVector ODE_Rcpp (double t, NumericVector y, NumericVector p){

                // Initialize ydot filled with zeros
                NumericVector ydot(y.length());

                ydot[0] = -p[0]*y[0] + p[1]*y[1]*y[2];
                ydot[1] = p[0]*y[0] - p[1]*y[1]*y[2] - p[2]*y[1]*y[1];
                ydot[2] = p[2]*y[1]*y[1];

                return ydot;

                }')



# R code to genrate time vector, IC and solve the equations
time_vec <- c(0.0, 0.4, 4.0, 40.0, 4E2, 4E3, 4E4, 4E5, 4E6, 4E7, 4E8, 4E9, 4E10)
IC <- c(1,0,0)
params <- c(0.04, 10000, 30000000)
reltol <- 1e-04
abstol <- c(1e-8,1e-14,1e-6)

## Solving the ODEs and Sensitivities using cvodes function
df1 <- cvodes(time_vec, IC, ODE_R , params, reltol, abstol,"STG",FALSE)           ## using R
df2 <- cvodes(time_vec, IC, ODE_Rcpp , params, reltol, abstol,"STG",FALSE)        ## using Rcpp

## Check that both solutions are identical
# identical(df1, df2)
