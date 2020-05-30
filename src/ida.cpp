//   Copyright (c) 2020, Satyaprakash Nayak
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions are
//   met:
//
//   Redistributions of source code must retain the above copyright
//   notice, this list of conditions and the following disclaimer.
//
//   Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in
//   the documentation and/or other materials provided with the
//   distribution.
//
//   Neither the name of the <ORGANIZATION> nor the names of its
//   contributors may be used to endorse or promote products derived
//   from this software without specific prior written permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <Rcpp.h>

#include <ida/ida.h>                          /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>           /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h>        /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h>        /* access to dense SUNLinearSolver      */
#include <sunnonlinsol/sunnonlinsol_newton.h> /* access to Newton SUNNonlinearSolver  */
#include <sundials/sundials_types.h>          /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>           /* defs. of SUNRabs, SUNRexp, etc.      */

#include <check_retval.h>

using namespace Rcpp;

// struct to use if R or Rcpp function is input as RHS residual function
struct res_func{
  Function res_eqn;
  NumericVector params;
};


// function called by IDAInit if user inputs R function
int res_function(realtype t, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data){

  // convert y to NumericVector y1
  int yy_len = NV_LENGTH_S(yy);

  // NumericVector analog of yy
  NumericVector yy1(yy_len);                       // filled with zeros
  realtype *yy_ptr = N_VGetArrayPointer(yy);
  for (int i = 0; i < yy_len; i++){
    yy1[i] = yy_ptr[i];                           // Ith(y,i+1);
  }

  int yp_len = NV_LENGTH_S(yp);
  // NumericVector analog of yp
  NumericVector yp1(yp_len);                      // filled with zeros
  realtype *yp_ptr = N_VGetArrayPointer(yp);
  for (int i = 0; i < yp_len; i++){
    yp1[i] = yp_ptr[i];                          // Ith(y,i+1);
  }

  // NumericVector analog of rr
  int rr_len = NV_LENGTH_S(rr);
  NumericVector rr1(rr_len);                     // filled with zeros

  if(!user_data){
    stop("Something went wrong, stopping!");
  }

  // cast void pointer to pointer to struct and assign rhs residual to a Function
  struct res_func *my_res_fun = (struct res_func*)user_data;

  if(my_res_fun){

    Function res_fun = (*my_res_fun).res_eqn;           // function
    NumericVector p_values = (*my_res_fun).params;      // rate parameters

    if (res_fun && (TYPEOF(res_fun) == CLOSXP)){
      // use the function to calculate value of RHS ----
      rr1 = res_fun(t, yy1, yp1, p_values);
    }
    else{
      stop("Something went wrong, stopping!");
    }
  }
  else {
    stop("Something went wrong, stopping!");
  }

  // convert NumericVector rr1 to N_Vector rr
  realtype *rr_ptr = N_VGetArrayPointer(rr);
  for (int i = 0; i <  rr_len; i++){
    rr_ptr[i] = rr1[i];
  }

  // everything went smoothly
  return(0);
}
//---RHS residual function definition ends ----------------------------------------------

//------------------------------------------------------------------------------
//'ida
//'
//' IDA solver to solve stiff DAEs
//'@param time_vector time vector
//'@param IC Initial Value of y
//'@param IRes Inital Value of ydot
//'@param input_function Right Hand Side function of DAEs
//'@param Parameters Parameters input to ODEs
//'@param reltolerance Relative Tolerance (a scalar, default value  = 1e-04)
//'@param abstolerance Absolute Tolerance (a scalar or vector with length equal to ydot, default = 1e-04)
// [[Rcpp::export]]
NumericMatrix ida(NumericVector time_vector, NumericVector IC, NumericVector IRes, SEXP input_function,
                  NumericVector Parameters,
                  double reltolerance = 0.0001, NumericVector abstolerance = 0.0001){


  int time_vec_len = time_vector.length();
  int y_len = IC.length();

  if(y_len != IRes.length()){
    stop("IC and IRes should be of same length");
  }

  int abstol_len = abstolerance.length();

  int flag;
  realtype reltol = reltolerance;

  realtype T0 = RCONST(time_vector[0]);     //RCONST(0.0);  // Initial Time

  double time;
  int NOUT = time_vec_len;

  // Set the vector absolute tolerance -----------------------------------------
  // abstol must be same length as IC
  N_Vector abstol = N_VNew_Serial(abstol_len);
  realtype *abstol_ptr = N_VGetArrayPointer(abstol);
  if(abstol_len == 1){
    // if a scalar is provided - use it to make a vector with same values
    for (int i = 0; i<y_len; i++){
      abstol_ptr[i] = abstolerance[0];
    }
  }
  else if (abstol_len == y_len){

    for (int i = 0; i<abstol_len; i++){
      abstol_ptr[i] = abstolerance[i];
    }
  }
  else if(abstol_len != 1 || abstol_len != y_len){

    stop("Absolute tolerance must be a scalar or a vector of same length as IC \n");
  }
  //----------------------------------------------------------------------------
  // // Set the initial values of y -----------------------------------------------
  N_Vector yy0 = N_VNew_Serial(y_len);          // declared as yy0 to be consistent with example C code
  realtype *yy0_ptr = N_VGetArrayPointer(yy0);
  for (int i = 0; i<y_len; i++){
    yy0_ptr[i] = IC[i];
  }

  // // Set the initial values of ydot --------------------------------------------
  N_Vector yp0 = N_VNew_Serial(y_len);         // declared as yp0 to be consistent with example C code
  realtype *yp0_ptr = N_VGetArrayPointer(yp0);
  for (int i = 0; i<y_len; i++){
    yp0_ptr[i] = IRes[i];
  }

  // //----------------------------------------------------------------------------
  void *ida_mem;
  ida_mem = NULL;

  /* Call IDACreate and IDAInit to initialize IDA memory */
  ida_mem = IDACreate();

  if(check_retval((void *)ida_mem, "IDACreate", 0)) {
    stop("Stopping IDA, something went wrong in allocating memory!");
    }

  // -- assign user input to the struct based on SEXP type of input_function
  if(!input_function){
    stop("Something is wrong with the input function, stopping!");
  }

  switch(TYPEOF(input_function)){

  case CLOSXP:{

    struct res_func my_res_function = {input_function, Parameters};

    // setting the user data in the rhs residual function
    flag = IDASetUserData(ida_mem, (void*)&my_res_function);
    if (check_retval(&flag, "CVodeSetUserData", 1)) { stop("Stopping IDA, something went wrong in setting user data!"); }

    flag = IDAInit(ida_mem, res_function, T0, yy0, yp0);
    if(check_retval(&flag, "IDAInit", 1)) { stop("Stopping, something went wrong in initializing IDA!"); };

    /* Call IDASVtolerances to set tolerances */
    flag = IDASVtolerances(ida_mem, reltol, abstol);
    if(check_retval(&flag, "IDASVtolerances", 1)) { stop("Stopping, something went wrong in setting tolerances!"); };

    /* Call IDARootInit to specify the root function grob with 2 components */
    // retval = IDARootInit(mem, 2, grob);
    // if (check_retval(&retval, "IDARootInit", 1)) return(1);

    /* Create dense SUNMatrix for use in linear solves */
    sunindextype y_len_M = y_len;
    SUNMatrix SM = SUNDenseMatrix(y_len_M, y_len_M);
    if(check_retval((void *)SM, "SUNDenseMatrix", 0)) { stop("Stopping IDA, something went wrong in setting the dense matrix!"); }


    // Create dense SUNLinearSolver object for use by IDA
    SUNLinearSolver LS = SUNLinSol_Dense(yy0, SM);
    if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) { stop("Stopping IDA, something went wrong in setting the linear solver!"); }

    /* Attach the matrix and linear solver */
    flag = IDASetLinearSolver(ida_mem, LS, SM);
    if(check_retval(&flag, "IDASetLinearSolver", 1)) return(1);

    // /* Set the user-supplied Jacobian routine */
    // retval = IDASetJacFn(mem, jacrob);
    // if(check_retval(&retval, "IDASetJacFn", 1)) return(1);

    /* Create Newton SUNNonlinearSolver object. IDA uses a
     * Newton SUNNonlinearSolver by default, so it is unecessary
     * to create it and attach it. */
    SUNNonlinearSolver NLS = SUNNonlinSol_Newton(yy0);
    if(check_retval((void *)NLS, "SUNNonlinSol_Newton", 0)) return(1);

    /* Attach the nonlinear solver */
    flag = IDASetNonlinearSolver(ida_mem, NLS);
    if(check_retval(&flag, "IDASetNonlinearSolver", 1)) return(1);

    /* In loop, call IDASolve, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached. */

    realtype tout;  // For output times

    int y_len_1 = y_len + 1; // remove later
    NumericMatrix soln(Dimension(time_vec_len,y_len_1));  // remove later


    // nothing to do with Events - single initialization of the ODE system
    // First row for initial conditions, First column is for time
    // int y_len_1 = y_len + 1;
    // NumericMatrix soln(Dimension(time_vec_len,y_len_1));

    // fill the first row of soln matrix with Initial Conditions
    soln(0,0) = time_vector[0];   // get the first time value
    for(int i = 0; i<y_len; i++){
      soln(0,i+1) = IC[i];
    }

    for(int iout = 0; iout < NOUT-1; iout++) {

      // output times start from the index after initial time
      tout = time_vector[iout+1];

      flag = IDASolve(ida_mem, tout, &time, yy0, yp0, IDA_NORMAL);

      if (check_retval(&flag, "CVode", 1)) {
        stop("Stopping IDA, something went wrong in solving the system of DAEs!"); break;
        } // Something went wrong in solving it!

      if (flag == IDA_SUCCESS) {

        // store results in soln matrix
        soln(iout+1, 0) = time;           // first column is for time
        for (int i = 0; i<y_len; i++){
          soln(iout+1, i+1) = yy0_ptr[i];
        }
      }
    }

    /* Free memory */
    IDAFree(&ida_mem);
    SUNNonlinSolFree(NLS);
    SUNLinSolFree(LS);
    SUNMatDestroy(SM);
    N_VDestroy(abstol);
    N_VDestroy(yy0);
    N_VDestroy(yp0);

    return soln;
    break;

  }

  default: {
    stop("Incorrect input function type - input function can be an R or Rcpp function");
  }

 }

}
