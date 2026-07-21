//   Copyright (c) 2016-2026, Satyaprakash Nayak
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
//   Neither sundialr nor the names of its
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
#include <nvector/nvector_serial.h>    /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h>   /* definition of type realtype */

#include <rhs_func.h>

using namespace Rcpp;

#define Ith(v,i)    NV_Ith_S(v,i-1)         /* i-th vector component i=1..NEQ */

// function called by CVodeInit if user inputs R function
//
// SUNDIALS calls this from its own C code, so an exception must not be allowed
// to escape: the body runs under sundials_callback_guard, which records the
// message and returns SUNDIALS' unrecoverable-error code instead. The solver
// then reports that message. See sundials_err_record.h.
int rhs_function(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data){

  // cast void pointer to pointer to struct and assign rhs to a Function
  struct rhs_func *my_rhs_fun = (struct rhs_func*)user_data;

  // nothing to record against, so just report the failure to SUNDIALS
  if(!my_rhs_fun){ return(-1); }

  return sundials_callback_guard(my_rhs_fun->err, [&]() -> int {

    // convert y to NumericVector y1
    int y_len = NV_LENGTH_S(y);

    NumericVector y1(y_len);    // filled with zeros
    for (int i = 0; i < y_len; i++){
      y1[i] = Ith(y,i+1);
    }

    NumericVector ydot1(y_len);    // filled with zeros

    Function rhs_fun = (*my_rhs_fun).rhs_eqn;           // function
    NumericVector p_values = (*my_rhs_fun).params;      // rate parameters

    if (!rhs_fun || (TYPEOF(rhs_fun) != CLOSXP)){
      stop("The RHS function is not an R or Rcpp function, stopping!");
    }

    // use the function to calculate value of RHS ----
    ydot1 = rhs_fun(t, y1, p_values);

    // The values below are copied out element by element, so a short return
    // would read past the end of ydot1 instead of being reported.
    if (ydot1.length() != y_len){
      stop("The RHS function must return a vector of the same length as the state vector: expected %d, got %d",
           y_len, ydot1.length());
    }

    // convert NumericVector ydot1 to N_Vector ydot
    sunrealtype *ydot_ptr = N_VGetArrayPointer(ydot);
    for (int i = 0; i<  y_len; i++){
      ydot_ptr[i] = ydot1[i];
    }

    // everything went smoothly
    return(0);
  });
}
//---RHS function definition ends ----------------------------------------------

