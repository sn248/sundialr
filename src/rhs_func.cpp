#include <Rcpp.h>
#include <nvector/nvector_serial.h>    /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h>   /* definition of type realtype */

#include <rhs_func.h>

using namespace Rcpp;

#define Ith(v,i)    NV_Ith_S(v,i-1)         /* i-th vector component i=1..NEQ */

// function called by CVodeInit if user inputs R function
int rhs_function(realtype t, N_Vector y, N_Vector ydot, void* user_data){

  // convert y to NumericVector y1
  int y_len = NV_LENGTH_S(y);

  NumericVector y1(y_len);    // filled with zeros
  // realtype *y_ptr = N_VGetArrayPointer(y);
  for (int i = 0; i < y_len; i++){
    y1[i] = Ith(y,i+1); // y_ptr[i];
  }

  NumericVector ydot1(y_len);    // filled with zeros

  if(!user_data){
    stop("Something went wrong, stopping!");
  }

  // cast void pointer to pointer to struct and assign rhs to a Function
  struct rhs_func *my_rhs_fun = (struct rhs_func*)user_data;

  if(my_rhs_fun){
    Function rhs_fun = (*my_rhs_fun).rhs_eqn;           // function
    NumericVector p_values = (*my_rhs_fun).params;      // rate parameters

    if (rhs_fun && (TYPEOF(rhs_fun) == CLOSXP)){
      // use the function to calculate value of RHS ----
      ydot1 = rhs_fun(t, y1, p_values);
    }
    else{
      stop("Something went wrong, stopping!");
    }
  }
  else {
    stop("Something went wrong, stopping!");
  }

  // convert NumericVector ydot1 to N_Vector ydot
  realtype *ydot_ptr = N_VGetArrayPointer(ydot);
  for (int i = 0; i<  y_len; i++){
    ydot_ptr[i] = ydot1[i];
  }

  // everything went smoothly
  return(0);
}
//---RHS function definition ends ----------------------------------------------

