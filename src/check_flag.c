#include <Rcpp.h>

// --- Function definition for check_flag() ------------------------------------
int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL)
  {
    // fprintf(stderr, "\n SUNDIALS_ERROR: %s() failed - returned NULL pointer \n\n", funcname);
    Rprintf("SUNDIALS_ERROR: %s() failed - returned NULL pointer", funcname);
    return(1);
  }

  /* Check if flag < 0 */
  else if (opt == 1)
  {
    errflag = (int *) flagvalue;
    if (*errflag < 0)
    {
      //fprintf(stderr, "\n SUNDIALS_ERROR: %s() failed with flag = %d \n\n", funcname, *errflag);
      Rprintf("SUNDIALS_ERROR: %s() failed with flag = %d", funcname, *errflag);
      return(1);
    }
  }
  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    // fprintf(stderr, "\n MEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    Rprintf("MEMORY_ERROR: %s() failed - returned NULL pointer", funcname);
    return (1);
  }

  return (0); // so returns 0 if everything is okay

}
// --- Function definition ends ------------------------------------------------
