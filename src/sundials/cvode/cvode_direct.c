/*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *-----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *-----------------------------------------------------------------
 * Implementation file for the deprecated direct linear solver interface in
 * CVODE; these routines now just wrap the updated CVODE generic
 * linear solver interface in cvode_ls.h.
 *-----------------------------------------------------------------*/

#include <cvode/cvode_ls.h>
#include <cvode/cvode_direct.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*=================================================================
  Exported Functions (wrappers for equivalent routines in cvode_ls.h)
  =================================================================*/

int CVDlsSetLinearSolver(void *cvode_mem, SUNLinearSolver LS, SUNMatrix A)
{ return(CVodeSetLinearSolver(cvode_mem, LS, A)); }

int CVDlsSetJacFn(void *cvode_mem, CVDlsJacFn jac)
{ return(CVodeSetJacFn(cvode_mem, jac)); }

int CVDlsGetWorkSpace(void *cvode_mem, long int *lenrwLS, long int *leniwLS)
{ return(CVodeGetLinWorkSpace(cvode_mem, lenrwLS, leniwLS)); }

int CVDlsGetNumJacEvals(void *cvode_mem, long int *njevals)
{ return(CVodeGetNumJacEvals(cvode_mem, njevals)); }

int CVDlsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS)
{ return(CVodeGetNumLinRhsEvals(cvode_mem, nfevalsLS)); }

int CVDlsGetLastFlag(void *cvode_mem, long int *flag)
{ return(CVodeGetLastLinFlag(cvode_mem, flag)); }

char *CVDlsGetReturnFlagName(long int flag)
{ return(CVodeGetLinReturnFlagName(flag)); }


#ifdef __cplusplus
}
#endif
