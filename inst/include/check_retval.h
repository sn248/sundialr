// File: check_flag.h

#ifndef CHECK_RETVAL_H
#define CHECK_RETVAL_H

// Two overloads so the compiler dispatches on the argument's type and a
// mismatched call stops compiling. See src/check_retval.cpp for the rationale.
//   - integer return codes: fails when flag < 0
//   - allocating calls returning a pointer: fails when the pointer is NULL
int check_retval(int flag, const char *funcname);
int check_retval(const void *returnvalue, const char *funcname);

#endif /* CHECK_RETVAL */
