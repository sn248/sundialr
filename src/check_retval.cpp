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

// --- Function definition for check_retval() ----------------------------------
/*
 * Check a SUNDIALS call's outcome. Two overloads, one per kind of return value:
 *
 *   check_retval(int flag,          const char*)  fails when flag < 0
 *                                                 (functions returning a code)
 *   check_retval(const void *ptr,   const char*)  fails when ptr == NULL
 *                                                 (allocating functions)
 *
 * Both return 1 on failure and 0 otherwise. Neither prints anything: the caller
 * raises the failure with sundials_stop(), which reports the message recorded by
 * the SUNDIALS error handler - naming the function, file, line and cause - so
 * printing here as well would duplicate that. funcname is kept in the signature
 * because it documents the call site.
 *
 * The compiler dispatches on the argument's type, so the two kinds of check can
 * no longer be confused. This replaced an earlier single function that took an
 * `int opt` selecting the check; opt had to match the TYPE of the value passed
 * and nothing enforced it, so a wrong pairing failed silently in both
 * directions - passing a pointer where a code was expected dereferenced it, and
 * passing &flag for the NULL check tested an address that is never NULL. Pass
 * the flag or the pointer itself; do not take its address or cast it.
 */

int check_retval(int flag, const char *funcname)
{
  /* An integer return code signals failure when negative. */
  if (flag < 0) { return(1); }
  return(0);
}

int check_retval(const void *returnvalue, const char *funcname)
{
  /* An allocating call signals failure by returning a NULL pointer. */
  if (returnvalue == NULL) { return(1); }
  return(0);
}

// --- Function definition ends ------------------------------------------------
