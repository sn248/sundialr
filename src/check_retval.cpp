//   Copyright (c) 2024, Satyaprakash Nayak
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

// --- Function definition for check_flag() ------------------------------------
/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 *
 * This only reports whether the call failed; it does not print anything. The
 * caller raises the failure with sundials_stop(), which reports the message
 * recorded by the SUNDIALS error handler - naming the function, file, line and
 * cause - so printing here as well would duplicate that. funcname is kept in
 * the signature because it documents the call site.
 *
 * CAUTION when adding a call: opt must match the TYPE of the value passed, and
 * nothing enforces that. Getting it wrong fails silently, in both directions:
 *
 *   check_retval(&flag, "X", 0)      asks whether &flag is NULL. It never is,
 *                                    so a genuine failure is ignored.
 *   check_retval((void *)SM, "X", 1) reads the first bytes of a SUNMatrix as a
 *                                    return code, and dereferences NULL in
 *                                    exactly the case the check exists to catch.
 *
 * Both are one character away from the correct call. Every call site in this
 * package was audited and is correct: 31 pass &flag with opt = 1, 15 pass an
 * object pointer with opt = 0, and opt = 2 is never used. Keep it that way.
 *
 * The durable fix, if this is revisited, is to drop opt in favour of two
 * overloads - check_retval(int flag, const char*) and
 * check_retval(const void*, const char*) - so the compiler dispatches and a
 * mismatch stops compiling. It is deferred only because it touches every call
 * site and there is no live defect.
 */

int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) { return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) { return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) { return(1); }

  return(0);
}

// --- Function definition ends ------------------------------------------------
