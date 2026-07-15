#! /bin/sh
# CRAN-compliance patches for the bundled SUNDIALS source (v7.8.0).
#
# CRAN ("Writing R Extensions" §1.6.4) forbids compiled code that can call
# abort()/exit() or write to stdout/stderr. SUNDIALS errors are redirected to
# Rf_error() via a registered SUNContext_PushErrHandler (see
# inst/include/sundials_err_handler.h), which makes these code paths
# unreachable at runtime — but the symbols must not appear in the linked
# library, so they are patched out of the source before compilation.
#
# Usage: cran_patches.sh <path-to-extracted-sundials-source>
#   (run from src/; cmake_call.sh calls it with "sundials-src")
#
# NOTE: All patches use 'perl -pi -e' / python3 instead of 'sed -i' for
#   portability. macOS BSD sed requires 'sed -i ""' while GNU sed uses
#   'sed -i'; perl handles in-place editing identically on both platforms.
#
# Patches applied (verified against SUNDIALS 7.8.0 — re-verify each pattern
# on any SUNDIALS upgrade; the patterns match exact source lines):
#
# sundials_errors.c:
#   - SUNAbortErrHandlerFn: remove abort() (replaced by our Rf_error handler)
#   - SUNGlobalFallbackErrHandler: remove fprintf(stderr,...) calls
#   - sunCreateLogMessage: remove "FATAL LOGGER ERROR" fprintf(stderr,...)
# sundials_logger.c:
#   - sunCreateLogPayload: remove "FATAL LOGGER ERROR" fprintf(stderr,...)
#   - sunOpenLogFile: remove "stdout"/"stderr" filename-to-FILE* mappings
#     (new in 7.8.0: the logger has its own copy of the futils mapping)
#   - sunCloseLogFile: remove stdout/stderr comparisons (safe since no fp can
#     be stdout/stderr after the sunOpenLogFile patch)
#   - SUNLogger_Create: set error_fp/warning_fp/debug_fp/info_fp defaults to
#     NULL instead of stderr/stdout (7.8.0 sends logs to stdout by default)
# sundials_futils.c:
#   - SUNFileOpen: remove "stdout"/"stderr" filename-to-FILE* mappings
#   - SUNFileClose: remove stdout/stderr pointer comparisons before fclose
# sundials_nvector.c:
#   - N_VPrint: remove printf("NULL...\n") calls that GCC optimizes to puts()
# nvector_serial.c:
#   - N_VPrint_Serial: remove direct stdout argument to N_VPrintFile_Serial
# cvode/cvodes/ida/idas *_io.c, *_ls.c, *_diag.c:
#   - replace sprintf(name, "LITERAL") with strcpy(name, "LITERAL") to remove
#     the sprintf (___sprintf_chk on macOS) symbol from the linked libraries

set -e

SRC=${1:-sundials-src}

if [ ! -d "${SRC}/src/sundials" ]; then
    echo "cran_patches.sh: SUNDIALS source not found at ${SRC}" >&2
    exit 1
fi

## ---- sundials_errors.c -----------------------------------------------------

perl -pi -e 's|  abort\(\);|  /* CRAN: abort() removed; Rf_error raised via SUNContext_PushErrHandler */|' \
    "${SRC}/src/sundials/sundials_errors.c"

perl -pi -e 's|  fprintf\(stderr, "%s", log_msg\);|  /* CRAN: fprintf(stderr) removed */|g' \
    "${SRC}/src/sundials/sundials_errors.c"

## ---- "FATAL LOGGER ERROR" fprintf(stderr) block ----------------------------
# The same 3-line block appears in sundials_errors.c (sunCreateLogMessage)
# and sundials_logger.c (sunCreateLogPayload).

python3 -c "
import re, sys
for fname in ['${SRC}/src/sundials/sundials_errors.c',
              '${SRC}/src/sundials/sundials_logger.c']:
    with open(fname) as f: c = f.read()
    c, n = re.subn(
        r'    char\* fileAndLine = sunCombineFileAndLine\([^\n]+\);\n'
        r'    fprintf\(stderr,[^;]+FATAL LOGGER ERROR[^;]+;\n'
        r'    free\(fileAndLine\);',
        '    /* CRAN: fprintf(stderr) removed; fatal logger error silenced */',
        c)
    if n != 1:
        sys.exit('cran_patches.sh: FATAL LOGGER patch matched %d times in %s (expected 1)' % (n, fname))
    with open(fname, 'w') as f: f.write(c)
"

## ---- sundials_logger.c ------------------------------------------------------
# sunOpenLogFile: delete the stdout/stderr mapping branches, then strip the
# leading "else " from the fopen branch so no dangling else is left.

perl -ni -e 'print unless m{    if \(!strcmp\(fname, "stdout"\)\) \{ fp = stdout; \}}' \
    "${SRC}/src/sundials/sundials_logger.c"

perl -ni -e 'print unless m{    else if \(!strcmp\(fname, "stderr"\)\) \{ fp = stderr; \}}' \
    "${SRC}/src/sundials/sundials_logger.c"

perl -pi -e 's|    else \{ fp = fopen\(fname, mode\); \}|    fp = fopen(fname, mode); /* CRAN: stdout/stderr mappings removed */|' \
    "${SRC}/src/sundials/sundials_logger.c"

# sunCloseLogFile: safe since no fp can be stdout/stderr after the patch above
perl -pi -e 's|  if \(fp && fp != stdout && fp != stderr\) \{ fclose\(\(FILE\*\)fp\); \}|  if (fp) { fclose((FILE*)fp); } /* CRAN: stdout/stderr refs removed */|' \
    "${SRC}/src/sundials/sundials_logger.c"

# SUNLogger_Create: default output FILE* handles must not be stdout/stderr
perl -pi -e 's{logger->(\w+_fp)(\s*)= (?:stderr|stdout);}{logger->$1$2= NULL; /* CRAN: stdout/stderr not allowed */}' \
    "${SRC}/src/sundials/sundials_logger.c"

## ---- sundials_futils.c ------------------------------------------------------
# SUNFileOpen: same treatment as sunOpenLogFile (variable here is "filename")

perl -ni -e 'print unless m{    if \(!strcmp\(filename, "stdout"\)\) \{ fp = stdout; \}}' \
    "${SRC}/src/sundials/sundials_futils.c"

perl -ni -e 'print unless m{    else if \(!strcmp\(filename, "stderr"\)\) \{ fp = stderr; \}}' \
    "${SRC}/src/sundials/sundials_futils.c"

perl -pi -e 's|    else \{ fp = fopen\(filename, mode\); \}|    fp = fopen(filename, mode); /* CRAN: stdout/stderr mappings removed */|' \
    "${SRC}/src/sundials/sundials_futils.c"

perl -pi -e 's|  if \(fp && \(fp != stdout\) && \(fp != stderr\)\) \{ fclose\(fp\); \}|  if (fp) { fclose(fp); } /* CRAN: stdout/stderr refs removed */|' \
    "${SRC}/src/sundials/sundials_futils.c"

## ---- sundials_nvector.c -----------------------------------------------------
# N_VPrint uses printf("NULL...\n") which GCC -O2 optimises to puts().
# Replace with no-ops; these are debug-only print paths.

perl -pi -e 's|  if \(v == NULL\) \{ printf\("NULL Vector\\n"\); \}|  if (v == NULL) { /* CRAN: puts removed */ }|' \
    "${SRC}/src/sundials/sundials_nvector.c"

perl -pi -e 's|  else if \(v->ops->nvprint == NULL\) \{ printf\("NULL Print Op\\n"\); \}|  else if (v->ops->nvprint == NULL) { /* CRAN: puts removed */ }|' \
    "${SRC}/src/sundials/sundials_nvector.c"

## ---- nvector_serial.c -------------------------------------------------------
# N_VPrint_Serial passes stdout directly to N_VPrintFile_Serial. Replace with
# a no-op; users who need vector printing can call N_VPrintFile_Serial with an
# explicit FILE* handle.

perl -pi -e 's|  N_VPrintFile_Serial\(x, stdout\);|  /* CRAN: stdout removed; use N_VPrintFile_Serial with explicit FILE* */|' \
    "${SRC}/src/nvector/serial/nvector_serial.c"

## ---- sprintf -> strcpy ------------------------------------------------------
# All sprintf calls in these files use a plain string literal with no format
# specifiers, so strcpy is a safe and equivalent replacement. cvode/ida are
# included as well as cvodes/idas even though only the latter are linked, so
# the compiled sources committed to src/sundials/ are uniformly clean.

for f in \
    "${SRC}/src/cvode/cvode_io.c" \
    "${SRC}/src/cvode/cvode_ls.c" \
    "${SRC}/src/cvode/cvode_diag.c" \
    "${SRC}/src/cvodes/cvodes_io.c" \
    "${SRC}/src/cvodes/cvodes_ls.c" \
    "${SRC}/src/cvodes/cvodes_diag.c" \
    "${SRC}/src/ida/ida_io.c" \
    "${SRC}/src/ida/ida_ls.c" \
    "${SRC}/src/idas/idas_io.c" \
    "${SRC}/src/idas/idas_ls.c"; do
  perl -pi -e 's/sprintf\((\w+), ("(?:[^"\\]|\\.)*")\)/strcpy($1, $2)/g' "$f"
done

## ---- newline termination ----------------------------------------------------
# R CMD check notes C sources/headers not terminated with a newline (e.g.
# sundomeigest_arnoldi.c in SUNDIALS 7.8.0). Append one where missing.

find "${SRC}/src" \( -name '*.c' -o -name '*.h' \) -print | while IFS= read -r f; do
  if [ -n "$(tail -c1 "$f")" ]; then echo >> "$f"; fi
done

## ---- verification guard -----------------------------------------------------
# Fail loudly if any CRAN-flagged call survives in sources compiled into the
# libraries linked by sundialr (core, cvodes, idas, nvecserial,
# sunlinsoldense, sunmatrixdense) plus the patched cvode/ida sources.
# Excluded: fmod_* dirs (Fortran interfaces, not compiled) and
# sundials_profiler.c (its printf is inside #if SUNDIALS_MPI_ENABLED, off).

GUARD_DIRS="${SRC}/src/sundials ${SRC}/src/cvode ${SRC}/src/cvodes \
${SRC}/src/ida ${SRC}/src/idas ${SRC}/src/nvector/serial \
${SRC}/src/sunlinsol/dense ${SRC}/src/sunmatrix/dense"

guard_fail=0
for pat in \
    'abort *\(' \
    'fprintf *\( *stderr' \
    'fprintf *\( *stdout' \
    '(^|[^a-zA-Z_])printf *\(' \
    '(^|[^a-zA-Z_n])sprintf *\('; do
  hits=`grep -rEn --include='*.c' --exclude-dir=fmod_int32 --exclude-dir=fmod_int64 \
        "$pat" ${GUARD_DIRS} 2>/dev/null | grep -v 'CRAN:' | grep -v 'sundials_profiler\.c' || true`
  if [ -n "${hits}" ]; then
    echo "cran_patches.sh: CRAN-flagged pattern '$pat' survives patching:" >&2
    echo "${hits}" >&2
    guard_fail=1
  fi
done

if [ ${guard_fail} -ne 0 ]; then
    echo "cran_patches.sh: patch patterns need re-anchoring against this SUNDIALS version" >&2
    exit 1
fi

echo "cran_patches.sh: all CRAN patches applied and verified"
