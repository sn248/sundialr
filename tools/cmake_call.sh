#! /bin/sh

: ${R_HOME=$(R RHOME)}
RSCRIPT_BIN=${R_HOME}/bin/Rscript
NCORES=`${RSCRIPT_BIN} -e "cat(min(2, parallel::detectCores(logical = FALSE)))"`

cd src

#### CMAKE CONFIGURATION ####
. ./scripts/cmake_config.sh

# compile sundials from source ###################################################
sh ./scripts/sundials_download.sh ${RSCRIPT_BIN}

# CRAN fix: patch SUNDIALS source to remove abort(), fprintf(stderr,...), and
# stdout/stderr symbol references, which CRAN flags as non-compliant.
# See "Writing R Extensions" §1.6.4. Our registered R error handler (via
# SUNContext_PushErrHandler) makes these code paths unreachable at runtime.
#
# sundials_errors.c:
#   - SUNAbortErrHandlerFn: remove abort() (replaced by our Rf_error handler)
#   - SUNGlobalFallbackErrHandler: remove fprintf(stderr,...) calls
# sundials_logger.c:
#   - SUNLogger_Create: set error_fp/warning_fp to NULL instead of stderr/stdout
#   - sunCloseLogFile: remove stdout/stderr comparisons (safe since fps are NULL)
#   - remove fprintf(stderr,...) fatal logger error message
# sundials_futils.c:
#   - SUNDIALSFileOpen: remove "stdout"/"stderr" filename-to-FILE* mappings
#   - SUNDIALSFileClose: remove stdout/stderr pointer comparisons before fclose
# sundials_nvector.c:
#   - N_VPrint: remove printf("NULL...\n") calls that GCC optimizes to puts()
# nvector_serial.c:
#   - N_VPrint_Serial: remove direct stdout argument to N_VPrintFile_Serial

sed -i 's|  abort();|  /* CRAN: abort() removed; Rf_error raised via SUNContext_PushErrHandler */|' \
    sundials-src/src/sundials/sundials_errors.c

sed -i 's|  fprintf(stderr, "%s", log_msg);|  /* CRAN: fprintf(stderr) removed */|g' \
    sundials-src/src/sundials/sundials_errors.c

sed -i 's|  logger->error_fp   = stderr;|  logger->error_fp   = NULL; /* CRAN: stderr not allowed */|' \
    sundials-src/src/sundials/sundials_logger.c

sed -i 's|  logger->warning_fp = stdout;|  logger->warning_fp = NULL; /* CRAN: stdout not allowed */|' \
    sundials-src/src/sundials/sundials_logger.c

sed -i 's|  if (fp && fp != stdout && fp != stderr)|  if (fp) /* CRAN: stdout/stderr refs removed */|' \
    sundials-src/src/sundials/sundials_logger.c

sed -i 's|    fprintf(stderr, "\[FATAL LOGGER ERROR\].*|    /* CRAN: fprintf(stderr) removed */|' \
    sundials-src/src/sundials/sundials_logger.c

# sundials_futils.c: SUNDIALSFileOpen maps "stdout"/"stderr" names to FILE* handles.
# Remove those branches; callers passing those names will get NULL (no output).
sed -i 's|    if (!strcmp(filename, "stdout")) { fp = stdout; }|    /* CRAN: stdout mapping removed */|' \
    sundials-src/src/sundials/sundials_futils.c

sed -i 's|    else if (!strcmp(filename, "stderr")) { fp = stderr; }|    /* CRAN: stderr mapping removed */|' \
    sundials-src/src/sundials/sundials_futils.c

# SUNDIALSFileClose guards fclose with fp != stdout and fp != stderr comparisons.
# Safe to remove since we no longer assign stdout/stderr to any fp above.
sed -i 's|  if (fp && (fp != stdout) && (fp != stderr)) { fclose(fp); }|  if (fp) { fclose(fp); } /* CRAN: stdout/stderr refs removed */|' \
    sundials-src/src/sundials/sundials_futils.c

# sundials_nvector.c: N_VPrint uses printf("NULL...\n") which GCC -O2 optimises
# to puts(). Replace with no-ops; these are debug-only print paths.
sed -i 's|  if (v == NULL) { printf("NULL Vector\\n"); }|  if (v == NULL) { /* CRAN: puts removed */ }|' \
    sundials-src/src/sundials/sundials_nvector.c

sed -i 's|  else if (v->ops->nvprint == NULL) { printf("NULL Print Op\\n"); }|  else if (v->ops->nvprint == NULL) { /* CRAN: puts removed */ }|' \
    sundials-src/src/sundials/sundials_nvector.c

# nvector_serial.c: N_VPrint_Serial passes stdout directly to N_VPrintFile_Serial.
# Replace with a no-op; users who need vector printing can call N_VPrintFile_Serial
# with an explicit FILE* handle.
sed -i 's|  N_VPrintFile_Serial(x, stdout);|  /* CRAN: stdout removed; use N_VPrintFile_Serial with explicit FILE* */|' \
    sundials-src/src/nvector/serial/nvector_serial.c

dot() { file=$1; shift; . "$file"; }
dot ./scripts/r_config.sh ""
${RSCRIPT_BIN} --vanilla -e 'getRversion() > "4.0.0"' | grep TRUE > /dev/null
if [ $? -eq 0 ]; then
  CMAKE_ADD_AR="-D CMAKE_AR=${AR}"
  CMAKE_ADD_RANLIB="-D CMAKE_RANLIB=${RANLIB}"
else
  CMAKE_ADD_AR=""
  CMAKE_ADD_RANLIB=""
fi
mkdir sundials-build
mkdir sundials
cd sundials-build
${CMAKE_BIN} \
    -D CMAKE_BUILD_TYPE=Release \
    -D BUILD_STATIC_LIBS=ON \
    -D BUILD_SHARED_LIBS=OFF \
    -D CMAKE_INSTALL_PREFIX=../../inst \
    -D BUILD_TESTING=OFF \
    -D CMAKE_C_STANDARD=11 \
    -D EXAMPLES_ENABLE_C=OFF \
    -D EXAMPLES_ENABLE_CXX=OFF \
    -D EXAMPLES_INSTALL=OFF \
    -D SUNDIALS_LOGGING_LEVEL=0 \
  ${CMAKE_ADD_AR} ${CMAKE_ADD_RANLIB} ../sundials-src
  # CRAN fix: SUNDIALS_LOGGING_LEVEL=0 disables all stdout/stderr logging output
  # which CRAN flags as non-compliant (compiled code must not write to stdout/stderr)
make -j${NCORES}
make install
cd ..
##mv sundials/lib* sundials/lib
mv sundials-src/src/* ./sundials
rm -fr sundials-src sundials-build
