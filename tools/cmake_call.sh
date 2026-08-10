#! /bin/sh

: ${R_HOME=$(R RHOME)}
RSCRIPT_BIN=${R_HOME}/bin/Rscript
NCORES=`${RSCRIPT_BIN} -e "cat(min(2, parallel::detectCores(logical = FALSE)))"`

cd src

#### CMAKE CONFIGURATION ####
. ./scripts/cmake_config.sh

# compile sundials from source ###################################################
sh ./scripts/sundials_download.sh ${RSCRIPT_BIN}

# CRAN fix: patch SUNDIALS source to remove abort(), fprintf(stderr,...),
# stdout/stderr symbol references, and sprintf, which CRAN flags as
# non-compliant ("Writing R Extensions" §1.6.4). The patches live in
# scripts/cran_patches.sh, which also verifies that no flagged call survives
# and fails the build otherwise (e.g. when a SUNDIALS upgrade shifts the
# patched lines).
sh ./scripts/cran_patches.sh sundials-src
if [ $? -ne 0 ]; then
    echo "Applying CRAN patches to SUNDIALS source failed!"
    exit 1
fi

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
    -D CMAKE_INSTALL_LIBDIR=lib \
    -D CMAKE_INSTALL_PREFIX=../../inst \
    -D BUILD_TESTING=OFF \
    -D CMAKE_C_STANDARD=99 \
    -D EXAMPLES_ENABLE_C=OFF \
    -D EXAMPLES_ENABLE_CXX=OFF \
    -D SUNDIALS_LOGGING_LEVEL=0 \
    -D SUNDIALS_ENABLE_ARKODE=OFF \
    -D SUNDIALS_ENABLE_KINSOL=OFF \
    -D CMAKE_C_FLAGS="${CFLAGS} -Wno-deprecated-declarations" \
  ${CMAKE_ADD_AR} ${CMAKE_ADD_RANLIB} ../sundials-src
  # CRAN fixes:
  #   SUNDIALS_LOGGING_LEVEL=0 disables all stdout/stderr logging output
  #   -Wno-deprecated-declarations suppresses N_VSpace/SUNMatSpace deprecation
  #   warnings from SUNDIALS 7.x (these functions are removed in 8.0 but still
  #   functional; CRAN treats installation warnings as errors)
make -j${NCORES}
make install
if [ $? -ne 0 ]; then
    echo "Make install failed!"
    exit 1
fi
cd ..

# WebAssembly fix: SUNDIALS 7.x bakes the shared object libraries (nvecserial,
# sunmatrix*, sunlinsol*, ...) into every package archive, so libsundials_idas.a
# and libsundials_cvodes.a ship duplicate copies of the same objects. Native ld
# extracts archive members on demand and ignores the extras, but Emscripten links
# the package .so as a side module with --whole-archive and every copy collides
# ("wasm-ld: error: duplicate symbol: N_VNewEmpty_Serial"). Strip the duplicates
# so each object is defined exactly once.
#
# Its exit status is deliberately NOT checked, and the script always exits 0
# anyway: it runs after `make install` but before the cleanup at the end of this
# file, so aborting here strands sundials-src/ and sundials-build/ in the source
# tree, which then ship inside the tarball and break the *next* build in a way
# that points at cran_patches.sh (this is what commit 050c575 did to Windows and
# macOS). A skipped dedupe leaves only the wasm build broken.
#
# AR/RANLIB are set by scripts/r_config.sh above but not exported, so pass them
# through explicitly - a cross build needs its own archiver (emar/emranlib), as
# host ranlib cannot index wasm objects.
AR="${AR}" RANLIB="${RANLIB}" sh ./scripts/dedupe_static_libs.sh ../inst/lib

# Reproducibility fix: cmake stamps the generated sundials_config.h with the
# compiler, its full flag list and a build timestamp. Those differ on every
# machine and every run, so this header - which is committed, like the rest of
# inst/include, because LinkingTo consumers need it - would otherwise carry one
# machine's build paths and show up as a spurious diff after every build. The
# strings are provenance metadata only (they feed SUNDIALS' optional Caliper
# instrumentation, not compilation), so blanking them makes the generated
# header identical from run to run.
SUNDIALS_CONFIG_H="../inst/include/sundials/sundials_config.h"
if [ -f "${SUNDIALS_CONFIG_H}" ]; then
    perl -i -pe 's/^(#define SUN_(?:C|CXX|FORTRAN)_COMPILER(?:_VERSION|_FLAGS)?)\s+".*"\s*$/$1 ""\n/;
                 s/^(#define SUN_JOB_(?:ID|START_TIME))\s+".*"\s*$/$1 ""\n/;' "${SUNDIALS_CONFIG_H}"
    if [ $? -ne 0 ]; then
        echo "Normalising sundials_config.h failed!"
        exit 1
    fi
fi

# Line-ending fix: sundials_config.h and sundials_export.h are the only two
# headers cmake *generates*. The other 42 in inst/include/sundials are copied
# verbatim out of the bundled tarball and so keep its LF endings, which is why
# these two alone came out CRLF on Windows and were reported by
#   checking line endings in C/C++/Fortran sources/headers ... NOTE
# Force both to LF so the installed headers are byte-identical whatever host
# built them, as the reproducibility fix above already requires.
#
# -Mopen=IO,:raw matters: without it Strawberry Perl applies its default :crlf
# layer to both ends, stripping the CR on read and adding it back on write, so
# the substitution would silently do nothing on the one platform that needs it.
for SUNDIALS_GEN_H in "../inst/include/sundials/sundials_config.h" \
                      "../inst/include/sundials/sundials_export.h"; do
    if [ -f "${SUNDIALS_GEN_H}" ]; then
        perl -i -Mopen=IO,:raw -pe 's/\r\n/\n/g' "${SUNDIALS_GEN_H}"
        if [ $? -ne 0 ]; then
            echo "Normalising line endings in ${SUNDIALS_GEN_H} failed!"
            exit 1
        fi
    fi
done

##mv sundials/lib* sundials/lib
mv sundials-src/src/* ./sundials
rm -fr sundials-src sundials-build


