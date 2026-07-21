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

##mv sundials/lib* sundials/lib
mv sundials-src/src/* ./sundials
rm -fr sundials-src sundials-build


