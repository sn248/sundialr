#! /bin/sh

: ${R_HOME=$(R RHOME)}
RSCRIPT_BIN=${R_HOME}/bin/Rscript
NCORES=`${RSCRIPT_BIN} -e "cat(min(2, parallel::detectCores(logical = FALSE)))"`

cd src

#### CMAKE CONFIGURATION ####
. ./scripts/cmake_config.sh

# compile sundials from source ###################################################
sh ./scripts/sundials_download.sh ${RSCRIPT_BIN}
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
    -D EXAMPLES_INSTALL=OFF \
  ${CMAKE_ADD_AR} ${CMAKE_ADD_RANLIB} ../sundials-src
make #-j${NCORES}
if [ $? -ne 0 ]; then
    echo "Make failed!"
    exit 1
fi

make install
if [ $? -ne 0 ]; then
    echo "Make install failed!"
    exit 1
fi
cd ..
##mv sundials/lib* sundials/lib
mv sundials-src/src/* ./sundials
rm -fr sundials-src sundials-build
