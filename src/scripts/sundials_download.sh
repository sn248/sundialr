#! /bin/sh

RSCRIPT_BIN=$1

## Uncompress sundials source
# ${RSCRIPT_BIN} -e "download.file(url = 'https://github.com/LLNL/sundials/releases/download/v7.2.0/sundials-7.2.0.tar.gz', destfile = 'sundials-7.2.0.tar.gz')"
#
# if [ $? -ne 0 ]; then
#     echo "Downloading sundials from github failed!"
#     exit 1
# fi
#
# ${RSCRIPT_BIN} -e "utils::untar(tarfile = 'sundials-7.2.0.tar.gz')"
# mv sundials-7.2.0 sundials-src


# Remove any tree left by an earlier configure run before extracting. Without
# this, the mv below moves the fresh sundials-7.8.0 *inside* the existing
# sundials-src instead of replacing it, so the stale (already CRAN-patched)
# sources survive; cran_patches.sh then re-patches them, its verification guard
# sees the patterns match zero times and the build dies with a message that
# points at the patch script rather than at the leftover directory.
#
# R CMD build installs the package to build vignettes, which runs configure in
# the source directory, so a failure anywhere in cmake_call.sh leaves exactly
# that state behind - and .Rbuildignore now keeps it out of the tarball too.
rm -fr sundials-src sundials-7.8.0

${RSCRIPT_BIN} -e "utils::untar(tarfile = 'sundials-mod-7.8.0.tar.gz')"
if [ $? -ne 0 ]; then
    echo "Could not extract the sundials tar file"
    exit 1
fi
#
mv sundials-7.8.0 sundials-src
