#! /bin/sh

RSCRIPT_BIN=$1

## Uncompress sundials source
${RSCRIPT_BIN} -e "download.file(url = 'https://github.com/LLNL/sundials/releases/download/v7.1.1/sundials-7.1.1.tar.gz', destfile = 'sundials-7.1.1.tar.gz')"
${RSCRIPT_BIN} -e "utils::untar(tarfile = 'sundials-7.1.1.tar.gz')"
mv sundials-7.1.1 sundials-src
