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


${RSCRIPT_BIN} -e "utils::untar(tarfile = 'sundials-mod-7.2.1.tar.gz')"
if [ $? -ne 0 ]; then
    echo "Could not extract the sundials tar file"
    exit 1
fi
#
mv sundials-7.2.1 sundials-src
