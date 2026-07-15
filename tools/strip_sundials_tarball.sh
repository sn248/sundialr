#! /bin/sh
# Strip an official SUNDIALS release tarball down to the subset needed to
# build the static libraries for sundialr, and repack it as
# sundials-mod-<version>.tar.gz (placed in the current directory).
#
# This removes examples, docs, tests, benchmarks, external bindings, and
# tooling to keep the R package tarball small (~2 MB vs ~25 MB upstream).
# Only the following are kept, mirroring the layout cmake needs:
#
#   cmake/  include/  src/  CMakeLists.txt
#   LICENSE  NOTICE  README.md  CHANGELOG.md  CITATIONS.md  CONTRIBUTING.md
#
# Usage:
#   tools/strip_sundials_tarball.sh <version> [path-to-sundials-<version>.tar.gz]
#
# If the tarball path is omitted, the official release is downloaded from
# GitHub. Example:
#   tools/strip_sundials_tarball.sh 7.8.0
#   mv sundials-mod-7.8.0.tar.gz src/

set -e

VERSION=$1
TARBALL=$2

if [ -z "${VERSION}" ]; then
    echo "Usage: $0 <version> [path-to-sundials-<version>.tar.gz]" >&2
    exit 1
fi

OUTDIR=$(pwd)
WORKDIR=$(mktemp -d)
trap 'rm -rf "${WORKDIR}"' EXIT

if [ -z "${TARBALL}" ]; then
    TARBALL="${WORKDIR}/sundials-${VERSION}.tar.gz"
    echo "Downloading sundials-${VERSION}.tar.gz from GitHub ..."
    curl -sL -o "${TARBALL}" \
        "https://github.com/LLNL/sundials/releases/download/v${VERSION}/sundials-${VERSION}.tar.gz"
fi

echo "Extracting ${TARBALL} ..."
tar -xzf "${TARBALL}" -C "${WORKDIR}"

SRCDIR="${WORKDIR}/sundials-${VERSION}"
if [ ! -d "${SRCDIR}" ]; then
    echo "Expected top-level directory sundials-${VERSION} not found in tarball" >&2
    exit 1
fi

# Delete everything not in the keep list.
KEEP="cmake include src CMakeLists.txt LICENSE NOTICE README.md CHANGELOG.md CITATIONS.md CONTRIBUTING.md"
for entry in "${SRCDIR}"/* "${SRCDIR}"/.[!.]*; do
    [ -e "${entry}" ] || continue
    name=$(basename "${entry}")
    keep_it=no
    for k in ${KEEP}; do
        if [ "${name}" = "${k}" ]; then keep_it=yes; break; fi
    done
    if [ "${keep_it}" = "no" ]; then
        echo "  removing ${name}"
        rm -rf "${entry}"
    fi
done

echo "Repacking sundials-mod-${VERSION}.tar.gz ..."
tar -czf "${OUTDIR}/sundials-mod-${VERSION}.tar.gz" -C "${WORKDIR}" "sundials-${VERSION}"

SIZE=$(du -h "${OUTDIR}/sundials-mod-${VERSION}.tar.gz" | cut -f1)
COUNT=$(tar -tzf "${OUTDIR}/sundials-mod-${VERSION}.tar.gz" | wc -l)
echo "Done: sundials-mod-${VERSION}.tar.gz (${SIZE}, ${COUNT} entries)"
