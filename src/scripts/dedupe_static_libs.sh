#!/bin/sh

# WebAssembly fix: SUNDIALS 7.x static libraries are self-contained. cmake bakes
# every entry of a module's OBJECT_LIBRARIES list (sunmemsys, nvecserial,
# sunmatrix*, sunlinsol*, sunnonlinsol*) into the archive itself - see
# $<TARGET_OBJECTS:...> in cmake/macros/SundialsAddLibrary.cmake - so
# libsundials_idas.a and libsundials_cvodes.a each carry their own private copy
# of nvector_serial.c.o, sundials_system_memory.c.o and friends.
#
# Native links tolerate that: ld extracts archive members on demand, the first
# definition wins, and the later copies are never pulled in. Emscripten links the
# package .so as a side module (-s SIDE_MODULE=1) and wraps every input in
# --whole-archive, because a side module's needed symbols are not known at link
# time. That force-includes all copies and they collide:
#
#   wasm-ld: error: duplicate symbol: N_VNewEmpty_Serial
#   >>> defined in ../inst/lib/libsundials_idas.a(nvector_serial.c.o)
#   >>> defined in ../inst/lib/libsundials_cvodes.a(nvector_serial.c.o)
#
# So delete each duplicated member from every archive except the LAST one that
# holds it in link order. Order matters and is not arbitrary: a later archive can
# still satisfy an earlier one in a single left-to-right pass, but not the
# reverse, so keeping the copy in the earlier archive instead would leave the
# later archive's own members with undefined references on native builds.
#
# THIS SCRIPT IS DELIBERATELY NON-FATAL. It always exits 0. An earlier version
# aborted the build when ar or ranlib failed, which broke Windows and macOS
# (commit 050c575, reverted in 5fa787c): it sits between `make install` and the
# `rm -fr sundials-src sundials-build` at the end of cmake_call.sh, so exiting
# non-zero here left both directories in the source tree, they were shipped
# inside the tarball, and the next run's cran_patches.sh died re-patching
# already-patched sources. Skipping the dedupe only leaves the wasm build broken
# - exactly as it is today - while every native platform links fine, because
# native ld ignores the duplicate members anyway. Degrading is always better than
# aborting here.
#
# It is also deliberately chatty: why ar/ranlib failed on Windows and macOS is
# still unknown, because rcmdcheck suppresses the build-phase output and the
# failure appeared in no log. Everything below is echoed so the next failure
# lands in 00install.out where it can be read.
#
# Usage: sh dedupe_static_libs.sh <path-to-lib-dir>

LIBDIR=$1
TAG="dedupe_static_libs.sh:"

warn_and_exit() {
    echo "${TAG} WARNING: $1"
    echo "${TAG} WARNING: skipping de-duplication; native builds are unaffected,"
    echo "${TAG} WARNING: the WebAssembly (webR) build will still fail to link."
    exit 0
}

if test -z "${LIBDIR}"; then
    warn_and_exit "no library directory given"
fi
if test ! -d "${LIBDIR}"; then
    warn_and_exit "library directory ${LIBDIR} does not exist"
fi

# AR/RANLIB come from scripts/r_config.sh (R CMD config), so a cross build gets
# its own toolchain - emar/emranlib under Emscripten. Plain ranlib cannot index
# wasm objects, so falling back to the host tools is a last resort only.
AR=${AR:-ar}
RANLIB=${RANLIB:-ranlib}

echo "${TAG} AR=${AR}"
echo "${TAG} RANLIB=${RANLIB}"

LINK_ORDER="libsundials_idas.a libsundials_cvodes.a libsundials_core.a"

REVERSED=""
for lib in ${LINK_ORDER}; do
    REVERSED="${lib} ${REVERSED}"
done

seen=""
for lib in ${REVERSED}; do
    if test ! -f "${LIBDIR}/${lib}"; then
        echo "${TAG} ${lib} not present, skipping"
        continue
    fi

    members=`${AR} t "${LIBDIR}/${lib}" 2>&1`
    if test $? -ne 0; then
        warn_and_exit "'${AR} t ${LIBDIR}/${lib}' failed: ${members}"
    fi

    dupes=""
    for member in ${members}; do
        case " ${seen} " in
            *" ${member} "*) dupes="${dupes} ${member}" ;;
            *)               seen="${seen} ${member}"   ;;
        esac
    done

    if test -z "${dupes}"; then
        echo "${TAG} ${lib}: no duplicated members"
        continue
    fi

    echo "${TAG} running: ${AR} d ${LIBDIR}/${lib}${dupes}"
    output=`${AR} d "${LIBDIR}/${lib}" ${dupes} 2>&1`
    if test $? -ne 0; then
        warn_and_exit "'${AR} d' failed on ${lib}: ${output}"
    fi

    echo "${TAG} running: ${RANLIB} ${LIBDIR}/${lib}"
    output=`${RANLIB} "${LIBDIR}/${lib}" 2>&1`
    if test $? -ne 0; then
        warn_and_exit "'${RANLIB}' failed on ${lib}: ${output}"
    fi

    echo "${TAG} ${lib}: removed${dupes}"
done

# Verification: no member name may be defined by two linked archives. A SUNDIALS
# upgrade can change which objects are bundled, so report rather than assume.
all_members=""
for lib in ${LINK_ORDER}; do
    if test -f "${LIBDIR}/${lib}"; then
        all_members="${all_members} `${AR} t ${LIBDIR}/${lib} 2>/dev/null`"
    fi
done

surviving=`echo "${all_members}" | tr ' ' '\n' | grep -v '^$' | sort | uniq -d`
if test -n "${surviving}"; then
    echo "${TAG} WARNING: duplicate members survive in ${LIBDIR}:"
    echo "${surviving}"
    warn_and_exit "de-duplication did not converge (SUNDIALS bundling may have changed)"
fi

echo "${TAG} OK: every object is defined exactly once across ${LINK_ORDER}"
exit 0
