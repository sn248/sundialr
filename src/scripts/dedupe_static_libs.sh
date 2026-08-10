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
# This runs on the archives that configure.ac actually links; LINK_ORDER below
# must stay in sync with sundialr_libs there. A SUNDIALS upgrade can change which
# objects are bundled, so the verification guard at the end fails the build if
# any duplicate survives, rather than letting it resurface as a wasm-only error.
#
# Usage: sh dedupe_static_libs.sh <path-to-lib-dir>

LIBDIR=$1

if test -z "${LIBDIR}"; then
    echo "dedupe_static_libs.sh: no library directory given" >&2
    exit 1
fi

# AR/RANLIB come from scripts/r_config.sh (R CMD config), so a cross build gets
# its own toolchain - emar/emranlib under Emscripten. Plain ranlib cannot index
# wasm objects, so falling back to the host tools is a last resort only.
AR=${AR:-ar}
RANLIB=${RANLIB:-ranlib}

LINK_ORDER="libsundials_idas.a libsundials_cvodes.a libsundials_core.a"

REVERSED=""
for lib in ${LINK_ORDER}; do
    REVERSED="${lib} ${REVERSED}"
done

seen=""
for lib in ${REVERSED}; do
    if test ! -f "${LIBDIR}/${lib}"; then
        continue
    fi
    dupes=""
    for member in `${AR} t "${LIBDIR}/${lib}"`; do
        case " ${seen} " in
            *" ${member} "*) dupes="${dupes} ${member}" ;;
            *)               seen="${seen} ${member}"   ;;
        esac
    done
    if test -n "${dupes}"; then
        ${AR} d "${LIBDIR}/${lib}" ${dupes}
        if test $? -ne 0; then
            echo "dedupe_static_libs.sh: ${AR} d failed on ${lib}" >&2
            exit 1
        fi
        ${RANLIB} "${LIBDIR}/${lib}"
        if test $? -ne 0; then
            echo "dedupe_static_libs.sh: ${RANLIB} failed on ${lib}" >&2
            exit 1
        fi
        echo "dedupe_static_libs.sh: removed${dupes} from ${lib}"
    fi
done

# Verification guard: no member name may be defined by two linked archives.
all_members=""
for lib in ${LINK_ORDER}; do
    if test -f "${LIBDIR}/${lib}"; then
        all_members="${all_members} `${AR} t ${LIBDIR}/${lib}`"
    fi
done

surviving=`echo "${all_members}" | tr ' ' '\n' | grep -v '^$' | sort | uniq -d`
if test -n "${surviving}"; then
    echo "dedupe_static_libs.sh: duplicate members survive in ${LIBDIR}:" >&2
    echo "${surviving}" >&2
    echo "dedupe_static_libs.sh: a whole-archive link (webR) would fail" >&2
    exit 1
fi

exit 0
