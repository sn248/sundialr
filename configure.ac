## setup configure.ac for the package
AC_PREREQ([2.70])
AC_INIT([sundialr],[0.2.0])
## Set R_HOME, respecting an environment variable if set
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
	AC_MSG_ERROR([Could not determine R_HOME.])
fi
CXX=`"${R_HOME}/bin/R" CMD config CXX`
if test -z "${CXX}"; then
	AC_MSG_ERROR([No C++ compiler available])
fi
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
CXX11FLAGS=`"${R_HOME}/bin/R" CMD config CXX11FLAGS`
AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX
## Default to build from source
need_to_build="yes"
## So do we need to build
if test x"${need_to_build}" != x"no"; then
      	AC_PATH_PROG(have_cmake, cmake, no)
	# CRAN fix: use POSIX-compliant single = in test; == is a bashism
	if test x"${have_cmake}" = x"no"; then
		. src/scripts/cmake_config.sh
	fi
  tools/cmake_call.sh
  ## cmake_call.sh exits non-zero on every failure it detects, but its status was
  ## previously ignored, so a failed SUNDIALS build let configure report success
  ## and surfaced much later as "cannot find -lsundials_idas" at link time - which
  ## reads like a bad link line rather than a build that never happened. Fail here
  ## instead, where the real error is still on screen.
  if test $? -ne 0; then
    AC_MSG_ERROR([building the bundled SUNDIALS libraries failed; see the output above])
  fi
  sundialr_include=""
  ## WebAssembly fix: link only the self-contained archives. SUNDIALS 7.x bakes
  ## the shared object libraries (nvecserial, sunmatrix*, sunlinsol*, ...) into
  ## every package archive, so libsundials_cvodes.a already provides nvecserial,
  ## sunlinsoldense and sunmatrixdense. Naming them again added a second copy of
  ## each object to the link, which only whole-archive linkers reject. Keep this
  ## in sync with LINK_ORDER in src/scripts/dedupe_static_libs.sh.
  sundialr_libs="-lsundials_idas -lsundials_cvodes -lsundials_core -lm"
  ## tools/remove_static_libs.sh
fi
## Now use all the values
AC_SUBST([SUNDIALR_CPPFLAGS],["$sundialr_include"])
AC_SUBST([SUNDIALR_LIBS],["$sundialr_libs"])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
