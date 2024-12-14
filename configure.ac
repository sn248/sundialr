## setup configure.ac for the package
AC_PREREQ([2.70])
AC_INIT([sundialr],[0.1.5])
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
	if test x"${have_cmake}" == x"no"; then
		. src/scripts/cmake_config.sh
	fi
  tools/cmake_call.sh
  sundialr_include=""
  sundialr_libs="-lsundials_idas -lsundials_cvodes -lsundials_nvecserial -lsundials_sunlinsoldense -lsundials_sunmatrixdense -lsundials_core -lm"
fi
## Now use all the values
AC_SUBST([SUNDIALR_CPPFLAGS],["$sundialr_include"])
AC_SUBST([SUNDIALR_LIBS],["$sundialr_libs"])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
