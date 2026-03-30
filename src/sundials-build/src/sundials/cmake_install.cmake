# Install script for directory: /home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/src/sundials

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/quantilogy/Documents/GitHub/sn248/sundialr/inst")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  MESSAGE("
Install shared components
")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-build/src/sundials/libsundials_core.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/sundials" TYPE FILE FILES
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_adaptcontroller.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_band.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_base.hpp"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_context.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_context.hpp"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_convertibleto.hpp"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_core.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_core.hpp"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_dense.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_direct.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_errors.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_futils.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_iterative.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_linearsolver.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_linearsolver.hpp"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_logger.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_math.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_matrix.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_matrix.hpp"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_memory.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_memory.hpp"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_mpi_types.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_nonlinearsolver.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_nonlinearsolver.hpp"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_nvector.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_nvector.hpp"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_profiler.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_profiler.hpp"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_stepper.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_types_deprecated.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_types.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/sundials_version.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/sundials/priv" TYPE FILE FILES
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/priv/sundials_context_impl.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/sundials/priv/sundials_errors_impl.h"
    )
endif()

