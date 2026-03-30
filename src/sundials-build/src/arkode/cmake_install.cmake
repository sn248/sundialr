# Install script for directory: /home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/src/arkode

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
Install ARKODE
")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-build/src/arkode/libsundials_arkode.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/arkode" TYPE FILE FILES
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/arkode/arkode.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/arkode/arkode_arkstep.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/arkode/arkode_bandpre.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/arkode/arkode_bbdpre.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/arkode/arkode_butcher.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/arkode/arkode_butcher_dirk.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/arkode/arkode_butcher_erk.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/arkode/arkode_erkstep.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/arkode/arkode_forcingstep.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/arkode/arkode_ls.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/arkode/arkode_lsrkstep.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/arkode/arkode_mristep.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/arkode/arkode_splittingstep.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/arkode/arkode_sprk.h"
    "/home/quantilogy/Documents/GitHub/sn248/sundialr/src/sundials-src/include/arkode/arkode_sprkstep.h"
    )
endif()

