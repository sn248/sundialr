# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2024, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was SUNDIALSConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

include(CMakeFindDependencyMacro)

### ------- Set FOUND status for SUNDIALS components

set(_installed_components "kinsol;idas;ida;cvodes;cvode;arkode;sunnonlinsolfixedpoint;sunnonlinsolnewton;sunlinsolsptfqmr;sunlinsolspgmr;sunlinsolspfgmr;sunlinsolspbcgs;sunlinsolpcg;sunlinsoldense;sunlinsolband;sunmatrixsparse;sunmatrixdense;sunmatrixband;nvecmanyvector;nvecserial;core;")

set(_comp_not_found "")
foreach(_comp ${SUNDIALS_FIND_COMPONENTS})
  if(_comp IN_LIST _installed_components)
    set(SUNDIALS_${_comp}_FOUND TRUE)
  else()
    set(SUNDIALS_${_comp}_FOUND FALSE)
    set(_comp_not_found "${_comp} ${_comp_not_found}")
  endif()
endforeach()

if(_comp_not_found)
  set(SUNDIALS_NOT_FOUND_MESSAGE "Component(s) not found: ${_comp_not_found}")
endif()

### ------- Import SUNDIALS targets

include("${CMAKE_CURRENT_LIST_DIR}/SUNDIALSTargets.cmake")

### ------- Alias targets

set(_SUNDIALS_ALIAS_TARGETS "sundials_kinsol->sundials_kinsol_static;sundials_idas->sundials_idas_static;sundials_ida->sundials_ida_static;sundials_cvodes->sundials_cvodes_static;sundials_cvode->sundials_cvode_static;sundials_arkode->sundials_arkode_static;sundials_sunnonlinsolfixedpoint->sundials_sunnonlinsolfixedpoint_static;sundials_sunnonlinsolnewton->sundials_sunnonlinsolnewton_static;sundials_sunlinsolsptfqmr->sundials_sunlinsolsptfqmr_static;sundials_sunlinsolspgmr->sundials_sunlinsolspgmr_static;sundials_sunlinsolspfgmr->sundials_sunlinsolspfgmr_static;sundials_sunlinsolspbcgs->sundials_sunlinsolspbcgs_static;sundials_sunlinsolpcg->sundials_sunlinsolpcg_static;sundials_sunlinsoldense->sundials_sunlinsoldense_static;sundials_sunlinsolband->sundials_sunlinsolband_static;sundials_sunmatrixsparse->sundials_sunmatrixsparse_static;sundials_sunmatrixdense->sundials_sunmatrixdense_static;sundials_sunmatrixband->sundials_sunmatrixband_static;sundials_nvecmanyvector->sundials_nvecmanyvector_static;sundials_nvecserial->sundials_nvecserial_static;sundials_core->sundials_core_static;")
foreach(ptr ${_SUNDIALS_ALIAS_TARGETS})
  string(REGEX REPLACE "sundials_" "" ptr "${ptr}")
  string(REGEX MATCHALL "([A-Za-z_]+)->([A-Za-z_]+)"
         _matches "${ptr}")
  set(_pointer ${CMAKE_MATCH_1})
  set(_pointee ${CMAKE_MATCH_2})
  if(NOT TARGET SUNDIALS::${_pointer})
    add_library(SUNDIALS::${_pointer} INTERFACE IMPORTED)
    target_link_libraries(SUNDIALS::${_pointer} INTERFACE SUNDIALS::${_pointee})
  endif()
endforeach()

### ------- Create TPL imported targets

if("OFF" AND NOT TARGET MPI::MPI_C)
  set(MPI_C_COMPILER "")
  find_dependency(MPI)
endif()

if("OFF" AND NOT TARGET OpenMP::OpenMP_C)
  find_dependency(OpenMP)
endif()

if("" AND NOT TARGET caliper)
  find_dependency(CALIPER PATHS "")
endif()

if("" AND NOT TARGET adiak::adiak)
  find_dependency(adiak PATHS "")
endif()

if("OFF" AND NOT (TARGET CUDA::cudart AND TARGET CUDA::cublas
   AND TARGET CUDA::cusparse AND TARGET CUDA::cusolver))
  find_dependency(CUDAToolkit)
endif()

if("OFF" AND NOT TARGET Ginkgo::ginkgo)
  if(NOT TARGET hwloc AND NOT HWLOC_DIR)
    set(HWLOC_DIR "")
  endif()
  find_dependency(Ginkgo PATHS "")
endif()

if("OFF" AND NOT TARGET SUNDIALS::HYPRE)
  add_library(SUNDIALS::HYPRE INTERFACE IMPORTED)
  target_link_libraries(SUNDIALS::HYPRE INTERFACE "")
  set_target_properties(SUNDIALS::HYPRE PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "")
endif()

if("OFF" AND NOT TARGET SUNDIALS::KLU)
  if("")
    find_dependency(KLU)
  else()
    add_library(SUNDIALS::KLU INTERFACE IMPORTED)
    target_link_libraries(SUNDIALS::KLU INTERFACE "")
    set_target_properties(SUNDIALS::KLU PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "")
  endif()
endif()

if("OFF" AND NOT TARGET Kokkos::kokkos)
  find_dependency(Kokkos PATHS "")
endif()

if("OFF" AND NOT TARGET Kokkos::kokkoskernels)
  find_dependency(KokkosKernels PATHS "")
endif()

if("OFF" AND NOT TARGET SUNDIALS::PETSC)
  add_library(SUNDIALS::PETSC INTERFACE IMPORTED)
  target_link_libraries(SUNDIALS::PETSC INTERFACE "")
  set_target_properties(SUNDIALS::PETSC PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "")

  # for backwards compatibility
  foreach(suffix SYS VEC MAT DM KSP SNES TS ALL)
    if(NOT TARGET SUNDIALS::PETSC_${suffix})
      add_library(SUNDIALS::PETSC_${suffix} INTERFACE IMPORTED)
      set_target_properties (SUNDIALS::PETSC_${suffix} PROPERTIES
        INTERFACE_LINK_LIBRARIES "SUNDIALS::PETSC")
    endif()
  endforeach()

  if("" MATCHES "Kokkos::kokkos")
    if(NOT TARGET Kokkos::kokkoskernels)
      find_dependency(KokkosKernels PATHS "")
    endif()
    if(NOT TARGET Kokkos::kokkos)
      find_dependency(Kokkos PATHS "")
    endif()
  endif()
endif()

if("OFF" AND NOT TARGET SUNDIALS::MAGMA)
  add_library(SUNDIALS::MAGMA INTERFACE IMPORTED)
  target_link_libraries(SUNDIALS::MAGMA INTERFACE "")
  set_target_properties(SUNDIALS::MAGMA PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "")
endif()

if("OFF" AND NOT TARGET MKL)
  find_dependency(MKL PATHS "")
endif()

if("OFF" AND NOT TARGET SUNDIALS::SUPERLUDIST)
  add_library(SUNDIALS::SUPERLUDIST INTERFACE IMPORTED)
  target_link_libraries(SUNDIALS::SUPERLUDIST INTERFACE "")
  set_target_properties(SUNDIALS::SUPERLUDIST PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "")
endif()

if("OFF" AND NOT TARGET SUNDIALS::SUPERLUMT)
  add_library(SUNDIALS::SUPERLUMT INTERFACE IMPORTED)
  target_link_libraries(SUNDIALS::SUPERLUMT INTERFACE "")
  set_target_properties(SUNDIALS::SUPERLUMT PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "")
endif()

if("OFF" AND NOT TARGET RAJA)
  if(NOT TARGET camp AND NOT camp_DIR)
    set(camp_DIR "")
  endif()
  find_dependency(RAJA PATHS "")
endif()

if("OFF" AND NOT TARGET SUNDIALS::TRILINOS)
  add_library(SUNDIALS::TRILINOS INTERFACE IMPORTED)
  target_link_libraries(SUNDIALS::TRILINOS INTERFACE "")
  set_target_properties(SUNDIALS::TRILINOS PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "")
endif()

if("OFF" AND NOT TARGET SUNDIALS::XBRAID)
  add_library(SUNDIALS::XBRAID INTERFACE IMPORTED)
  target_link_libraries(SUNDIALS::XBRAID INTERFACE "")
  set_target_properties(SUNDIALS::XBRAID PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "")
endif()

### ------- Check if required components were found

check_required_components(SUNDIALS)
