#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "SUNDIALS::core_static" for configuration "Release"
set_property(TARGET SUNDIALS::core_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::core_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_core.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::core_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::core_static "${_IMPORT_PREFIX}/lib/libsundials_core.a" )

# Import target "SUNDIALS::nvecserial_static" for configuration "Release"
set_property(TARGET SUNDIALS::nvecserial_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::nvecserial_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_nvecserial.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::nvecserial_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::nvecserial_static "${_IMPORT_PREFIX}/lib/libsundials_nvecserial.a" )

# Import target "SUNDIALS::nvecmanyvector_static" for configuration "Release"
set_property(TARGET SUNDIALS::nvecmanyvector_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::nvecmanyvector_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_nvecmanyvector.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::nvecmanyvector_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::nvecmanyvector_static "${_IMPORT_PREFIX}/lib/libsundials_nvecmanyvector.a" )

# Import target "SUNDIALS::sunmatrixband_static" for configuration "Release"
set_property(TARGET SUNDIALS::sunmatrixband_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::sunmatrixband_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_sunmatrixband.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunmatrixband_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunmatrixband_static "${_IMPORT_PREFIX}/lib/libsundials_sunmatrixband.a" )

# Import target "SUNDIALS::sunmatrixdense_static" for configuration "Release"
set_property(TARGET SUNDIALS::sunmatrixdense_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::sunmatrixdense_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_sunmatrixdense.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunmatrixdense_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunmatrixdense_static "${_IMPORT_PREFIX}/lib/libsundials_sunmatrixdense.a" )

# Import target "SUNDIALS::sunmatrixsparse_static" for configuration "Release"
set_property(TARGET SUNDIALS::sunmatrixsparse_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::sunmatrixsparse_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_sunmatrixsparse.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunmatrixsparse_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunmatrixsparse_static "${_IMPORT_PREFIX}/lib/libsundials_sunmatrixsparse.a" )

# Import target "SUNDIALS::sunlinsolband_static" for configuration "Release"
set_property(TARGET SUNDIALS::sunlinsolband_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::sunlinsolband_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_sunlinsolband.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsolband_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsolband_static "${_IMPORT_PREFIX}/lib/libsundials_sunlinsolband.a" )

# Import target "SUNDIALS::sunlinsoldense_static" for configuration "Release"
set_property(TARGET SUNDIALS::sunlinsoldense_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::sunlinsoldense_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_sunlinsoldense.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsoldense_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsoldense_static "${_IMPORT_PREFIX}/lib/libsundials_sunlinsoldense.a" )

# Import target "SUNDIALS::sunlinsolpcg_static" for configuration "Release"
set_property(TARGET SUNDIALS::sunlinsolpcg_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::sunlinsolpcg_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_sunlinsolpcg.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsolpcg_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsolpcg_static "${_IMPORT_PREFIX}/lib/libsundials_sunlinsolpcg.a" )

# Import target "SUNDIALS::sunlinsolspbcgs_static" for configuration "Release"
set_property(TARGET SUNDIALS::sunlinsolspbcgs_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::sunlinsolspbcgs_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_sunlinsolspbcgs.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsolspbcgs_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsolspbcgs_static "${_IMPORT_PREFIX}/lib/libsundials_sunlinsolspbcgs.a" )

# Import target "SUNDIALS::sunlinsolspfgmr_static" for configuration "Release"
set_property(TARGET SUNDIALS::sunlinsolspfgmr_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::sunlinsolspfgmr_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_sunlinsolspfgmr.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsolspfgmr_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsolspfgmr_static "${_IMPORT_PREFIX}/lib/libsundials_sunlinsolspfgmr.a" )

# Import target "SUNDIALS::sunlinsolspgmr_static" for configuration "Release"
set_property(TARGET SUNDIALS::sunlinsolspgmr_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::sunlinsolspgmr_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_sunlinsolspgmr.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsolspgmr_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsolspgmr_static "${_IMPORT_PREFIX}/lib/libsundials_sunlinsolspgmr.a" )

# Import target "SUNDIALS::sunlinsolsptfqmr_static" for configuration "Release"
set_property(TARGET SUNDIALS::sunlinsolsptfqmr_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::sunlinsolsptfqmr_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_sunlinsolsptfqmr.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsolsptfqmr_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsolsptfqmr_static "${_IMPORT_PREFIX}/lib/libsundials_sunlinsolsptfqmr.a" )

# Import target "SUNDIALS::sunnonlinsolnewton_static" for configuration "Release"
set_property(TARGET SUNDIALS::sunnonlinsolnewton_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::sunnonlinsolnewton_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_sunnonlinsolnewton.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunnonlinsolnewton_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunnonlinsolnewton_static "${_IMPORT_PREFIX}/lib/libsundials_sunnonlinsolnewton.a" )

# Import target "SUNDIALS::sunnonlinsolfixedpoint_static" for configuration "Release"
set_property(TARGET SUNDIALS::sunnonlinsolfixedpoint_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::sunnonlinsolfixedpoint_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_sunnonlinsolfixedpoint.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunnonlinsolfixedpoint_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunnonlinsolfixedpoint_static "${_IMPORT_PREFIX}/lib/libsundials_sunnonlinsolfixedpoint.a" )

# Import target "SUNDIALS::arkode_static" for configuration "Release"
set_property(TARGET SUNDIALS::arkode_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::arkode_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_arkode.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::arkode_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::arkode_static "${_IMPORT_PREFIX}/lib/libsundials_arkode.a" )

# Import target "SUNDIALS::cvode_static" for configuration "Release"
set_property(TARGET SUNDIALS::cvode_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::cvode_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_cvode.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::cvode_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::cvode_static "${_IMPORT_PREFIX}/lib/libsundials_cvode.a" )

# Import target "SUNDIALS::cvodes_static" for configuration "Release"
set_property(TARGET SUNDIALS::cvodes_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::cvodes_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_cvodes.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::cvodes_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::cvodes_static "${_IMPORT_PREFIX}/lib/libsundials_cvodes.a" )

# Import target "SUNDIALS::ida_static" for configuration "Release"
set_property(TARGET SUNDIALS::ida_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::ida_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_ida.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::ida_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::ida_static "${_IMPORT_PREFIX}/lib/libsundials_ida.a" )

# Import target "SUNDIALS::idas_static" for configuration "Release"
set_property(TARGET SUNDIALS::idas_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::idas_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_idas.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::idas_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::idas_static "${_IMPORT_PREFIX}/lib/libsundials_idas.a" )

# Import target "SUNDIALS::kinsol_static" for configuration "Release"
set_property(TARGET SUNDIALS::kinsol_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(SUNDIALS::kinsol_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsundials_kinsol.a"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::kinsol_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::kinsol_static "${_IMPORT_PREFIX}/lib/libsundials_kinsol.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
