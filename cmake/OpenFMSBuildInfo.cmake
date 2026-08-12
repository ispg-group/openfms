function(openfms_build_info_environment output_variable)
  # Build information environment variables
  set(_openfms_fortran_flags)
  openfms_collect_fortran_flags(_openfms_fortran_flags)

  set(_openfms_environment
    "OPENFMS_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
    "OPENFMS_VERSION=${OPENFMS_VERSION}"
    "OPENFMS_BACKEND=${OPENFMS_BACKEND}"
    "OPENFMS_BUILD_SYSTEM=CMake"
    "OPENFMS_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
    "OPENFMS_FC=${CMAKE_Fortran_COMPILER}"
    "OPENFMS_FC_ID=${CMAKE_Fortran_COMPILER_ID}"
    "OPENFMS_FC_VERSION=${CMAKE_Fortran_COMPILER_VERSION}"
    "OPENFMS_FFLAGS=${_openfms_fortran_flags}"
  )

  set("${output_variable}" ${_openfms_environment} PARENT_SCOPE)
endfunction()

function(openfms_add_build_info output_variable)
  # Directs generator script to create a Fortran source file (build_info.F90)
  # that can print the build information
  set(_openfms_build_info_source "${CMAKE_CURRENT_BINARY_DIR}/generated/build_info.F90")
  set(_openfms_build_info_stamp "${_openfms_build_info_source}.stamp")
  set(_openfms_build_info_script "${CMAKE_CURRENT_SOURCE_DIR}/generate_build_info.sh")
  openfms_build_info_environment(_openfms_build_info_environment)

  # Rerun build command when dependencies are newer than stamp 
  add_custom_command(
    OUTPUT "${_openfms_build_info_stamp}"
    BYPRODUCTS "${_openfms_build_info_source}"
    COMMAND
      "${CMAKE_COMMAND}" -E env
      ${_openfms_build_info_environment}
      "${_openfms_build_info_script}"
      "${_openfms_build_info_source}"
    COMMAND
      "${CMAKE_COMMAND}" -E touch
      "${_openfms_build_info_stamp}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    DEPENDS "${_openfms_build_info_script}"
    COMMENT "Generating OpenFMS build information"
    VERBATIM
  )

  add_custom_target(openfms_build_info
    DEPENDS "${_openfms_build_info_stamp}"
  )

  set_source_files_properties("${_openfms_build_info_source}" PROPERTIES
    GENERATED TRUE
  )

  set("${output_variable}" "${_openfms_build_info_source}" PARENT_SCOPE)
endfunction()

function(openfms_refresh_build_info_with source)
  add_custom_command(
    OUTPUT "${source}.stamp"
    APPEND
    DEPENDS ${ARGN}
  )
endfunction()

function(openfms_configure_build_info_source source)
  string(TIMESTAMP _openfms_build_date "%Y-%m-%d %H:%M:%S %z")
  set_source_files_properties("${source}" PROPERTIES
    COMPILE_DEFINITIONS "OPENFMS_BUILD_DATE=\"${_openfms_build_date}\""
  )
endfunction()
