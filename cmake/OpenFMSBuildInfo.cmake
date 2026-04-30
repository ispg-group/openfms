function(openfms_build_info_environment output_variable)
  # The generated Fortran build-info source prints the same settings that were
  # used for this CMake build. Collect them as environment variables so the
  # existing shell generator can stay build-system agnostic.
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
  # OpenFMS reports version, compiler, backend, and git state at runtime. Those
  # values are not known until configure/build time, so CMake asks the project
  # generator script to create a Fortran source file in the build tree.
  set(_openfms_build_info_source "${CMAKE_CURRENT_BINARY_DIR}/generated/build_info.F90")
  set(_openfms_build_info_script "${CMAKE_CURRENT_SOURCE_DIR}/generate_build_info.sh")
  openfms_build_info_environment(_openfms_build_info_environment)

  # The custom target gives the executable a concrete dependency on the
  # generated source, while the command itself avoids rewriting the file when
  # the generated contents are unchanged.
  add_custom_command(
    OUTPUT "${_openfms_build_info_source}"
    COMMAND
      "${CMAKE_COMMAND}" -E env
      ${_openfms_build_info_environment}
      "${_openfms_build_info_script}"
      "${_openfms_build_info_source}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    DEPENDS "${_openfms_build_info_script}"
    COMMENT "Generating OpenFMS build information"
    VERBATIM
  )

  add_custom_target(openfms_build_info
    DEPENDS "${_openfms_build_info_source}"
  )

  set("${output_variable}" "${_openfms_build_info_source}" PARENT_SCOPE)
endfunction()

function(openfms_configure_build_info_source source)
  # Keep the build date fresh for each CMake configure without baking it into
  # the generated file. This lets generate_build_info.sh remain responsible for
  # stable metadata while CMake injects the timestamp as a compile definition.
  string(TIMESTAMP _openfms_build_date "%Y-%m-%d %H:%M:%S %z")
  set_source_files_properties("${source}" PROPERTIES
    COMPILE_DEFINITIONS "OPENFMS_BUILD_DATE=\"${_openfms_build_date}\""
  )
endfunction()
