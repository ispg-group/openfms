include_guard(GLOBAL)

function(openfms_add_core_library target)
  add_library("${target}" OBJECT
    ${OPENFMS_CORE_SOURCES}
  )

  set_target_properties("${target}" PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
  )

  target_include_directories("${target}" PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/modules"
    "${CMAKE_CURRENT_BINARY_DIR}/modules"
  )

  openfms_configure_fortran_target("${target}")
endfunction()

function(openfms_add_driver target core_target build_info_source)
  openfms_configure_build_info_source("${build_info_source}")

  add_executable("${target}"
    src/openfms.F90
    "${build_info_source}"
  )
  add_dependencies("${target}" openfms_build_info)
  target_link_libraries("${target}" PRIVATE "${core_target}")

  set_target_properties("${target}" PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )

  openfms_configure_fortran_target("${target}")
endfunction()
