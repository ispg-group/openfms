include_guard(GLOBAL)

function(openfms_add_unit_tests openfms_core_target)

  # Function runs python script (discover_openfms_unit_tests.py)
  # to discover defined unit tests and generates a CMake file 
  # (OpenFMSDiscoveredUnitTests.cmake) that adds each unit test 
  # to the ctest suite as a separate entry

  find_package(Python3 REQUIRED COMPONENTS Interpreter)

  file(GLOB OPENFMS_UNIT_TEST_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/*.F90"
  )
  file(GLOB OPENFMS_UNIT_TEST_SUITE_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/test_*.F90"
  )
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/main.F90"
    ${OPENFMS_UNIT_TEST_SUITE_SOURCES}
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/discover_openfms_unit_tests.py"
  )

  add_executable(openfms_unit_tests ${OPENFMS_UNIT_TEST_SOURCES})
  target_link_libraries(openfms_unit_tests PRIVATE "${openfms_core_target}")
  set_target_properties(openfms_unit_tests PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/unit_test_modules"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_include_directories(openfms_unit_tests PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/unit_tests"
    "${CMAKE_CURRENT_BINARY_DIR}/unit_test_modules"
  )
  openfms_configure_fortran_target(openfms_unit_tests)

  set(_openfms_discovered_unit_tests
    "${CMAKE_CURRENT_BINARY_DIR}/OpenFMSDiscoveredUnitTests.cmake"
  )
  execute_process(
    COMMAND
      "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/discover_openfms_unit_tests.py"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    RESULT_VARIABLE _openfms_unit_test_discovery_result
  )
  if(NOT _openfms_unit_test_discovery_result EQUAL 0)
    message(FATAL_ERROR "OpenFMS unit test discovery failed")
  endif()
  include("${_openfms_discovered_unit_tests}")
  
endfunction()

function(openfms_add_discovered_unit_test suite_name test_name)

  # Function adds unit test to the ctest suite 
  # Called from the auto-generated CMake file OpenFMSDiscoveredUnitTests.cmake

  set(_openfms_ctest_name "unit.${suite_name}.${test_name}")

  add_test(
    NAME "${_openfms_ctest_name}"
    COMMAND "$<TARGET_FILE:openfms_unit_tests>" "${suite_name}" "${test_name}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  )
  set_tests_properties("${_openfms_ctest_name}" PROPERTIES
    LABELS "unit;${suite_name}"
  )

endfunction()
