include_guard(GLOBAL)

function(openfms_add_unit_tests openfms_core_target)
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

  openfms_discover_unit_tests(OPENFMS_UNIT_CTEST_ENTRIES)
  foreach(OPENFMS_UNIT_CTEST_ENTRY IN LISTS OPENFMS_UNIT_CTEST_ENTRIES)
    openfms_add_discovered_unit_test("${OPENFMS_UNIT_CTEST_ENTRY}")
  endforeach()
endfunction()

function(openfms_add_discovered_unit_test ctest_entry)
  string(REGEX REPLACE "^([^|]+)\\|(.+)$" "\\1" _openfms_unit_suite_name "${ctest_entry}")
  string(REGEX REPLACE "^([^|]+)\\|(.+)$" "\\2" _openfms_unit_test_name "${ctest_entry}")
  set(_openfms_ctest_name "unit.${_openfms_unit_suite_name}.${_openfms_unit_test_name}")

  add_test(
    NAME "${_openfms_ctest_name}"
    COMMAND "$<TARGET_FILE:openfms_unit_tests>" "${_openfms_unit_suite_name}" "${_openfms_unit_test_name}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  )
  set_tests_properties("${_openfms_ctest_name}" PROPERTIES
    LABELS "unit;${_openfms_unit_suite_name}"
  )
endfunction()

function(openfms_discover_unit_suites output_variable)
  # Unit suites are discovered from:
  #   unit_tests/main.F90: new_testsuite("SuiteName", collect_suite)
  # The resulting entries are stored as collect_suite|SuiteName.
  file(READ "${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/main.F90" _openfms_unit_test_main)
  string(REGEX MATCHALL
    "new_testsuite[ \t\r\n]*\\([ \t\r\n]*\"[^\"]+\"[ \t\r\n]*,[ \t\r\n]*[A-Za-z_][A-Za-z0-9_]*"
    _openfms_unit_suite_matches
    "${_openfms_unit_test_main}"
  )
  if(NOT _openfms_unit_suite_matches)
    message(FATAL_ERROR "No unit test suites were discovered in unit_tests/main.F90")
  endif()

  set(_openfms_unit_suite_entries)
  foreach(_openfms_unit_suite_match IN LISTS _openfms_unit_suite_matches)
    string(REGEX REPLACE
      ".*\"([^\"]+)\"[ \t\r\n]*,[ \t\r\n]*([A-Za-z_][A-Za-z0-9_]*).*"
      "\\2|\\1"
      _openfms_unit_suite_entry
      "${_openfms_unit_suite_match}"
    )
    list(APPEND _openfms_unit_suite_entries "${_openfms_unit_suite_entry}")
  endforeach()

  set("${output_variable}" ${_openfms_unit_suite_entries} PARENT_SCOPE)
endfunction()

function(openfms_lookup_unit_suite suite_entries collect_name output_variable)
  foreach(_openfms_unit_suite_entry IN LISTS suite_entries)
    string(REGEX REPLACE "^([^|]+)\\|(.+)$" "\\1" _openfms_collect_name "${_openfms_unit_suite_entry}")
    string(REGEX REPLACE "^([^|]+)\\|(.+)$" "\\2" _openfms_suite_name "${_openfms_unit_suite_entry}")
    if("${_openfms_collect_name}" STREQUAL "${collect_name}")
      set("${output_variable}" "${_openfms_suite_name}" PARENT_SCOPE)
      return()
    endif()
  endforeach()

  set("${output_variable}" "" PARENT_SCOPE)
endfunction()

function(openfms_discover_unit_source_tests test_source suite_entries output_variable)
  file(READ "${test_source}" _openfms_unit_test_source_contents)
  string(REGEX MATCH
    "subroutine[ \t]+(collect_[A-Za-z0-9_]+)[ \t\r\n]*\\("
    _openfms_unit_collect_match
    "${_openfms_unit_test_source_contents}"
  )
  if(NOT _openfms_unit_collect_match)
    message(FATAL_ERROR "No collect_*_suite subroutine found in ${test_source}")
  endif()
  set(_openfms_unit_collect_name "${CMAKE_MATCH_1}")

  openfms_lookup_unit_suite("${suite_entries}" "${_openfms_unit_collect_name}" _openfms_unit_suite_name)
  if(NOT _openfms_unit_suite_name)
    message(FATAL_ERROR
      "Unit test collect routine ${_openfms_unit_collect_name} from ${test_source} "
      "is not registered in unit_tests/main.F90"
    )
  endif()

  string(REGEX MATCHALL
    "new_unittest[ \t\r\n]*\\([ \t\r\n]*\"[^\"]+\""
    _openfms_unit_test_matches
    "${_openfms_unit_test_source_contents}"
  )
  if(NOT _openfms_unit_test_matches)
    message(FATAL_ERROR "No new_unittest entries found in ${test_source}")
  endif()

  set(_openfms_unit_ctest_entries)
  foreach(_openfms_unit_test_match IN LISTS _openfms_unit_test_matches)
    string(REGEX REPLACE
      ".*\"([^\"]+)\".*"
      "\\1"
      _openfms_unit_test_name
      "${_openfms_unit_test_match}"
    )
    list(APPEND _openfms_unit_ctest_entries "${_openfms_unit_suite_name}|${_openfms_unit_test_name}")
  endforeach()

  set("${output_variable}" ${_openfms_unit_ctest_entries} PARENT_SCOPE)
endfunction()

function(openfms_discover_unit_tests output_variable)
  file(GLOB OPENFMS_UNIT_TEST_SUITE_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/test_*.F90"
  )

  openfms_discover_unit_suites(_openfms_unit_suite_entries)
  set(OPENFMS_UNIT_CTEST_NAMES)
  set(OPENFMS_UNIT_CTEST_ENTRIES)
  foreach(OPENFMS_UNIT_TEST_SOURCE IN LISTS OPENFMS_UNIT_TEST_SUITE_SOURCES)
    openfms_discover_unit_source_tests(
      "${OPENFMS_UNIT_TEST_SOURCE}"
      "${_openfms_unit_suite_entries}"
      _openfms_unit_source_ctest_entries
    )
    foreach(OPENFMS_UNIT_CTEST_ENTRY IN LISTS _openfms_unit_source_ctest_entries)
      set(OPENFMS_UNIT_CTEST_NAME "unit.${OPENFMS_UNIT_CTEST_ENTRY}")
      if(OPENFMS_UNIT_CTEST_NAME IN_LIST OPENFMS_UNIT_CTEST_NAMES)
        message(FATAL_ERROR "Duplicate unit test name discovered: ${OPENFMS_UNIT_CTEST_NAME}")
      endif()
      list(APPEND OPENFMS_UNIT_CTEST_NAMES "${OPENFMS_UNIT_CTEST_NAME}")
      list(APPEND OPENFMS_UNIT_CTEST_ENTRIES "${OPENFMS_UNIT_CTEST_ENTRY}")
    endforeach()
  endforeach()

  set("${output_variable}" ${OPENFMS_UNIT_CTEST_ENTRIES} PARENT_SCOPE)
endfunction()
