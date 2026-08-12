function(openfms_add_integration_tests openfms_target)
  # Loop through tests/ and add all directories with capitalized first letter
  # (note: capitalization is a convention in project for excluding tests)
  file(GLOB OPENFMS_TEST_DIRS
    LIST_DIRECTORIES true
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/*"
  )

  foreach(OPENFMS_TEST_DIR IN LISTS OPENFMS_TEST_DIRS)
    if(IS_DIRECTORY "${OPENFMS_TEST_DIR}")
      get_filename_component(OPENFMS_TEST_NAME "${OPENFMS_TEST_DIR}" NAME)
      # Only add test for capitalized first letter
      if(OPENFMS_TEST_NAME MATCHES "^[A-Z]") 
        add_test(
          NAME "${OPENFMS_TEST_NAME}"
          COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/tests/test.sh" "$<TARGET_FILE:${openfms_target}>" "${OPENFMS_TEST_NAME}" test
          WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        )
        set_tests_properties("${OPENFMS_TEST_NAME}" PROPERTIES
          LABELS "integration"
        )
      endif()
    endif()
  endforeach()
endfunction()
