function(openfms_add_unit_tests openfms_core_target)
  set(OPENFMS_UNIT_TEST_SOURCES
    unit_tests/main.F90
    unit_tests/test_bundle.F90
    unit_tests/test_particle.F90
    unit_tests/test_trajectory.F90
    unit_tests/test_overlap.F90
    unit_tests/testdrive.F90
    unit_tests/testutils.F90
  )

  add_executable(openfms_unit_tests ${OPENFMS_UNIT_TEST_SOURCES})
  target_link_libraries(openfms_unit_tests PRIVATE "${openfms_core_target}")

  # Keep test .mod files separate from the production module directory
  set_target_properties(openfms_unit_tests PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/unit_test_modules"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_include_directories(openfms_unit_tests PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/unit_tests"
    "${CMAKE_CURRENT_BINARY_DIR}/unit_test_modules"
  )
  openfms_configure_fortran_target(openfms_unit_tests)

  set(OPENFMS_UNIT_TEST_SUITES
    BundleModule
    ParticleModule
    TrajectoryModule
  )

  foreach(suite_name IN LISTS OPENFMS_UNIT_TEST_SUITES)
    set(_openfms_ctest_name "unit.${suite_name}")

    add_test(
      NAME "${_openfms_ctest_name}"
      COMMAND "$<TARGET_FILE:openfms_unit_tests>" "${suite_name}"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    set_tests_properties("${_openfms_ctest_name}" PROPERTIES
      LABELS "unit;${suite_name}"
    )
  endforeach()

endfunction()
