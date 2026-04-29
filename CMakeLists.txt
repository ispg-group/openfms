cmake_minimum_required(VERSION 3.18)

project(OpenFMS LANGUAGES Fortran)

# Project options

set(OPENFMS_VERSION "v0.1-alpha")
option(OPENFMS_ENABLE_NATIVE_OPTIMIZATION "Optimize release builds for the current CPU" OFF)

if(NOT CMAKE_CONFIGURATION_TYPES)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
  endif()
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    Debug Release
  )
endif()

# CMake currently builds the standalone analytical-potential backend only.
# The Makefile remains the supported path for external interfaces.

set(OPENFMS_BACKEND "zero" CACHE STRING "OpenFMS backend")
set_property(CACHE OPENFMS_BACKEND PROPERTY STRINGS zero)
if(NOT OPENFMS_BACKEND STREQUAL "zero")
  message(FATAL_ERROR
    "CMake currently supports only OPENFMS_BACKEND=zero. "
    "Use the Makefile build for TeraChem or Quantics."
  )
endif()

# Source and compiler configuration

include(cmake/OpenFMSSources.cmake)
include(cmake/OpenFMSCompilerOptions.cmake)
include(cmake/OpenFMSBuildInfo.cmake)
include(cmake/OpenFMSTargets.cmake)

# Targets

openfms_add_build_info(OPENFMS_BUILD_INFO_SOURCE)

openfms_add_core_library(openfms_core)
openfms_add_driver(openfms openfms_core "${OPENFMS_BUILD_INFO_SOURCE}")

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # EISPACK is old fixed-form Fortran and needs GNU's legacy mode.
  set_source_files_properties(src/modules/EispackModule.f PROPERTIES
    COMPILE_OPTIONS "-std=legacy"
  )
endif()

# Tests

include(CTest)
if(BUILD_TESTING)
  include(cmake/OpenFMSUnitTests.cmake)
  openfms_add_unit_tests(openfms_core)
  include(cmake/OpenFMSIntegrationTests.cmake)
  openfms_add_integration_tests(openfms)
endif()
