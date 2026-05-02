cmake_minimum_required(VERSION 3.18)

project(OpenFMS LANGUAGES Fortran)

set(OPENFMS_VERSION "v0.1-alpha")

option(OPENFMS_ENABLE_NATIVE_OPTIMIZATION "Optimize release builds for the current CPU" OFF)

# Set standard build type to debug if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

set(OPENFMS_BACKEND "zero" CACHE STRING "OpenFMS backend")

# CMake currently only supports the standalone analytical-potential
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

# Add targets
openfms_add_build_info(OPENFMS_BUILD_INFO_SOURCE)
openfms_add_core_library(openfms_core)
openfms_add_executable(openfms openfms_core "${OPENFMS_BUILD_INFO_SOURCE}")

# Compiler flag exception: EISPACK is old fixed-form Fortran
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set_source_files_properties(
    src/modules/EispackModule.f PROPERTIES COMPILE_OPTIONS "-std=legacy"
  )
endif()

# Add test suite
include(CTest)

include(cmake/OpenFMSUnitTests.cmake)
openfms_add_unit_tests(openfms_core)

include(cmake/OpenFMSIntegrationTests.cmake)
openfms_add_integration_tests(openfms)
