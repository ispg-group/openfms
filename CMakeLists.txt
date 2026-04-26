cmake_minimum_required(VERSION 3.18)

project(OpenFMS LANGUAGES Fortran)

option(OPENFMS_ENABLE_TERACHEM "Build the TeraChem backend" OFF)
option(OPENFMS_ENABLE_QUANTICS "Build the Quantics backend" OFF)
option(OPENFMS_ENABLE_NATIVE_OPTIMIZATION "Optimize release builds for the current CPU" OFF)

# Make CMake build default to "Debug" and provide valid options
if(NOT CMAKE_CONFIGURATION_TYPES)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
  endif()
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    Debug Release RelWithDebInfo MinSizeRel
  )
endif()

# TeraChem and Quantics not (yet) supported via CMake
if(OPENFMS_ENABLE_TERACHEM)
  message(FATAL_ERROR "The CMake TeraChem backend is not wired yet. Use the existing Makefile build for now.")
endif()

if(OPENFMS_ENABLE_QUANTICS)
  message(FATAL_ERROR "The CMake Quantics backend is not wired yet. Use the existing Makefile build for now.")
endif()

include(cmake/OpenFMSSources.cmake)

# Match the Makefile build by embedding the current git revision in the banner.
find_package(Git QUIET)
if(Git_FOUND)
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse HEAD
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE OPENFMS_GIT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" status --porcelain
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE OPENFMS_GIT_STATUS
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if(OPENFMS_GIT_STATUS)
    string(APPEND OPENFMS_GIT_VERSION "+")
  endif()
endif()

if(NOT OPENFMS_GIT_VERSION)
  set(OPENFMS_GIT_VERSION "unknown")
endif()

add_executable(openfms ${OPENFMS_SOURCES})

# Keep generated Fortran module files in the build tree.
set_target_properties(openfms PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
)

target_include_directories(openfms PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/modules"
  "${CMAKE_CURRENT_BINARY_DIR}/modules"
)

target_compile_definitions(openfms PRIVATE
  GIT_VER="${OPENFMS_GIT_VERSION}"
)

# Compiler flags are grouped by compiler and build type, then attached below.
set(OPENFMS_GNU_WARNING_FLAGS
  -Wno-unused-dummy-argument
  -Wno-unused-function
  -Wno-maybe-uninitialized
)

set(OPENFMS_GNU_DEBUG_FLAGS
  -Og
  -g
  -fcheck=all
  -ffpe-trap=zero,overflow
  -fimplicit-none
)
set(OPENFMS_GNU_RELEASE_FLAGS -O3)
set(OPENFMS_GNU_RELWITHDEBINFO_FLAGS -O2 -g)
set(OPENFMS_GNU_MINSIZEREL_FLAGS -Os)

set(OPENFMS_INTEL_DEBUG_FLAGS
  -O0
  -g
  -traceback
  -check
  bounds
)
set(OPENFMS_INTEL_RELEASE_FLAGS -O3)
set(OPENFMS_INTEL_RELWITHDEBINFO_FLAGS -O2 -g)
set(OPENFMS_INTEL_MINSIZEREL_FLAGS -Os)

# These helpers hide the generator expressions needed for target-local flags.
function(openfms_add_fortran_options compiler_ids)
  target_compile_options(openfms PRIVATE
    "$<$<COMPILE_LANG_AND_ID:Fortran,${compiler_ids}>:${ARGN}>"
  )
endfunction()

function(openfms_add_fortran_config_options compiler_ids config)
  target_compile_options(openfms PRIVATE
    "$<$<AND:$<COMPILE_LANG_AND_ID:Fortran,${compiler_ids}>,$<CONFIG:${config}>>:${ARGN}>"
  )
endfunction()

openfms_add_fortran_options("GNU" ${OPENFMS_GNU_WARNING_FLAGS})
openfms_add_fortran_config_options("GNU" Debug ${OPENFMS_GNU_DEBUG_FLAGS})
openfms_add_fortran_config_options("GNU" Release ${OPENFMS_GNU_RELEASE_FLAGS})
openfms_add_fortran_config_options("GNU" RelWithDebInfo ${OPENFMS_GNU_RELWITHDEBINFO_FLAGS})
openfms_add_fortran_config_options("GNU" MinSizeRel ${OPENFMS_GNU_MINSIZEREL_FLAGS})

openfms_add_fortran_config_options("Intel,IntelLLVM" Debug ${OPENFMS_INTEL_DEBUG_FLAGS})
openfms_add_fortran_config_options("Intel,IntelLLVM" Release ${OPENFMS_INTEL_RELEASE_FLAGS})
openfms_add_fortran_config_options("Intel,IntelLLVM" RelWithDebInfo ${OPENFMS_INTEL_RELWITHDEBINFO_FLAGS})
openfms_add_fortran_config_options("Intel,IntelLLVM" MinSizeRel ${OPENFMS_INTEL_MINSIZEREL_FLAGS})

if(OPENFMS_ENABLE_NATIVE_OPTIMIZATION)
  target_compile_options(openfms PRIVATE
    $<$<AND:$<COMPILE_LANG_AND_ID:Fortran,GNU>,$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>>:-march=native>
    $<$<AND:$<COMPILE_LANG_AND_ID:Fortran,Intel,IntelLLVM>,$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>>:-xHost>
  )
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # EISPACK is old fixed-form Fortran and needs GNU's legacy mode.
  set_source_files_properties(src/modules/EispackModule.f PROPERTIES
    COMPILE_OPTIONS "-std=legacy"
  )
endif()

enable_testing()

# Register each integration-test directory as its own CTest entry.
file(GLOB OPENFMS_TEST_DIRS
  LIST_DIRECTORIES true
  CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/tests/*"
)

foreach(OPENFMS_TEST_DIR IN LISTS OPENFMS_TEST_DIRS)
  if(IS_DIRECTORY "${OPENFMS_TEST_DIR}")
    get_filename_component(OPENFMS_TEST_NAME "${OPENFMS_TEST_DIR}" NAME)
    add_test(
      NAME "${OPENFMS_TEST_NAME}"
      COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/tests/test.sh" "$<TARGET_FILE:openfms>" "${OPENFMS_TEST_NAME}" test
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    set_tests_properties("${OPENFMS_TEST_NAME}" PROPERTIES
      LABELS "integration"
    )
  endif()
endforeach()
