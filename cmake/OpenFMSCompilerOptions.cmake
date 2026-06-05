# GNU flags
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

# Intel flags
set(OPENFMS_INTEL_DEBUG_FLAGS
  -O0
  -g
  -traceback
  -check
  bounds
)
set(OPENFMS_INTEL_RELEASE_FLAGS -O3)

function(openfms_collect_fortran_flags output_variable)
  # Function collects the set of flags for the given configuration (Debug, Release)
  # and the given compiler (GNU, Intel) into output_variable
  set(_openfms_fortran_flags ${CMAKE_Fortran_FLAGS})

  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    list(APPEND _openfms_fortran_flags ${OPENFMS_GNU_WARNING_FLAGS})
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      list(APPEND _openfms_fortran_flags ${OPENFMS_GNU_DEBUG_FLAGS})
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
      list(APPEND _openfms_fortran_flags ${OPENFMS_GNU_RELEASE_FLAGS})
      if(OPENFMS_ENABLE_NATIVE_OPTIMIZATION)
        list(APPEND _openfms_fortran_flags -march=native)
      endif()
    endif()
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      list(APPEND _openfms_fortran_flags ${OPENFMS_INTEL_DEBUG_FLAGS})
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
      list(APPEND _openfms_fortran_flags ${OPENFMS_INTEL_RELEASE_FLAGS})
      if(OPENFMS_ENABLE_NATIVE_OPTIMIZATION)
        list(APPEND _openfms_fortran_flags -xHost)
      endif()
    endif()
  endif()

  # Remove semi colons between list items before returning full string 
  # of flags
  string(REPLACE ";" " " _openfms_fortran_flags "${_openfms_fortran_flags}")
  set("${output_variable}" "${_openfms_fortran_flags}" PARENT_SCOPE)
endfunction()

function(openfms_configure_fortran_target target)
  # Apply the Fortran options to a target that compiles project sources
  # The syntax $<$<... makes sure to add flags conditionally ( e.g., add debug flags if Debug config
  # is activated but otherwise not)
  target_compile_options("${target}" PRIVATE
    "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:${OPENFMS_GNU_WARNING_FLAGS}>"
    "$<$<AND:$<COMPILE_LANG_AND_ID:Fortran,GNU>,$<CONFIG:Debug>>:${OPENFMS_GNU_DEBUG_FLAGS}>"
    "$<$<AND:$<COMPILE_LANG_AND_ID:Fortran,GNU>,$<CONFIG:Release>>:${OPENFMS_GNU_RELEASE_FLAGS}>"
    "$<$<AND:$<COMPILE_LANG_AND_ID:Fortran,Intel,IntelLLVM>,$<CONFIG:Debug>>:${OPENFMS_INTEL_DEBUG_FLAGS}>"
    "$<$<AND:$<COMPILE_LANG_AND_ID:Fortran,Intel,IntelLLVM>,$<CONFIG:Release>>:${OPENFMS_INTEL_RELEASE_FLAGS}>"
  )

  if(OPENFMS_ENABLE_NATIVE_OPTIMIZATION)
    # Add native CPU optimization (note: makes executable non-transferrable to other machines)
    target_compile_options("${target}" PRIVATE
      $<$<AND:$<COMPILE_LANG_AND_ID:Fortran,GNU>,$<CONFIG:Release>>:-march=native>
      $<$<AND:$<COMPILE_LANG_AND_ID:Fortran,Intel,IntelLLVM>,$<CONFIG:Release>>:-xHost>
    )
  endif()
endfunction()
