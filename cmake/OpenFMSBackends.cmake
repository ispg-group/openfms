include_guard(GLOBAL)

set(_OPENFMS_DEFAULT_BACKEND zero)

if(NOT DEFINED CACHE{OPENFMS_BACKEND})
  set(_OPENFMS_LEGACY_BACKENDS)
  if(DEFINED OPENFMS_ENABLE_ZERO AND OPENFMS_ENABLE_ZERO)
    list(APPEND _OPENFMS_LEGACY_BACKENDS zero)
  endif()
  if(DEFINED OPENFMS_ENABLE_QUANTICS AND OPENFMS_ENABLE_QUANTICS)
    list(APPEND _OPENFMS_LEGACY_BACKENDS quantics)
  endif()
  if(DEFINED OPENFMS_ENABLE_TERACHEM AND OPENFMS_ENABLE_TERACHEM)
    list(APPEND _OPENFMS_LEGACY_BACKENDS terachem)
  endif()

  list(LENGTH _OPENFMS_LEGACY_BACKENDS _OPENFMS_LEGACY_BACKEND_COUNT)
  if(_OPENFMS_LEGACY_BACKEND_COUNT EQUAL 1)
    list(GET _OPENFMS_LEGACY_BACKENDS 0 _OPENFMS_DEFAULT_BACKEND)
  elseif(_OPENFMS_LEGACY_BACKEND_COUNT GREATER 1)
    message(FATAL_ERROR
      "Enable exactly one OpenFMS backend. "
      "Use OPENFMS_BACKEND=zero, terachem, or quantics."
    )
  endif()
endif()

set(OPENFMS_BACKEND "${_OPENFMS_DEFAULT_BACKEND}" CACHE STRING "OpenFMS backend")
set_property(CACHE OPENFMS_BACKEND PROPERTY STRINGS
  zero
  terachem
  quantics
)

function(openfms_select_backend output_variable)
  set(_openfms_supported_backends zero terachem quantics)

  if(NOT OPENFMS_BACKEND IN_LIST _openfms_supported_backends)
    message(FATAL_ERROR
      "Unsupported OpenFMS backend: ${OPENFMS_BACKEND}. "
      "Use OPENFMS_BACKEND=zero, terachem, or quantics."
    )
  endif()

  if(OPENFMS_BACKEND STREQUAL "terachem")
    message(FATAL_ERROR "The CMake TeraChem backend is not wired yet. Use the existing Makefile build for now.")
  elseif(OPENFMS_BACKEND STREQUAL "quantics")
    message(FATAL_ERROR "The CMake Quantics backend is not wired yet. Use the existing Makefile build for now.")
  endif()

  set("${output_variable}" "${OPENFMS_BACKEND}" PARENT_SCOPE)
endfunction()
