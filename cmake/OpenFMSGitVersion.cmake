include_guard(GLOBAL)

function(openfms_get_git_version output_variable)
  set(OPENFMS_VERSION "unknown")

  find_package(Git QUIET)
  if(Git_FOUND)
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" rev-parse HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE OPENFMS_GIT_REVISION
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
    if(OPENFMS_GIT_REVISION)
      set(OPENFMS_VERSION "${OPENFMS_GIT_REVISION}")
      if(OPENFMS_GIT_STATUS)
        string(APPEND OPENFMS_VERSION "+")
      endif()
    endif()
  endif()

  set("${output_variable}" "${OPENFMS_VERSION}" PARENT_SCOPE)
endfunction()
