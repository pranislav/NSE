cmake_minimum_required(VERSION 3.13)

if(NOT DEFINED BUNDLE_EXECUTABLE OR NOT EXISTS "${BUNDLE_EXECUTABLE}")
  message(FATAL_ERROR "BUNDLE_EXECUTABLE must point to a built executable")
endif()

if(NOT DEFINED BUNDLE_DIR)
  message(FATAL_ERROR "BUNDLE_DIR must be set")
endif()

if(NOT DEFINED BUNDLE_SYSTEM_LIBRARY_REGEX OR BUNDLE_SYSTEM_LIBRARY_REGEX STREQUAL "")
  set(BUNDLE_SYSTEM_LIBRARY_REGEX
    "^(linux-vdso|ld-linux|ld64|libc|libm|libmvec|libpthread|librt|libdl|libresolv|libnsl|libutil|libanl)(-|\\.so|$)"
    )
endif()

set(bundle_lib_dir "${BUNDLE_DIR}/lib")

file(REMOVE_RECURSE "${BUNDLE_DIR}")
file(MAKE_DIRECTORY "${bundle_lib_dir}")
get_filename_component(executable_name "${BUNDLE_EXECUTABLE}" NAME)
set(bundled_binary_name ".${executable_name}.bin")
file(COPY "${BUNDLE_EXECUTABLE}" DESTINATION "${BUNDLE_DIR}")
file(RENAME "${BUNDLE_DIR}/${executable_name}" "${BUNDLE_DIR}/${bundled_binary_name}")

function(collect_ldd_dependencies binary_path output_variable)
  execute_process(
    COMMAND ldd "${binary_path}"
    RESULT_VARIABLE ldd_result
    OUTPUT_VARIABLE ldd_output
    ERROR_VARIABLE ldd_error
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE
    )

  if(NOT ldd_result EQUAL 0)
    message(FATAL_ERROR "ldd failed for ${binary_path}: ${ldd_error}")
  endif()

  set(found_dependencies "")
  string(REPLACE "\n" ";" ldd_lines "${ldd_output}")
  foreach(line IN LISTS ldd_lines)
    string(STRIP "${line}" line)

    if(line MATCHES "not found")
      message(FATAL_ERROR "Unresolved shared library dependency for ${binary_path}: ${line}")
    endif()

    set(candidate "")
    if(line MATCHES "^[^=]+=>[ \t]+([^ \t]+)[ \t]+\\(")
      set(candidate "${CMAKE_MATCH_1}")
    elseif(line MATCHES "^(/[^\t ]+)[ \t]+\\(")
      set(candidate "${CMAKE_MATCH_1}")
    endif()

    if(candidate STREQUAL "" OR NOT EXISTS "${candidate}")
      continue()
    endif()

    get_filename_component(candidate_name "${candidate}" NAME)
    if(candidate_name MATCHES "${BUNDLE_SYSTEM_LIBRARY_REGEX}")
      continue()
    endif()

    list(APPEND found_dependencies "${candidate}")
  endforeach()
  list(REMOVE_DUPLICATES found_dependencies)
  set("${output_variable}" "${found_dependencies}" PARENT_SCOPE)
endfunction()

function(copy_library_with_symlinks source_path destination_dir)
  set(current "${source_path}")
  set(link_names "")

  while(IS_SYMLINK "${current}")
    get_filename_component(current_name "${current}" NAME)
    list(APPEND link_names "${current_name}")

    file(READ_SYMLINK "${current}" next_link)
    if(NOT IS_ABSOLUTE "${next_link}")
      get_filename_component(current_dir "${current}" DIRECTORY)
      set(next_link "${current_dir}/${next_link}")
    endif()
    get_filename_component(current "${next_link}" ABSOLUTE)
  endwhile()

  if(NOT EXISTS "${current}")
    message(FATAL_ERROR "Broken shared library symlink chain from ${source_path}")
  endif()

  get_filename_component(real_name "${current}" NAME)
  file(COPY "${current}" DESTINATION "${destination_dir}")

  foreach(link_name IN LISTS link_names)
    if(NOT link_name STREQUAL real_name)
      file(REMOVE "${destination_dir}/${link_name}")
      execute_process(
        COMMAND "${CMAKE_COMMAND}" -E create_symlink "${real_name}" "${destination_dir}/${link_name}"
        RESULT_VARIABLE symlink_result
        )
      if(NOT symlink_result EQUAL 0)
        message(FATAL_ERROR "Failed to create symlink ${destination_dir}/${link_name} -> ${real_name}")
      endif()
    endif()
  endforeach()
endfunction()

collect_ldd_dependencies("${BUNDLE_EXECUTABLE}" dependencies)

set(flexiblas_backend_libraries "")
foreach(dependency IN LISTS dependencies)
  get_filename_component(dependency_name "${dependency}" NAME)
  if(dependency_name MATCHES "^libflexiblas(64)?\\.so")
    get_filename_component(flexiblas_lib_dir "${dependency}" DIRECTORY)
    if(dependency_name MATCHES "^libflexiblas64")
      set(flexiblas_module_dir "${flexiblas_lib_dir}/flexiblas64")
    else()
      set(flexiblas_module_dir "${flexiblas_lib_dir}/flexiblas")
    endif()

    if(EXISTS "${flexiblas_module_dir}")
      file(GLOB flexiblas_modules "${flexiblas_module_dir}/libflexiblas*.so")
      list(APPEND flexiblas_backend_libraries ${flexiblas_modules})
    endif()
  endif()
endforeach()
list(REMOVE_DUPLICATES flexiblas_backend_libraries)

foreach(flexiblas_backend IN LISTS flexiblas_backend_libraries)
  collect_ldd_dependencies("${flexiblas_backend}" flexiblas_backend_dependencies)
  list(APPEND dependencies ${flexiblas_backend_dependencies})
endforeach()
list(REMOVE_DUPLICATES dependencies)

foreach(dependency IN LISTS dependencies)
  copy_library_with_symlinks("${dependency}" "${bundle_lib_dir}")
endforeach()

if(flexiblas_backend_libraries)
  set(bundle_flexiblas_dir "${bundle_lib_dir}/flexiblas")
  file(MAKE_DIRECTORY "${bundle_flexiblas_dir}")

  foreach(flexiblas_backend IN LISTS flexiblas_backend_libraries)
    copy_library_with_symlinks("${flexiblas_backend}" "${bundle_flexiblas_dir}")
  endforeach()

  file(WRITE "${BUNDLE_DIR}/flexiblasrc"
    "default = openblas-openmp\n"
    "\n"
    "[openblas-openmp]\n"
    "library = libflexiblas_openblas-openmp.so\n"
    "\n"
    "[netlib]\n"
    "library = libflexiblas_netlib.so\n"
    )
endif()

file(WRITE "${BUNDLE_DIR}/${executable_name}"
  "#!/usr/bin/env sh\n"
  "set -eu\n"
  "self_dir=$(CDPATH= cd -- \"$(dirname -- \"$0\")\" && pwd)\n"
  "if [ -d \"$self_dir/lib/flexiblas\" ]; then\n"
  "  export FLEXIBLAS_CONFIG=\"$self_dir/flexiblasrc\"\n"
  "  if [ -n \"\${FLEXIBLAS_LIBRARY_PATH:-}\" ]; then\n"
  "    export FLEXIBLAS_LIBRARY_PATH=\"$self_dir/lib/flexiblas:\$FLEXIBLAS_LIBRARY_PATH\"\n"
  "  else\n"
  "    export FLEXIBLAS_LIBRARY_PATH=\"$self_dir/lib/flexiblas\"\n"
  "  fi\n"
  "fi\n"
  "exec \"$self_dir/${bundled_binary_name}\" \"$@\"\n"
  )
execute_process(
  COMMAND chmod +x "${BUNDLE_DIR}/${executable_name}"
  RESULT_VARIABLE chmod_result
  )
if(NOT chmod_result EQUAL 0)
  message(FATAL_ERROR "Failed to make ${BUNDLE_DIR}/${executable_name} executable")
endif()

list(LENGTH dependencies dependency_count)
message(STATUS "Created ${BUNDLE_DIR} with ${executable_name} and ${dependency_count} shared libraries")
