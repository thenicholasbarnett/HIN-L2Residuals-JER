# Included from CMakeLists.txt before project() -- fixes CMSSW GCC's broken
# annobin plugin, points CMake at CMSSW ROOT (system ROOT on lxplus lacks
# libCling.so). No-op unless CMSSW_BASE is set.
if(DEFINED ENV{CMSSW_BASE})
    execute_process(COMMAND g++ -dumpspecs
        OUTPUT_VARIABLE _specs OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    string(REGEX REPLACE "-fplugin=annobin" "" _specs "${_specs}")
    string(REGEX REPLACE "-fplugin-arg-annobin-[^ \n]*" "" _specs "${_specs}")
    set(_specs_file "${CMAKE_CURRENT_LIST_DIR}/.noannobin.specs")
    file(WRITE "${_specs_file}" "${_specs}")
    set(CMAKE_CXX_FLAGS_INIT "-specs=${_specs_file}")
    set(CMAKE_C_FLAGS_INIT "-specs=${_specs_file}")

    # Stripping the specs can drop CMSSW GCC's own libstdc++ dir, breaking
    # linking ("undefined reference to GLIBCXX_3.4.31") -- add it back explicitly.
    execute_process(COMMAND g++ -print-file-name=libstdc++.so.6
        OUTPUT_VARIABLE _libstdcxx OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(_libstdcxx AND NOT "${_libstdcxx}" STREQUAL "libstdc++.so.6")
        get_filename_component(_gxx_lib_dir "${_libstdcxx}" DIRECTORY)
        set(CMAKE_EXE_LINKER_FLAGS
            "${CMAKE_EXE_LINKER_FLAGS} -L${_gxx_lib_dir} -Wl,-rpath,${_gxx_lib_dir}")
        set(CMAKE_SHARED_LINKER_FLAGS
            "${CMAKE_SHARED_LINKER_FLAGS} -L${_gxx_lib_dir} -Wl,-rpath,${_gxx_lib_dir}")
    endif()

    # System ROOT lacks libCling.so; CMSSW ROOT has it and matches the GCC ABI.
    execute_process(COMMAND root-config --prefix
        OUTPUT_VARIABLE _r OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    foreach(_subdir cmake lib/cmake/ROOT)
        if(EXISTS "${_r}/${_subdir}/ROOTConfig.cmake")
            set(ROOT_DIR "${_r}/${_subdir}" CACHE PATH "" FORCE)
            break()
        endif()
    endforeach()

    # Make DT_RPATH explicit (EL9's ld.bfd default, but not guaranteed) so a
    # future linker can't switch to DT_RUNPATH, which LD_LIBRARY_PATH could
    # override and load the wrong ROOT on a Condor worker.
    set(CMAKE_EXE_LINKER_FLAGS
        "${CMAKE_EXE_LINKER_FLAGS} -Wl,--disable-new-dtags")
    set(CMAKE_SHARED_LINKER_FLAGS
        "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--disable-new-dtags")
endif()
