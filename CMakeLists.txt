cmake_minimum_required(VERSION 3.10)

execute_process(COMMAND hostname
    OUTPUT_VARIABLE L2RESIDUALS_HOSTNAME
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
if(L2RESIDUALS_HOSTNAME MATCHES "lxplus" AND NOT DEFINED ENV{CMSSW_BASE})
    message(FATAL_ERROR
        "This looks like lxplus (${L2RESIDUALS_HOSTNAME}), but cmsenv is not active.\n"
        "Builds intended for HTCondor must use CMSSW GCC + CMSSW ROOT.\n"
        "Run something like:\n"
        "  source /cvmfs/cms.cern.ch/cmsset_default.sh\n"
        "  cd <CMSSW_RELEASE>/src\n"
        "  cmsenv\n"
        "  cd /path/to/L2Residuals-2024ppref\n"
        "  cmake -S . -B build\n"
        "  cmake --build build")
endif()

# Must be before project(): CMSSW GCC has a broken annobin plugin and
# system ROOT on lxplus lacks libCling.so — both fixed with CMSSW GCC + CMSSW ROOT.
if(DEFINED ENV{CMSSW_BASE})
    # CMSSW GCC has broken annobin plugin specs (the .so is missing from cvmfs).
    # Fix: dump its specs, strip annobin entries, compile with -specs=<cleaned>.
    execute_process(COMMAND g++ -dumpspecs
        OUTPUT_VARIABLE _specs OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    string(REGEX REPLACE "-fplugin=annobin" "" _specs "${_specs}")
    string(REGEX REPLACE "-fplugin-arg-annobin-[^ \n]*" "" _specs "${_specs}")
    set(_specs_file "${CMAKE_CURRENT_LIST_DIR}/.noannobin.specs")
    file(WRITE "${_specs_file}" "${_specs}")
    set(CMAKE_CXX_FLAGS_INIT "-specs=${_specs_file}")
    set(CMAKE_C_FLAGS_INIT   "-specs=${_specs_file}")

    # The stripped specs may not propagate CMSSW GCC's own lib dir to ld, causing
    # "undefined reference to GLIBCXX_3.4.31" when linking against CMSSW ROOT libs.
    # Explicitly add CMSSW GCC's libstdc++ dir to ld's search path and RPATH.
    execute_process(COMMAND g++ -print-file-name=libstdc++.so.6
        OUTPUT_VARIABLE _libstdcxx OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(_libstdcxx AND NOT "${_libstdcxx}" STREQUAL "libstdc++.so.6")
        get_filename_component(_gxx_lib_dir "${_libstdcxx}" DIRECTORY)
        set(CMAKE_EXE_LINKER_FLAGS
            "${CMAKE_EXE_LINKER_FLAGS} -L${_gxx_lib_dir} -Wl,-rpath,${_gxx_lib_dir}")
        set(CMAKE_SHARED_LINKER_FLAGS
            "${CMAKE_SHARED_LINKER_FLAGS} -L${_gxx_lib_dir} -Wl,-rpath,${_gxx_lib_dir}")
    endif()

    # CMSSW ROOT: system ROOT on lxplus lacks libCling.so; CMSSW ROOT has it.
    # Both compiler and ROOT must come from the same GCC 13 ABI.
    execute_process(COMMAND root-config --prefix
        OUTPUT_VARIABLE _r OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    foreach(_subdir cmake lib/cmake/ROOT)
        if(EXISTS "${_r}/${_subdir}/ROOTConfig.cmake")
            set(ROOT_DIR "${_r}/${_subdir}" CACHE PATH "" FORCE)
            break()
        endif()
    endforeach()

    # EL9's ld.bfd defaults to DT_RPATH already, but make it explicit so future
    # linker upgrades can't silently switch to DT_RUNPATH (which LD_LIBRARY_PATH
    # can override, allowing the wrong ROOT to load on Condor worker nodes).
    set(CMAKE_EXE_LINKER_FLAGS
        "${CMAKE_EXE_LINKER_FLAGS} -Wl,--disable-new-dtags")
    set(CMAKE_SHARED_LINKER_FLAGS
        "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--disable-new-dtags")
endif()

project(L2Residuals)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Arch ROOT 6.40.02 declares a Vdt CMake dependency that isn't installed.
# Prepend our cmake/ stub so it takes precedence over ROOT's FindVdt.cmake.
if(NOT DEFINED ENV{CMSSW_BASE})
    list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
endif()

# find ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Tree Hist MathCore)
include(${ROOT_USE_FILE})

# Embed the found ROOT lib path as RPATH so cmsenv's LD_LIBRARY_PATH cannot
# swap in a different ROOT version at runtime
set(CMAKE_BUILD_RPATH "${ROOT_LIBRARY_DIR}")

# nlohmann-json (required by JSON_handler.h for golden JSON support)
# Vendored single header takes priority — no network dependency
set(NLOHMANN_JSON_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external")

# shared library from src/ only. External vendored .cc files that need
# compiling are covered by thin .cxx wrappers in src/ (e.g. src/CommandLine.cxx,
# src/JetResolutionObject.cxx) so that both CMake and SCRAM pick them up via
# the same src/*.cxx glob without symlink or path-resolution ambiguity.
# Header-only vendored deps (toml++, nlohmann-json) need nothing here.
#
# JetResolutionObject.h self-defines STANDALONE (see external/jetmet_jer/) so
# every TU that includes it gets the TFormula-based standalone path rather than
# CMSSW's reco::FormulaEvaluator path. The PUBLIC flag below propagates the
# definition to tests and any downstream that links against l2residuals, keeping
# the class layout consistent across all translation units.
file(GLOB SOURCES CONFIGURE_DEPENDS src/*.cxx)
if(SOURCES)
    add_library(l2residuals SHARED ${SOURCES})
    target_include_directories(l2residuals PUBLIC
        include cfg
        ${NLOHMANN_JSON_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/external/jetmet
        ${CMAKE_CURRENT_SOURCE_DIR}/external/jetmet_jer)
    target_compile_definitions(l2residuals PUBLIC L2RESIDUALS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    target_compile_definitions(l2residuals PUBLIC STANDALONE)
    target_link_libraries(l2residuals ${ROOT_LIBRARIES})
else()
    # header-only placeholder until src/ is populated
    add_library(l2residuals INTERFACE)
    target_include_directories(l2residuals INTERFACE include cfg)
    target_compile_definitions(l2residuals INTERFACE L2RESIDUALS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    target_link_libraries(l2residuals INTERFACE ${ROOT_LIBRARIES})
endif()

# build each entry-point in macros/ as a standalone executable
file(GLOB MACROS CONFIGURE_DEPENDS macros/*.C macros/*.cxx)
foreach(macro ${MACROS})
    get_filename_component(name ${macro} NAME_WE)
    add_executable(${name} ${macro})
    if(SOURCES)
        target_link_libraries(${name} l2residuals)
    else()
        target_include_directories(${name} PRIVATE include cfg)
        target_link_libraries(${name} ${ROOT_LIBRARIES})
    endif()
endforeach()

# tests
enable_testing()
file(GLOB TESTS CONFIGURE_DEPENDS tests/*.cxx)
foreach(test ${TESTS})
    get_filename_component(name ${test} NAME_WE)
    add_executable(${name} ${test})
    if(SOURCES)
        target_link_libraries(${name} l2residuals)
    else()
        target_include_directories(${name} PRIVATE include cfg)
        target_link_libraries(${name} ${ROOT_LIBRARIES})
    endif()
    add_test(NAME ${name} COMMAND ${name})
    # Keep test CWD at the repo root so any test-created output paths stay
    # predictable. CONFIG no longer has an implicit fallback (see
    # AnalysisConfig.cxx), so tests must set L2RESIDUALS_CONFIG explicitly
    # rather than relying on one.
    set_tests_properties(${name} PROPERTIES
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        ENVIRONMENT "L2RESIDUALS_CONFIG=${PROJECT_SOURCE_DIR}/cfg/2024ppRef.toml")
endforeach()

add_test(NAME TestBuild
    COMMAND bash ${PROJECT_SOURCE_DIR}/tests/TestBuild.sh
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
