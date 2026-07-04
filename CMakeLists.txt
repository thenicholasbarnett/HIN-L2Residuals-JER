cmake_minimum_required(VERSION 3.10)

execute_process(COMMAND hostname
    OUTPUT_VARIABLE L2RESIDUALS_HOSTNAME
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
if(L2RESIDUALS_HOSTNAME MATCHES "lxplus" AND NOT DEFINED ENV{CMSSW_BASE})
    message(FATAL_ERROR
        "This looks like lxplus (${L2RESIDUALS_HOSTNAME}), but cmsenv is not active.\n"
        "HTCondor builds need CMSSW GCC + CMSSW ROOT -- run cmsenv first.")
endif()

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/CMSSWCompat.cmake)

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

# Shared library built from src/*.cxx only -- vendored .cc files needing
# compilation ride in via thin .cxx wrappers (e.g. src/CommandLine.cxx) so
# CMake and SCRAM share one glob. Header-only vendored deps need nothing here.
#
# STANDALONE is PUBLIC so every TU, including tests and downstream users,
# sees JetResolutionObject.h's TFormula-based path consistently.
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
    # WORKING_DIRECTORY keeps test-created paths predictable; L2RESIDUALS_CONFIG
    # is set explicitly since Config() has no implicit fallback (AnalysisConfig.cxx).
    set_tests_properties(${name} PROPERTIES
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        ENVIRONMENT "L2RESIDUALS_CONFIG=${PROJECT_SOURCE_DIR}/cfg/2024ppRef.toml")
endforeach()

add_test(NAME TestBuild
    COMMAND bash ${PROJECT_SOURCE_DIR}/tests/TestBuild.sh
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
