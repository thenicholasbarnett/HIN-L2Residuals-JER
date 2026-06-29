cmake_minimum_required(VERSION 3.10)

# Must be before project(): CMSSW's cvmfs GCC has annobin baked into its specs
# but the plugin .so is absent, causing every compile to fail. Use system GCC instead.
if(DEFINED ENV{CMSSW_BASE})
    # System GCC: CMSSW's GCC has broken annobin plugin specs on cvmfs
    set(CMAKE_C_COMPILER   /usr/bin/gcc CACHE PATH "" FORCE)
    set(CMAKE_CXX_COMPILER /usr/bin/g++ CACHE PATH "" FORCE)
    # CMSSW ROOT: system ROOT on lxplus lacks libCling.so; CMSSW ROOT has it.
    # Find CMSSW's ROOT via root-config so the binary links against a consistent ROOT+cling.
    execute_process(COMMAND root-config --prefix
        OUTPUT_VARIABLE _cmssw_root_prefix
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    foreach(_subdir cmake lib/cmake/ROOT)
        if(EXISTS "${_cmssw_root_prefix}/${_subdir}/ROOTConfig.cmake")
            set(ROOT_DIR "${_cmssw_root_prefix}/${_subdir}" CACHE PATH "" FORCE)
            break()
        endif()
    endforeach()
endif()

project(L2Residuals)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

# find ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Tree Hist MathCore)
include(${ROOT_USE_FILE})

# Embed the found ROOT lib path as RPATH so cmsenv's LD_LIBRARY_PATH cannot
# swap in a different ROOT version at runtime
set(CMAKE_BUILD_RPATH "${ROOT_LIBRARY_DIR}")

# nlohmann-json (required by JSON_handler.h for golden JSON support)
find_package(nlohmann_json QUIET)
if(NOT nlohmann_json_FOUND)
    # Fall back to pkg-config on systems where the CMake package isn't registered
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(NLOHMANN_JSON REQUIRED nlohmann_json)
    set(NLOHMANN_JSON_INCLUDE_DIR ${NLOHMANN_JSON_INCLUDE_DIRS})
else()
    set(NLOHMANN_JSON_INCLUDE_DIR "")  # header-only; imported target carries it
endif()

# shared library from src/
file(GLOB SOURCES src/*.cxx)
if(SOURCES)
    add_library(l2residuals SHARED ${SOURCES})
    target_include_directories(l2residuals PUBLIC include cfg ${NLOHMANN_JSON_INCLUDE_DIR})
    if(nlohmann_json_FOUND)
        target_link_libraries(l2residuals ${ROOT_LIBRARIES} nlohmann_json::nlohmann_json)
    else()
        target_link_libraries(l2residuals ${ROOT_LIBRARIES})
    endif()
else()
    # header-only placeholder until src/ is populated
    add_library(l2residuals INTERFACE)
    target_include_directories(l2residuals INTERFACE include cfg)
    target_link_libraries(l2residuals INTERFACE ${ROOT_LIBRARIES})
endif()

# build each entry-point in macros/ as an executable, output to bin/
file(GLOB MACROS macros/*.C macros/*.cxx)
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
file(GLOB TESTS tests/*.cxx)
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
endforeach()

add_test(NAME TestBuild
    COMMAND bash ${PROJECT_SOURCE_DIR}/tests/TestBuild.sh
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
