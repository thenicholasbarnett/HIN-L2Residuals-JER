cmake_minimum_required(VERSION 3.10)
project(L2Residuals)

set(CMAKE_CXX_STANDARD 17)

# find ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Tree Hist MathCore)
include(${ROOT_USE_FILE})

# shared library from src/
file(GLOB SOURCES src/*.cxx)
if(SOURCES)
    add_library(l2residuals SHARED ${SOURCES})
    target_include_directories(l2residuals PUBLIC include configs)
    target_link_libraries(l2residuals ${ROOT_LIBRARIES})
else()
    # header-only placeholder until src/ is populated
    add_library(l2residuals INTERFACE)
    target_include_directories(l2residuals INTERFACE include configs)
    target_link_libraries(l2residuals INTERFACE ${ROOT_LIBRARIES})
endif()

# build each macro as an executable
file(GLOB MACROS macros/*.C)
foreach(macro ${MACROS})
    get_filename_component(name ${macro} NAME_WE)
    add_executable(${name} ${macro})
    if(SOURCES)
        target_link_libraries(${name} l2residuals)
    else()
        target_include_directories(${name} PRIVATE include configs)
        target_link_libraries(${name} ${ROOT_LIBRARIES})
    endif()
endforeach()
