# Stub for Arch Linux where ROOT 6.40.02 declares a Vdt CMake dependency
# that isn't actually required for our build. ROOT on Arch was compiled
# without vdt; this stub satisfies find_dependency(Vdt) from ROOTConfig.cmake.
set(VDT_FOUND TRUE)
set(Vdt_FOUND TRUE)   # find_dependency checks ${dep}_FOUND (mixed-case)
set(VDT_INCLUDE_DIRS "")
set(VDT_LIBRARIES "")
if(NOT TARGET VDT::VDT)
    add_library(VDT::VDT INTERFACE IMPORTED)
endif()
