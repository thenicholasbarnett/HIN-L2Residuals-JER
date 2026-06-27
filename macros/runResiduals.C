// Compiled:    ./bin/runResiduals <data.root> <mc.root> <output.root>
// Interpreted: root -l -b -q 'macros/runResiduals.C("data.root","mc.root","out.root")'
//              (build the library first: cmake --build build)
//              (run from the repo root so relative paths resolve correctly)

#ifdef __CLING__
R__ADD_INCLUDE_PATH(include)
R__ADD_INCLUDE_PATH(cfg)
#if defined(__APPLE__)
R__LOAD_LIBRARY(lib/libl2residuals.dylib)
#else
R__LOAD_LIBRARY(lib/libl2residuals.so)
#endif
#endif

#include "ResidualsExtractor.h"

#ifndef __CLING__
#include <iostream>
int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: runResiduals <data.root> <mc.root> <output.root>\n";
        return 1;
    }
    runResiduals(argv[1], argv[2], argv[3]);
    return 0;
}
#endif
