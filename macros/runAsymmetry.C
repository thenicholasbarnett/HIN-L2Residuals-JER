// Compiled:    ./bin/runAsymmetry <input.root> <output.root> [--monte-carlo|-mc|--zero-bias|-zb|--hard-probes|-hp] [maxEvents]
// Interpreted: root -l -b -q 'macros/runAsymmetry.C("in.root","out.root")'
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

#include "RunAsymmetry.h"

#ifndef __CLING__
#include <iostream>
#include <cstdlib>
#include <exception>
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: runAsymmetry <input.root> <output.root>"
                     " [--mc|--zero-bias|--hard-probes] [maxEvents]\n";
        return 1;
    }
    Long64_t maxEvents = (argc > 4) ? std::atoll(argv[4]) : -1;
    try {
        runAsymmetry(argv[1], argv[2], argc > 3 ? argv[3] : "--hard-probes", maxEvents);
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
    return 0;
}
#endif
