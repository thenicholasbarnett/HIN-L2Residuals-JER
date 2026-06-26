// Compiled:    ./bin/runAsymmetry <input.root> <output.root> [--mc|--zero-bias|--hard-probes]
// Interpreted: root -l -b -q 'macros/runAsymmetry.C("in.root","out.root")'
//              (build the library first: cmake --build build)
//              (run from the repo root so relative paths resolve correctly)

#ifdef __CLING__
R__ADD_INCLUDE_PATH(include)
R__ADD_INCLUDE_PATH(cfg)
R__LOAD_LIBRARY(lib/libl2residuals.so)
#endif

#include "RunAsymmetry.h"

#ifndef __CLING__
#include <iostream>
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: runAsymmetry <input.root> <output.root> [--mc|--zero-bias|--hard-probes]\n";
        return 1;
    }
    runAsymmetry(argv[1], argv[2], argc > 3 ? argv[3] : "--hard-probes");
    return 0;
}
#endif
