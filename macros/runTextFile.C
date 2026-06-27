// Compiled:    ./bin/runTextFile <residuals.root> <output.txt> [method] [cone]
// Interpreted: root -l -b -q 'macros/runTextFile.C("residuals.root","out.txt")'
//              (build the library first: cmake --build build)
//              (run from the repo root so relative paths resolve correctly)
//
// method: gauss (default) | trunc90 | trunc95
// cone:   ak4PF (default) | ak2PF | ak3PF | ak5PF | ak6PF

#ifdef __CLING__
R__ADD_INCLUDE_PATH(include)
R__ADD_INCLUDE_PATH(cfg)
#if defined(__APPLE__)
R__LOAD_LIBRARY(lib/libl2residuals.dylib)
#else
R__LOAD_LIBRARY(lib/libl2residuals.so)
#endif
#endif

#include "TextFileWriter.h"

#ifndef __CLING__
#include <iostream>
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: runTextFile <residuals.root> <output.txt> [method] [cone]\n";
        return 1;
    }
    TString method = (argc > 3) ? argv[3] : "gauss";
    TString cone   = (argc > 4) ? argv[4] : "ak4PF";
    runTextFile(argv[1], argv[2], method, cone);
    return 0;
}
#endif
