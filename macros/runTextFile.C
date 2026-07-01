// Compiled:    ./bin/runTextFile <hp_residuals.root> <zb_residuals.root> <output.root> <output_text_prefix> [method] [direct]
// Interpreted: root -l -b -q 'macros/runTextFile.C("hp.root","zb.root","out.root","corrections/hp0_zb0")'
//              (build the library first: cmake --build build)
//              (run from the repo root so relative paths resolve correctly)
//
// Processes every cone in cfg.coneLabels. Per pT_avg slice, uses hp_residuals
// if the slice starts at or above cfg.hltJ80Thresh, otherwise zb_residuals.
// Writes "<output_text_prefix>_<cone>_abseta[_norm].txt" and "..._<cone>_eta[_norm].txt".
//
// method: gauss (default) | doubleGauss | trunc90 | trunc95
// [direct]: pass the literal word "direct" to use the non-normalized intercepts
//           instead of the kFSR-normalized ones (the standard/default method).

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
int main( int argc, char* argv[] ){
    if( argc < 5 ){
        std::cerr << "Usage: runTextFile <hp_residuals.root> <zb_residuals.root> <output.root> <output_text_prefix> [method] [direct]\n";
        return 1;
    }
    TString method = ( argc > 5 ) ? argv[5] : "gauss";
    bool useNorm = !( ( argc > 6 ) && TString( argv[6] ) == "direct" );
    runTextFile( argv[1], argv[2], argv[3], argv[4], method, useNorm );
    return 0;
}
#endif
