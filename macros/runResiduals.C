// Compiled:    ./bin/runResiduals <data.root> <mc.root> <output.root> [CONFIG=path]
// Interpreted: root -l -b -q 'macros/runResiduals.C("data.root","mc.root","out.root")'
//              (build the library first: cmake --build build)
//              (for interpreted ROOT, run from the repo root or set L2RESIDUALS_HOME)

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
#include "ConfigCli.h"

#ifndef __CLING__
#include <iostream>
int main( int argc, char* argv[] ){
    L2ConfigCli::ApplyConfigArgument( argc, argv );
    std::vector<std::string> args = L2ConfigCli::PositionalArgs( argc, argv );
    if( args.size() < 3 ){
        std::cerr << "Usage: runResiduals <data.root> <mc.root> <output.root>"
                  << L2ConfigCli::ConfigUsage() << "\n";
        return 1;
    }
    runResiduals( args[0], args[1], args[2] );
    return 0;
}
#endif
