// Compiled:    ./bin/runAsymmetry <input.root> <output.root> [--monte-carlo|-mc|--non-triggered|-nt|--triggered|-t] [maxEvents] [CONFIG=path]
// Interpreted: root -l -b -q 'macros/runAsymmetry.C("in.root","out.root")'
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

#include "RunAsymmetry.h"
#include "ConfigCli.h"

#ifndef __CLING__
#include <iostream>
#include <cstdlib>
#include <exception>
int main( int argc, char* argv[] ){
    L2ConfigCli::ApplyConfigArgument( argc, argv );
    std::vector<std::string> args = L2ConfigCli::PositionalArgs( argc, argv );
    if( args.size() < 2 ){
        std::cerr << "Usage: runAsymmetry <input.root> <output.root>"
                     " [--mc|--non-triggered|--triggered] [maxEvents]"
                  << L2ConfigCli::ConfigUsage() << "\n";
        return 1;
    }
    Long64_t maxEvents = ( args.size() > 3 ) ? std::atoll( args[3].c_str() ) : -1;
    try {
        runAsymmetry( args[0], args[1], args.size() > 2 ? args[2] : "--triggered", maxEvents );
    } catch( const std::exception& e ){
        std::cerr << e.what() << "\n";
        return 1;
    }
    return 0;
}
#endif
