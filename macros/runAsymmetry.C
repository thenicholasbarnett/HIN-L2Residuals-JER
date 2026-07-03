// CMake:       ./build/bin/runAsymmetry -input in.root -output out.root -mode triggered|non-triggered|mc [-maxevents n] -config path
//              ./build/bin/runAsymmetry args.config  (with lines like: input = in.root)
// Interpreted: export L2RESIDUALS_CONFIG=/path/to/cfg/2024ppRef.toml  (required -- no implicit default)
//              root -l -b -q 'macros/runAsymmetry.C("in.root","out.root")'
//              (build the library first: cmake --build build)
//              (for interpreted ROOT, run from the repo root or set L2RESIDUALS_HOME)
//
// Compiled arguments accept JetMET-style "-key value", config files with
// "key = value", and the original KEY=value shell-token form. Unknown keys,
// malformed options, and missing required values are immediate CLI errors.
// config/CONFIG is always required; there is no default TOML.

#ifdef __CLING__
R__ADD_INCLUDE_PATH(include)
R__ADD_INCLUDE_PATH(cfg)
#if defined(__APPLE__)
R__LOAD_LIBRARY(build/lib/libl2residuals.dylib)
#else
R__LOAD_LIBRARY(build/lib/libl2residuals.so)
#endif
#endif

#include "RunAsymmetry.h"
#include "CliTokens.h"

#ifndef __CLING__
#include <cstdlib>
#include <exception>
#include <iostream>
#include <set>
int main( int argc, char* argv[] ){
    static const char* const kUsage =
        "Usage: runAsymmetry [-input in.root] [-output out.root]"
        " [-mode triggered|non-triggered|mc] [-maxevents n] [-config path]\n"
        "       runAsymmetry args.config   # config file lines use: key = value\n"
        "       Legacy KEY=value tokens are also accepted.\n";

    const std::set<std::string> kKnownKeys = { "INPUT", "OUTPUT", "MODE", "MAXEVENTS", "CONFIG" };
    L2Cli::Tokens t = L2Cli::ParseTokens( argc, argv, kKnownKeys, kUsage );

    setenv( "L2RESIDUALS_CONFIG", t.Require( "CONFIG", kUsage ).c_str(), 1 );

    TString input  = t.Require( "INPUT", kUsage );
    TString output = t.Require( "OUTPUT", kUsage );
    TString mode   = t.Require( "MODE", kUsage );
    Long64_t maxEvents = t.Has( "MAXEVENTS" ) ? std::atoll( t.Get( "MAXEVENTS" ).c_str() ) : -1;

    try {
        runAsymmetry( input, output, mode, maxEvents );
    } catch( const std::exception& e ){
        std::cerr << e.what() << "\n";
        return 1;
    }
    return 0;
}
#endif
