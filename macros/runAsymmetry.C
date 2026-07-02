// CMake:       ./build/bin/runAsymmetry INPUT=in.root OUTPUT=out.root MODE=triggered|non-triggered|mc [MAXEVENTS=n] CONFIG=path
// Interpreted: export L2RESIDUALS_CONFIG=/path/to/cfg/2024ppRef.toml  (required -- no implicit default)
//              root -l -b -q 'macros/runAsymmetry.C("in.root","out.root")'
//              (build the library first: cmake --build build)
//              (for interpreted ROOT, run from the repo root or set L2RESIDUALS_HOME)
//
// Every argument is a KEY=value token -- there are no positional arguments,
// and none may be misspelled or omitted silently: an unknown token, a
// malformed token, or a missing required token is an immediate CLI error.
// CONFIG is always required; there is no default TOML.

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
        "Usage: runAsymmetry INPUT=in.root OUTPUT=out.root MODE=triggered|non-triggered|mc"
        " [MAXEVENTS=n] CONFIG=path\n";

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
