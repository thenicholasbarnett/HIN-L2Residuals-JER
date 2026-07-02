// CMake:       ./build/bin/runResiduals DATA=data.root MC=mc.root OUTPUT=out.root CONFIG=path
// Interpreted: export L2RESIDUALS_CONFIG=/path/to/cfg/2024ppRef.toml  (required -- no implicit default)
//              root -l -b -q 'macros/runResiduals.C("data.root","mc.root","out.root")'
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

#include "ResidualsExtractor.h"
#include "CliTokens.h"

#ifndef __CLING__
#include <cstdlib>
#include <iostream>
#include <set>
int main( int argc, char* argv[] ){
    static const char* const kUsage =
        "Usage: runResiduals DATA=data.root MC=mc.root OUTPUT=out.root CONFIG=path\n";

    const std::set<std::string> kKnownKeys = { "DATA", "MC", "OUTPUT", "CONFIG" };
    L2Cli::Tokens t = L2Cli::ParseTokens( argc, argv, kKnownKeys, kUsage );

    setenv( "L2RESIDUALS_CONFIG", t.Require( "CONFIG", kUsage ).c_str(), 1 );

    TString data   = t.Require( "DATA", kUsage );
    TString mc     = t.Require( "MC", kUsage );
    TString output = t.Require( "OUTPUT", kUsage );

    runResiduals( data, mc, output );
    return 0;
}
#endif
