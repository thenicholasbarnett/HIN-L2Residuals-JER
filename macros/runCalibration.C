// CMake:       ./build/bin/runCalibration -data data.root -mc mc.root -output out.root -config path
//              ./build/bin/runCalibration args.config  (with lines like: data = data.root)
// Interpreted: export L2RESIDUALS_CONFIG=/path/to/cfg/2024ppRef.toml  (required -- no implicit default)
//              root -l -b -q 'macros/runCalibration.C("data.root","mc.root","out.root")'
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

#include "CalibrationExtractor.h"
#include "CliTokens.h"

#ifndef __CLING__
#include <cstdlib>
#include <iostream>
#include <set>
int main( int argc, char* argv[] ){
    static const char* const kUsage =
        "Usage: runCalibration [-data data.root] [-mc mc.root] [-output out.root] [-config path]\n"
        "       runCalibration args.config   # config file lines use: key = value\n"
        "       Legacy KEY=value tokens are also accepted.\n";

    const std::set<std::string> kKnownKeys = { "DATA", "MC", "OUTPUT", "CONFIG" };
    L2Cli::Tokens t = L2Cli::ParseTokens( argc, argv, kKnownKeys, kUsage );

    setenv( "L2RESIDUALS_CONFIG", t.Require( "CONFIG", kUsage ).c_str(), 1 );

    TString data   = t.Require( "DATA", kUsage );
    TString mc     = t.Require( "MC", kUsage );
    TString output = t.Require( "OUTPUT", kUsage );

    runCalibration( data, mc, output );
    return 0;
}
#endif
