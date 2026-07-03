// CMake:       ./build/bin/runResponse -input mc_asymmetry.root -output response.root -config path
//              ./build/bin/runResponse args.config  (with lines like: input = mc_asymmetry.root)
// Interpreted: export L2RESIDUALS_CONFIG=/path/to/cfg/2024ppRef.toml  (required -- no implicit default)
//              root -l -b -q 'macros/runResponse.C("mc_asymmetry.root","response.root")'
//              (build the library first: cmake --build build)
//              (for interpreted ROOT, run from the repo root or set L2RESIDUALS_HOME)
//
// Compiled arguments accept JetMET-style "-key value", config files with
// "key = value", and the original KEY=value shell-token form. Unknown keys,
// malformed options, and missing required values are immediate CLI errors.
// config/CONFIG is always required; there is no default TOML.
//
// Reads a single hadded Step 1 MC file (runAsymmetry MODE=mc) -- no DATA=,
// no MODE=, since the MC-only response THnSparses this reads don't exist in
// a data-mode Step 1 file. See include/ResponseExtractor.h for what gets
// extracted and why (JES/JER binned by gen quantities, not reco).

#ifdef __CLING__
R__ADD_INCLUDE_PATH(include)
R__ADD_INCLUDE_PATH(cfg)
#if defined(__APPLE__)
R__LOAD_LIBRARY(build/lib/libl2residuals.dylib)
#else
R__LOAD_LIBRARY(build/lib/libl2residuals.so)
#endif
#endif

#include "ResponseExtractor.h"
#include "CliTokens.h"

#ifndef __CLING__
#include <cstdlib>
#include <iostream>
#include <set>
int main( int argc, char* argv[] ){
    static const char* const kUsage =
        "Usage: runResponse [-input mc_asymmetry.root] [-output response.root] [-config path]\n"
        "       runResponse args.config   # config file lines use: key = value\n"
        "       Legacy KEY=value tokens are also accepted.\n";

    const std::set<std::string> kKnownKeys = { "INPUT", "OUTPUT", "CONFIG" };
    L2Cli::Tokens t = L2Cli::ParseTokens( argc, argv, kKnownKeys, kUsage );

    setenv( "L2RESIDUALS_CONFIG", t.Require( "CONFIG", kUsage ).c_str(), 1 );

    TString input  = t.Require( "INPUT", kUsage );
    TString output = t.Require( "OUTPUT", kUsage );

    runResponse( input, output );
    return 0;
}
#endif
