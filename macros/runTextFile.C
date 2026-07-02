// CMake:       ./build/bin/runTextFile TRIGGERED=trig.root NONTRIGGERED=notrig.root OUTPUT=out.root [PREFIX=name] [METHOD=gauss] [NORM=true] CONFIG=path
//              ./build/bin/runTextFile SINGLE=residuals.root OUTPUT=out.root [PREFIX=name] [METHOD=gauss] [NORM=true] CONFIG=path
// Interpreted: export L2RESIDUALS_CONFIG=/path/to/cfg/2024ppRef.toml  (required -- no implicit default)
//              root -l -b -q 'macros/runTextFile.C("triggered.root","nontriggered.root","out.root")'
//              (build the library first: cmake --build build)
//              (for interpreted ROOT, run from the repo root or set L2RESIDUALS_HOME)
//
// Every argument is a KEY=value token -- there are no positional arguments,
// and none may be misspelled or omitted silently: an unknown token, a
// malformed token, or a missing required token is an immediate CLI error.
// CONFIG is always required; there is no default TOML.
//
// Processes every cone in cfg.coneLabels. Per pT_avg slice, uses the
// triggered residuals if the slice starts at or above cfg.hltJ80Thresh,
// otherwise the non-triggered residuals.
// Correction text files always go to data/jec/preliminary/ (relative to the
// repo root, created automatically) -- gitignored, meant for locally
// generated/preliminary corrections. Writes
// "data/jec/preliminary/<prefix>_<cone>_abseta[_norm].txt" and
// "..._<cone>_eta[_norm].txt".
//
// SINGLE=: for a dataset with no triggered/non-triggered split (e.g. one
//          min-bias or single-trigger sample) -- every pT slice reads from
//          the one file. Mutually exclusive with TRIGGERED=/NONTRIGGERED=.
// PREFIX=: a plain filename prefix, NOT a path -- must not contain '/'.
//          Optional; defaults to "L2Residual" when omitted.
// METHOD=: gauss | doubleGauss | trunc90 | trunc95 (default: cfg [step3] default_method)
// NORM=:   "true" (default) uses the kFSR-normalized intercepts (the standard
//          method); "false" uses the direct, non-normalized variant instead.
// [step3] eta_mode in the TOML ("both" | "abseta" | "eta") controls which of
// the two text files (per cone) actually get written; see TextFileWriter.cxx.

#ifdef __CLING__
R__ADD_INCLUDE_PATH(include)
R__ADD_INCLUDE_PATH(cfg)
#if defined(__APPLE__)
R__LOAD_LIBRARY(build/lib/libl2residuals.dylib)
#else
R__LOAD_LIBRARY(build/lib/libl2residuals.so)
#endif
#endif

#include "TextFileWriter.h"
#include "CliTokens.h"
#include "AnalysisConfig.h"

#ifndef __CLING__
#include <iostream>
#include <set>
int main( int argc, char* argv[] ){
    static const char* const kUsage =
        "Usage: runTextFile TRIGGERED=trig.root NONTRIGGERED=notrig.root OUTPUT=out.root"
        " [PREFIX=name] [METHOD=gauss] [NORM=true] CONFIG=path\n"
        "       runTextFile SINGLE=residuals.root OUTPUT=out.root"
        " [PREFIX=name] [METHOD=gauss] [NORM=true] CONFIG=path\n"
        "  PREFIX: plain filename prefix (no '/'), defaults to \"L2Residual\"\n";

    const std::set<std::string> kKnownKeys = {
        "TRIGGERED", "NONTRIGGERED", "SINGLE", "OUTPUT", "PREFIX", "METHOD", "NORM", "CONFIG"
    };
    L2Cli::Tokens t = L2Cli::ParseTokens( argc, argv, kKnownKeys, kUsage );

    setenv( "L2RESIDUALS_CONFIG", t.Require( "CONFIG", kUsage ).c_str(), 1 );

    TString output = t.Require( "OUTPUT", kUsage );
    TString prefix = t.Get( "PREFIX", "" );
    TString method = t.Has( "METHOD" ) ? TString( t.Get( "METHOD" ) ) : Config().defaultMethod;
    bool useNorm = t.Get( "NORM", "true" ) != "false";

    const bool hasSingle = t.Has( "SINGLE" );
    const bool hasSplit  = t.Has( "TRIGGERED" ) || t.Has( "NONTRIGGERED" );

    if( hasSingle && hasSplit ){
        std::cerr << "ERROR: pass either SINGLE=... or TRIGGERED=.../NONTRIGGERED=..., not both\n" << kUsage;
        return 1;
    }

    if( hasSingle ){
        runTextFile( t.Get( "SINGLE" ), output, prefix, method, useNorm );
        return 0;
    }

    TString triggered    = t.Require( "TRIGGERED", kUsage );
    TString nontriggered = t.Require( "NONTRIGGERED", kUsage );
    runTextFile( triggered, nontriggered, output, prefix, method, useNorm );
    return 0;
}
#endif
