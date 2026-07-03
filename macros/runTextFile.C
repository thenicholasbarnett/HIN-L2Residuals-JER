// CMake:       ./build/bin/runTextFile -triggered trig.root -nontriggered notrig.root -output out.root [-prefix name] [-method gauss] [-norm true] -config path
//              ./build/bin/runTextFile -triggered trig.root -output out.root [-prefix name] [-method gauss] [-norm true] -config path
//              ./build/bin/runTextFile -nontriggered notrig.root -output out.root [-prefix name] [-method gauss] [-norm true] -config path
//              ./build/bin/runTextFile args.config  (with lines like: triggered = trig.root)
// Interpreted: export L2RESIDUALS_CONFIG=/path/to/cfg/2024ppRef.toml  (required -- no implicit default)
//              root -l -b -q 'macros/runTextFile.C("triggered.root","nontriggered.root","out.root")'
//              (build the library first: cmake --build build)
//              (for interpreted ROOT, run from the repo root or set L2RESIDUALS_HOME)
//
// Compiled arguments accept JetMET-style "-key value", config files with
// "key = value", and the original KEY=value shell-token form. Unknown keys,
// malformed options, and missing required values are immediate CLI errors.
// config/CONFIG is always required; there is no default TOML.
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
// TRIGGERED=/NONTRIGGERED=: pass both for the merge (as above). Pass only
//          one for single-dataset mode -- the two are NOT interchangeable:
//          TRIGGERED= alone means this one dataset is itself trigger-biased,
//          so pT_avg slices below cfg.hltJ80Thresh are dropped entirely (no
//          non-triggered fallback exists to fill them in). NONTRIGGERED=
//          alone means the dataset is not trigger-biased, so every pT_avg
//          slice is used unconditionally, no threshold cut. At least one of
//          the two is required.
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
        "Usage: runTextFile [-triggered trig.root] [-nontriggered notrig.root]"
        " [-output out.root] [-prefix name] [-method gauss] [-norm true] [-config path]\n"
        "       runTextFile args.config   # config file lines use: key = value\n"
        "       Pass triggered, nontriggered, or both. Legacy KEY=value tokens are also accepted.\n"
        "  prefix: plain filename prefix (no '/'), defaults to \"L2Residual\"\n";

    const std::set<std::string> kKnownKeys = {
        "TRIGGERED", "NONTRIGGERED", "OUTPUT", "PREFIX", "METHOD", "NORM", "CONFIG"
    };
    L2Cli::Tokens t = L2Cli::ParseTokens( argc, argv, kKnownKeys, kUsage );

    setenv( "L2RESIDUALS_CONFIG", t.Require( "CONFIG", kUsage ).c_str(), 1 );

    TString output = t.Require( "OUTPUT", kUsage );
    TString prefix = t.Get( "PREFIX", "" );
    TString method = t.Has( "METHOD" ) ? TString( t.Get( "METHOD" ) ) : Config().defaultMethod;
    bool useNorm = t.Get( "NORM", "true" ) != "false";

    const bool hasTrig   = t.Has( "TRIGGERED" );
    const bool hasNoTrig = t.Has( "NONTRIGGERED" );

    if( !hasTrig && !hasNoTrig ){
        std::cerr << "ERROR: pass TRIGGERED=..., NONTRIGGERED=..., or both\n" << kUsage;
        return 1;
    }

    if( hasTrig && hasNoTrig ){
        runTextFile( t.Get( "TRIGGERED" ), t.Get( "NONTRIGGERED" ), output, prefix, method, useNorm );
        return 0;
    }

    if( hasTrig ){
        runTextFile( t.Get( "TRIGGERED" ), SingleDatasetKind::Triggered, output, prefix, method, useNorm );
    } else {
        runTextFile( t.Get( "NONTRIGGERED" ), SingleDatasetKind::NonTriggered, output, prefix, method, useNorm );
    }
    return 0;
}
#endif
