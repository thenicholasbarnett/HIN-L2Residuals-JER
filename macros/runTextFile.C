// Compiled:    ./bin/runTextFile <triggered_residuals.root> <nontriggered_residuals.root> <output.root> <output_text_prefix> [method] [direct] [CONFIG=path]
//              ./bin/runTextFile --single <residuals.root> <output.root> <output_text_prefix> [method] [direct] [CONFIG=path]
// Interpreted: root -l -b -q 'macros/runTextFile.C("triggered.root","nontriggered.root","out.root","corrections/hp0_zb0")'
//              (build the library first: cmake --build build)
//              (for interpreted ROOT, run from the repo root or set L2RESIDUALS_HOME)
//
// Processes every cone in cfg.coneLabels. Per pT_avg slice, uses the triggered
// residuals if the slice starts at or above cfg.hltJ80Thresh, otherwise the
// non-triggered residuals.
// Writes "<output_text_prefix>_<cone>_abseta[_norm].txt" and "..._<cone>_eta[_norm].txt".
//
// --single: for a dataset with no triggered/non-triggered split (e.g. one
//           min-bias or single-trigger sample) — every pT slice reads from
//           <residuals.root>.
// method: gauss | doubleGauss | trunc90 | trunc95 (default: cfg [step3] default_method)
// [direct]: pass the literal word "direct" to use the non-normalized intercepts
//           instead of the kFSR-normalized ones (the standard/default method).
// [step3] eta_mode in the TOML ("both" | "abseta" | "eta") controls which of
// the two text files (per cone) actually get written; see TextFileWriter.cxx.

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
#include "ConfigCli.h"
#include "AnalysisConfig.h"

#ifndef __CLING__
#include <iostream>
int main( int argc, char* argv[] ){
    static const char* const kUsage =
        "Usage: runTextFile <triggered_residuals.root> <nontriggered_residuals.root> <output.root> <output_text_prefix> [method] [direct] [CONFIG=path]\n"
        "       runTextFile --single <residuals.root> <output.root> <output_text_prefix> [method] [direct] [CONFIG=path]\n";

    L2ConfigCli::ApplyConfigArgument( argc, argv );
    std::vector<std::string> args = L2ConfigCli::PositionalArgs( argc, argv );

    if( !args.empty() && TString( args[0] ) == "--single" ){
        if( args.size() < 4 ){ std::cerr << kUsage; return 1; }
        TString method = ( args.size() > 4 ) ? args[4] : Config().defaultMethod;
        bool useNorm = !( ( args.size() > 5 ) && TString( args[5] ) == "direct" );
        runTextFile( args[1], args[2], args[3], method, useNorm );
        return 0;
    }

    if( args.size() < 4 ){ std::cerr << kUsage; return 1; }
    TString method = ( args.size() > 4 ) ? args[4] : Config().defaultMethod;
    bool useNorm = !( ( args.size() > 5 ) && TString( args[5] ) == "direct" );
    runTextFile( args[0], args[1], args[2], args[3], method, useNorm );
    return 0;
}
#endif
