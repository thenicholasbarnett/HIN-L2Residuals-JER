// CMake:       ./build/bin/runCalibration -data data.root -mc mc.root -output out.root -config path
//              ./build/bin/runCalibration args.config  (with lines like: data = data.root)
// Interpreted: export L2RESIDUALS_CONFIG=/path/to/cfg/2024ppRef.toml  (required -- no implicit default)
//              root -l -b -q 'macros/runCalibration.C("data.root","mc.root","out.root")'
//              (build the library first: cmake --build build)
//              (for interpreted ROOT, run from the repo root or set L2RESIDUALS_HOME)
//
// Compiled arguments use JetMET's own CommandLine parser (vendored under
// external/jetmet/): "-key value" on the shell, or a leading params.config
// file with "key = value" lines. Unknown/unused options and missing required
// values are immediate CLI errors, reported together by CommandLine::check().
// config is always required; there is no default TOML. Keys are matched
// exactly as written below (case-sensitive).

#ifdef __CLING__
R__ADD_INCLUDE_PATH(include)
R__ADD_INCLUDE_PATH(cfg)
R__ADD_INCLUDE_PATH(external)
#if defined(__APPLE__)
R__LOAD_LIBRARY(build/lib/libl2residuals.dylib)
#else
R__LOAD_LIBRARY(build/lib/libl2residuals.so)
#endif
#endif

#include "CalibrationExtractor.h"
#include "jetmet/CommandLine.h"

#ifndef __CLING__
#include <cstdlib>
#include <iostream>
int main( int argc, char* argv[] ){
    static const char* const kUsage =
        "Usage: runCalibration -data data.root -mc mc.root -output out.root -config path\n"
        "       runCalibration args.config   # config file lines use: key = value\n";

    CommandLine cl;
    if( !cl.parse( argc, argv ) ) return 1;

    std::string data   = cl.getValue<std::string>( "data" );
    std::string mc     = cl.getValue<std::string>( "mc" );
    std::string output = cl.getValue<std::string>( "output" );
    std::string config = cl.getValue<std::string>( "config" );

    if( !cl.check() ){
        std::cerr << kUsage;
        return 1;
    }

    setenv( "L2RESIDUALS_CONFIG", config.c_str(), 1 );

    runCalibration( data, mc, output );
    return 0;
}
#endif
