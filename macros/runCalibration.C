// CMake:       ./build/bin/runCalibration -data data.root -mc mc.root -output out.root -mode jec -config path
//              ./build/bin/runCalibration -data data.root -mc mc.root -output out.root -mode jer -config path
//              ./build/bin/runCalibration args.config  (with lines like: data = data.root)
// Interpreted: export L2RESIDUALS_CONFIG=/path/to/cfg/2024ppRef.toml  (required -- no implicit default)
//              root -l -b -q 'macros/runCalibration.C("data.root","mc.root","out.root","jec")'
//              (build the library first: cmake --build build)
//              (for interpreted ROOT, run from the repo root or set L2RESIDUALS_HOME)
//
// -mode jec  -- extract L2Residual correction (mean of A distribution)
// -mode jer  -- extract JER scale factor (stddev of A distribution)
// Mode is required -- no default. JEC and JER produce separate output files
// and are run as separate passes; the mode selection determines which
// intercept TH1Ds and R(alpha) graphs are written.
//
// Compiled arguments use JetMET's own CommandLine parser (vendored under
// external/jetmet/): "-key value" on the shell, or a leading params.config
// file with "key = value" lines. Unknown/unused options and missing required
// values are immediate CLI errors, reported together by CommandLine::check().
// config is always required; there is no default TOML. Keys are matched
// exactly as written below (case-sensitive).

#ifdef __CLING__
// clang-format off
R__ADD_INCLUDE_PATH(include)
R__ADD_INCLUDE_PATH(cfg)
R__ADD_INCLUDE_PATH(external)
#if defined(__APPLE__)
R__LOAD_LIBRARY(build/lib/libl2residuals.dylib)
#else
R__LOAD_LIBRARY(build/lib/libl2residuals.so)
#endif
// clang-format on
#endif

#include "CalibrationExtractor.h"
#include "jetmet/CommandLine.h"

#ifndef __CLING__
#include <cstdlib>
#include <iostream>
int main(int argc, char *argv[]) {
  static const char *const kUsage =
      "Usage: runCalibration -data data.root -mc mc.root -output out.root "
      "-mode jec|jer -config path\n"
      "       runCalibration args.config   # config file lines use: key = "
      "value\n"
      "  -mode jec  extract L2Residual correction (mean of A)\n"
      "  -mode jer  extract JER scale factor (stddev of A)\n";

  CommandLine cl;
  if (!cl.parse(argc, argv))
    return 1;

  std::string data = cl.getValue<std::string>("data");
  std::string mc = cl.getValue<std::string>("mc");
  std::string output = cl.getValue<std::string>("output");
  std::string config = cl.getValue<std::string>("config");
  std::string modeStr = cl.getValue<std::string>("mode");

  if (!cl.check()) {
    std::cerr << kUsage;
    return 1;
  }

  CalibrationMode mode;
  if (modeStr == "jec") {
    mode = CalibrationMode::JEC;
  } else if (modeStr == "jer") {
    mode = CalibrationMode::JER;
  } else {
    std::cerr << "ERROR: -mode must be \"jec\" or \"jer\", got: \"" << modeStr
              << "\"\n"
              << kUsage;
    return 1;
  }

  setenv("L2RESIDUALS_CONFIG", config.c_str(), 1);

  runCalibration(data, mc, output, mode);
  return 0;
}
#endif
