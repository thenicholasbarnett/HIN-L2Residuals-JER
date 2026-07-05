
// USAGE
//
// Binary:      
//
// ./build/bin/runCalibration
//   -data data.root
//   -mc mc.root
//   -output out.root
//   -mode jec
//   -config path
//
//./build/bin/runCalibration args.config  # config file lines: key = value
//
// Interpreted:
//
// export L2RESIDUALS_CONFIG=/path/to/cfg/2024ppRef.toml  # required
// root -l -b -q 'macros/runCalibration.C("data.root", "mc.root", "out.root", "jec")'

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
