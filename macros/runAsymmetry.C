
// USAGE
//
// Binary:
//
// ./build/bin/runAsymmetry
//   -input in.root
//   -output out.root
//   -mode triggered|non-triggered|mc
//   # -maxevents n
//   # -calibration jec|jer   (default jec)
//   # -closure true|false    (default false)
//   -config path
//
// ./build/bin/runAsymmetry args.config  # config file lines: key = value
//
// -calibration jer -closure true JER-smears MC jets (hybrid method,
// JetSmearer.h) before histogramming, using jer_closure.resolution_files/
// scale_factor_files from the config -- only valid combined with -mode mc,
// for producing the smeared-MC input to the JER SF closure check. Any other
// -calibration/-closure combination is a no-op, mirroring runPlotting's
// -closure flag having no effect outside the plot families it targets.
//
// Interpreted:
//
// export L2RESIDUALS_CONFIG=/path/to/cfg/2024ppRef.toml  # required
// root -l -b -q 'macros/runAsymmetry.C("in.root", "out.root")'

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

#include "RunAsymmetry.h"
#include "jetmet/CommandLine.h"

#ifndef __CLING__
#include <cstdlib>
#include <exception>
#include <iostream>

int main(int argc, char *argv[]) {
  static const char *const kUsage =
      "Usage: runAsymmetry -input in.root -output out.root"
      " -mode triggered|non-triggered|mc [-maxevents n]\n"
      "                     [-calibration jec|jer] [-closure true|false] "
      "-config path\n"
      "       runAsymmetry args.config   # config file lines use: key = value\n"
      "  -calibration jer -closure true: JER-smear MC jets before "
      "histogramming\n"
      "         (needs -mode mc and jer_closure.resolution_files/"
      "scale_factor_files\n"
      "         in the config). Any other combination is a no-op.\n";

  CommandLine cl;
  if (!cl.parse(argc, argv))
    return 1;

  std::string input = cl.getValue<std::string>("input");
  std::string output = cl.getValue<std::string>("output");
  std::string mode = cl.getValue<std::string>("mode");
  long long maxEvents = cl.getValue<long long>("maxevents", -1LL);
  std::string calibration =
      cl.getValue<std::string>("calibration", std::string("jec"));
  bool closure = cl.getValue<bool>("closure", false);
  std::string config = cl.getValue<std::string>("config");

  if (!cl.check()) {
    std::cerr << kUsage;
    return 1;
  }

  if (calibration != "jec" && calibration != "jer") {
    std::cerr << "ERROR: -calibration must be \"jec\" or \"jer\", got: \""
              << calibration << "\"\n"
              << kUsage;
    return 1;
  }

  setenv("L2RESIDUALS_CONFIG", config.c_str(), 1);

  const bool jerClosure = (calibration == "jer") && closure;

  try {
    runAsymmetry(input, output, mode, (Long64_t)maxEvents, jerClosure);
  } catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
  return 0;
}
#endif
