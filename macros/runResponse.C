// CMake:       ./build/bin/runResponse -input mc_asymmetry.root -output response.root -config path
//              ./build/bin/runResponse args.config  (with lines like: input = mc_asymmetry.root)
// Interpreted: export L2RESIDUALS_CONFIG=/path/to/cfg/2024ppRef.toml  (required -- no implicit default)
//              root -l -b -q 'macros/runResponse.C("mc_asymmetry.root","response.root")'
//              (build the library first: cmake --build build)
//              (for interpreted ROOT, run from the repo root or set L2RESIDUALS_HOME)
//
// Compiled arguments use JetMET's own CommandLine parser (vendored under
// external/jetmet/): "-key value" on the shell, or a leading params.config
// file with "key = value" lines. Unknown/unused options and missing required
// values are immediate CLI errors, reported together by CommandLine::check().
// config is always required; there is no default TOML. Keys are matched
// exactly as written below (case-sensitive).
//
// Reads a single hadded Step 1 MC file (runAsymmetry -mode mc) -- no -data,
// no -mode, since the MC-only response THnSparses this reads don't exist in
// a data-mode Step 1 file. See include/ResponseExtractor.h for what gets
// extracted and why (JES/JER binned by gen quantities, not reco).

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

#include "ResponseExtractor.h"
#include "jetmet/CommandLine.h"

#ifndef __CLING__
#include <cstdlib>
#include <iostream>
int main(int argc, char *argv[]) {
  static const char *const kUsage =
      "Usage: runResponse -input mc_asymmetry.root -output response.root "
      "-config path\n"
      "       runResponse args.config   # config file lines use: key = value\n";

  CommandLine cl;
  if (!cl.parse(argc, argv))
    return 1;

  std::string input = cl.getValue<std::string>("input");
  std::string output = cl.getValue<std::string>("output");
  std::string config = cl.getValue<std::string>("config");

  if (!cl.check()) {
    std::cerr << kUsage;
    return 1;
  }

  setenv("L2RESIDUALS_CONFIG", config.c_str(), 1);

  runResponse(input, output);
  return 0;
}
#endif
