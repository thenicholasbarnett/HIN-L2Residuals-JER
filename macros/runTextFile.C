
// USAGE
//
// Binary:
//
// ./build/bin/runTextFile
//   -triggered trig.root
//   -nontriggered notrig.root
//   -output out.root
//   # [-tag name]
//   # [-method gauss]
//   # [-norm true]
//   -config path
//
// ./build/bin/runTextFile args.config  # config file lines: key = value
//
// Interpreted: 
//
// export L2RESIDUALS_CONFIG=/path/to/cfg/2024ppRef.toml  # required
// root -l -b -q 'macros/runTextFile.C("triggered.root", "nontriggered.root", "out.root")'
//
// compatible with only triggered residuals, only non-trigger residuals, or both

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

#include "AnalysisConfig.h"
#include "TextFileWriter.h"
#include "jetmet/CommandLine.h"

#ifndef __CLING__
#include <cstdlib>
#include <iostream>
int main(int argc, char *argv[]) {
  static const char *const kUsage =
      "Usage: runTextFile -triggered trig.root -nontriggered notrig.root"
      " -output out.root [-tag name] [-method gauss] [-norm true] -config "
      "path\n"
      "       runTextFile -triggered trig.root -output out.root ... -config "
      "path\n"
      "       runTextFile -nontriggered notrig.root -output out.root ... "
      "-config path\n"
      "       runTextFile -resolution response.root [-tag name] -config path\n"
      "       runTextFile args.config   # config file lines use: key = value\n"
      "  -resolution: write JER pT resolution text file from a runResponse "
      "output\n"
      "               (independent of -triggered/-nontriggered; can be "
      "combined)\n"
      "  tag: plain filename prefix (no '/'), defaults to \"L2Residual\" "
      "(JEC/JER SF)\n"
      "          or \"JER_ptresolution\" (resolution-only run)\n";

  CommandLine cl;
  if (!cl.parse(argc, argv))
    return 1;

  // config is queried and exported first, ahead of check()'s aggregate
  // report: Config() (called just below for -method's default) is a
  // process-wide singleton that reads L2RESIDUALS_CONFIG, so it must be
  // set before any other default-value expression can safely touch it.
  std::string config = cl.getValue<std::string>("config");
  if (config.empty()) {
    std::cerr << "ERROR: missing required config=... argument\n" << kUsage;
    return 1;
  }
  setenv("L2RESIDUALS_CONFIG", config.c_str(), 1);

  std::string output = cl.getValue<std::string>("output", std::string(""));
  std::string tag = cl.getValue<std::string>("tag", std::string(""));
  std::string method = cl.getValue<std::string>(
      "method", std::string(Config().defaultMethod.Data()));
  bool useNorm = cl.getValue<bool>("norm", true);
  std::string trig = cl.getValue<std::string>("triggered", std::string(""));
  std::string noTrig =
      cl.getValue<std::string>("nontriggered", std::string(""));
  std::string resolution =
      cl.getValue<std::string>("resolution", std::string(""));

  if (!cl.check()) {
    std::cerr << kUsage;
    return 1;
  }

  const bool hasTrig = !trig.empty();
  const bool hasNoTrig = !noTrig.empty();
  const bool hasResolution = !resolution.empty();

  if (!hasTrig && !hasNoTrig && !hasResolution) {
    std::cerr << "ERROR: pass -triggered ..., -nontriggered ..., -resolution "
                 "..., or a combination\n"
              << kUsage;
    return 1;
  }

  if ((hasTrig || hasNoTrig) && output.empty()) {
    std::cerr
        << "ERROR: -output is required when using -triggered or -nontriggered\n"
        << kUsage;
    return 1;
  }

  if (hasResolution) {
    runTextFilePtResolution(resolution, tag);
  }

  if (hasTrig || hasNoTrig) {
    if (hasTrig && hasNoTrig) {
      runTextFile(trig, noTrig, output, tag, method, useNorm);
    } else if (hasTrig) {
      runTextFile(trig, SingleDatasetKind::Triggered, output, tag, method,
                  useNorm);
    } else {
      runTextFile(noTrig, SingleDatasetKind::NonTriggered, output, tag, method,
                  useNorm);
    }
  }
  return 0;
}
#endif
