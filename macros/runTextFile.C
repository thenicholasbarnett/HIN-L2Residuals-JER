
// USAGE
//
// Binary:
//
// ./build/bin/runTextFile
//   -triggered trig.root
//   -nontriggered notrig.root
//   -output out.root
//   -mode jec
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
// compatible with only triggered residuals, only non-trigger residuals, or both.
// -mode jec|jer is required whenever -triggered/-nontriggered is used, mirroring
// runCalibration's own -mode -- a Step 2 file only ever has one calibration kind.

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
      " -output out.root -mode jec|jer [-tag name] [-method gauss] [-norm "
      "true] -config path\n"
      "       runTextFile -triggered trig.root -output out.root -mode "
      "jec|jer ... -config path\n"
      "       runTextFile -nontriggered notrig.root -output out.root -mode "
      "jec|jer ... -config path\n"
      "       runTextFile -resolution response.root [-tag name] -config path\n"
      "       runTextFile args.config   # config file lines use: key = value\n"
      "  -mode jec|jer: required whenever -triggered or -nontriggered is "
      "given, mirroring\n"
      "                 runCalibration's own -mode -- a Step 2 file only "
      "ever has one\n"
      "                 calibration kind. jec writes the L2Residual "
      "correction text files;\n"
      "                 jer writes the JER SF text files.\n"
      "  -resolution: write JER pT resolution text file from a runResponse "
      "output\n"
      "               (independent of -triggered/-nontriggered/-mode; can "
      "be combined)\n"
      "  tag: plain filename prefix (no '/'), defaults to \"L2Residual\" "
      "(-mode jec) or\n"
      "       \"JER_SF\" (-mode jer)\n"
      "          or \"JER_ptresolution\" (resolution-only run)\n";

  CommandLine cl;
  if (!cl.parse(argc, argv))
    return 1;

  // config exported first, ahead of check() aggregate
  // report: config() is a process-wide singleton reads L2RESIDUALS_CONFIG 
  // must be set first
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
  std::string modeStr = cl.getValue<std::string>("mode", std::string(""));

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

  CalibrationMode mode = CalibrationMode::JEC;
  if (hasTrig || hasNoTrig) {
    if (modeStr == "jec") {
      mode = CalibrationMode::JEC;
    } else if (modeStr == "jer") {
      mode = CalibrationMode::JER;
    } else {
      std::cerr << "ERROR: -mode must be \"jec\" or \"jer\" when using "
                   "-triggered or -nontriggered, got: \""
                << modeStr << "\"\n"
                << kUsage;
      return 1;
    }
  }

  if (hasResolution) {
    runTextFilePtResolution(resolution, tag);
  }

  if (hasTrig || hasNoTrig) {
    if (hasTrig && hasNoTrig) {
      runTextFile(trig, noTrig, output, mode, tag, method, useNorm);
    } else if (hasTrig) {
      runTextFile(trig, SingleDatasetKind::Triggered, output, mode, tag,
                  method, useNorm);
    } else {
      runTextFile(noTrig, SingleDatasetKind::NonTriggered, output, mode, tag,
                  method, useNorm);
    }
  }
  return 0;
}
#endif
