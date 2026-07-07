
// USAGE
//
// Binary:
//
// ./build/bin/runPlotting
//   -input residuals.root
//   -config path
//   [-outdir dir]
//   [-closure true]
//   [-calibration JEC|JER]
//   [-tag name]
//   [-sample true]
//   [-flags "..."]
//
// ./build/bin/runPlotting args.config  # config file lines: key = value
//
// -flags is last on purpose: the vendored CommandLine parser greedily
// absorbs every following bare word into whichever flag precedes it until
// it hits a token that starts with '-' -- a single-value flag placed after
// -flags (or any typo'd flag missing its leading '-') silently gets eaten
// into -flags's value instead of raising an error. Put single-value flags
// like -config before -flags, or leave nothing after it.
//
// Interpreted:
//
// export L2RESIDUALS_CONFIG=/path/to/cfg/2024ppRef.toml  # required
// root -l -b -q 'macros/runPlotting.C("residuals.root")'

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

#include "TFile.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include "AnalysisConfig.h"
#include "Binning.h"
#include "ProgressBar.h"
#include "jetmet/CommandLine.h"

#include "plotting/AlphaExtrapolations.h"
#include "plotting/AsymmetryDistributions.h"
#include "plotting/Counters.h"
#include "plotting/EtaSymmetry.h"
#include "plotting/EventQA.h"
#include "plotting/FinalCorrections.h"
#include "plotting/Kinematics.h"
#include "plotting/MethodComparisons.h"
#include "plotting/NormComparisons.h"
#include "plotting/PtExtrapolations.h"
#include "plotting/ROverlays.h"
#include "plotting/ResponsePlots.h"
#include "plotting/Style.h"
#include "plotting/Utilities.h"

#include <iostream>
#include <vector>

// flags has three modes:
//   (empty, the default):
//     event + kinematics tag/probe
//     adist + roverlay + alpha
//     ptfit + finals
//
//   "all": literal firehose, every plot family below, unconditionally
//
//   "etasym"    : full-eta vs |eta| reflected
//   "methods"   : double gauss vs gauss vs trunc90 vs trunc95
//   "finals"    : R_MC/R_data at alpha→0, pT slices overlaid
//                 CLOSURE=true fixes y-range, adds guid lines
//   "normcomp"  : direct vs kFSR-norm correction factor overlay + ratio panel
//   "adist"     : asymmetry distributions
//   "roverlay"  : R_data and R_MC overlay with ratio panel per (alpha, pT_avg)
//   "alpha"     : kFSR extrapolations: R_MC/R_data(α)/R_MC/R_data(0.30) vs alpha
//   "ptfit"     : correction factor vs pT_avg per eta bin
//   "kinematics": inclusive/tag/probe jet kinematics
//   "event"     : vz, primary vertex filter, HLT trigger
//   "response"  : per-bin response distributions with gauss fit, 
//                 plus JES/JER vs eta_gen and vs pT_gen
//                 summary overlays (incl/tag/probe)
//
// CALIBRATION=JEC|JER (default JEC, orthogonal to FLAGS) 
// JEC: mean-derived
// JER SF: width-derived
//
// -sample true  # default false
// for previewing small collection of plots for cosmetic changes

// -sample key, only some representative plots for checking cosmetics
constexpr int kSampleMaxPlots = 1000;
constexpr int kSampleMaxSeconds = 60;

bool WantsFinalsByDefault(TFile *fIn, const TString &cone,
                          const BinningConfig &bins, bool useJer) {
  if (bins.ptavgSlices.empty())
    return false;
  const TString etaMode = L2Name::EtaModeKey(false);
  const TString interceptName = L2Name::ObjectName(
      cone, CalibKind("intercept", useJer),
      {etaMode, L2Name::PtKey(bins.ptavgSlices[0])}, {kMethodKeys[0]});
  if (HasHAny(fIn, {cone + "/" + interceptName}))
    return false;

  // no JER SF equivalent of the corrfinal grid exists yet -- see PlotFinals
  if (useJer)
    return false;

  const TString gridName =
      L2Name::ObjectName(cone, "corrfinal", {etaMode}, {kMethodKeys[0]});
  return HasHAny(fIn, {cone + "/" + gridName + "_norm", gridName + "_norm",
                       cone + "/" + gridName, gridName});
}

void runPlotting(TString residualsFile, TString outDir = "", TString flags = "",
                 bool isClosure = false, bool useJer = false,
                 TString tag = "", bool sample = false) {
  const AnalysisConfig &cfg = Config();
  PrintConfigSummary(cfg);

  SetupPlotStyle();

  TFile *fIn = TFile::Open(residualsFile, "read");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "Cannot open " << residualsFile << "\n";
    return;
  }

  if (outDir.IsNull())
    outDir = MakePlotDir("plots_residuals");
  if (!tag.IsNull()) {
    if (tag.Contains("/")) {
      std::cerr << "ERROR: -tag must not contain '/': " << tag << "\n";
      return;
    }
    outDir = outDir + "/" + tag;
  }
  if (gSystem->mkdir(outDir, true) < 0 && gSystem->AccessPathName(outDir)) {
    std::cerr << "Cannot create output directory: " << outDir << "\n";
    return;
  }

  BinningConfig bins;
  
  const int minEntriesPlot =
      (cfg.minEntriesPerBin > 0) ? cfg.minEntriesPerBin : 100;

  const bool doAllLiteral = (flags == "all");
  const bool doSmartDefault = flags.IsNull();

  const bool doEtaSym = doAllLiteral || flags.Contains("etasym");
  const bool doMethods = doAllLiteral || flags.Contains("methods");
  const bool doNormComp = doAllLiteral || flags.Contains("normcomp");
  const bool doAdist =
      doAllLiteral || doSmartDefault || flags.Contains("adist");
  const bool doRover =
      doAllLiteral || doSmartDefault || flags.Contains("roverlay");
  const bool doAlpha =
      doAllLiteral || doSmartDefault || flags.Contains("alpha");
  const bool doPtFit =
      doAllLiteral || doSmartDefault || flags.Contains("ptfit");
  const bool doKine =
      doAllLiteral || doSmartDefault || flags.Contains("kinematics");
  const bool doEvent =
      doAllLiteral || doSmartDefault || flags.Contains("event");
  const bool doResponse =
      doAllLiteral || doSmartDefault || flags.Contains("response");

  // kinematics
  const bool kineIncludeIncl = doAllLiteral || flags.Contains("kinematics");

  bool doFinals = doAllLiteral || flags.Contains("finals");
  if (doSmartDefault && !cfg.coneLabels.empty())
    doFinals = WantsFinalsByDefault(fIn, cfg.coneLabels[0], bins, useJer);

  int totalPlots = 0;
  if (doEvent)
    totalPlots += CountEventPlots(fIn);
  for (const TString &cone : cfg.coneLabels) {
    if (doEtaSym)
      totalPlots += CountEtaSymPlots(fIn, cone, bins, useJer);
    if (doMethods)
      totalPlots += CountMethodCompPlots(fIn, cone, bins, false, useJer);
    if (doMethods)
      totalPlots += CountMethodCompPlots(fIn, cone, bins, true, useJer);
    if (doFinals)
      totalPlots += CountFinalsPlots(fIn, cone, bins, useJer);
    if (doNormComp)
      totalPlots += CountNormCompPlots(fIn, cone, bins, false, useJer);
    if (doNormComp)
      totalPlots += CountNormCompPlots(fIn, cone, bins, true, useJer);
    if (doAdist)
      totalPlots += CountAsymDistPlots(fIn, cone, bins, minEntriesPlot);
    if (doRover)
      totalPlots += CountROverlayPlots(fIn, cone, bins, useJer);
    if (doAlpha)
      totalPlots += CountAlphaFitPlots(fIn, cone, bins, useJer);
    if (doPtFit)
      totalPlots += CountPtFitPlots(fIn, cone);
    if (doKine)
      totalPlots += CountKinematicsPlots(fIn, cone, kineIncludeIncl);
    if (doResponse)
      totalPlots += CountResponsePlots(fIn, cone, minEntriesPlot);
  }

  if (totalPlots == 0) {
    std::cout << "No matching plots found for " << residualsFile
              << " with flags \"" << flags << "\"\n";
    fIn->Close();
    return;
  }

  ProgressBar pb("Saving plots:", totalPlots);
  if (sample)
    pb.EnableSample(kSampleMaxPlots, kSampleMaxSeconds);

  bool sampleStopped = false;
  try {
    if (doEvent)
      PlotEvent(fIn, outDir, pb);

    for (const TString &cone : cfg.coneLabels) {
      if (doEtaSym)
        PlotEtaSym(fIn, outDir, cone, bins, pb, useJer);
      if (doMethods)
        PlotMethodComp(fIn, outDir, cone, bins, false, pb, useJer);
      if (doMethods)
        PlotMethodComp(fIn, outDir, cone, bins, true, pb, useJer);
      if (doFinals)
        PlotFinals(fIn, outDir, cone, bins, pb, isClosure, useJer);
      if (doNormComp)
        PlotNormComp(fIn, outDir, cone, bins, false, pb, useJer);
      if (doNormComp)
        PlotNormComp(fIn, outDir, cone, bins, true, pb, useJer);
      if (doAdist)
        PlotAsymDist(fIn, outDir, cone, bins, minEntriesPlot, pb);
      if (doRover)
        PlotROverlay(fIn, outDir, cone, bins, pb, useJer);
      if (doAlpha)
        PlotAlphaFit(fIn, outDir, cone, bins, pb, useJer);
      if (doPtFit)
        PlotPtFit(fIn, outDir, cone, pb);
      if (doKine)
        PlotKinematics(fIn, outDir, cone, kineIncludeIncl, pb);
      if (doResponse)
        PlotResponse(fIn, outDir, cone, cfg.responseGausFitHalfWidth,
                     minEntriesPlot, pb);
    }
  } catch (const SampleLimitReached &) {
    sampleStopped = true;
    std::cout << "\n-sample true: stopped after " << pb.current << "/"
              << totalPlots << " plots (cap: " << kSampleMaxPlots
              << " plots or " << kSampleMaxSeconds << "s, whichever first)\n";
  }

  // Finish() forces current=total before its final draw
  if (sampleStopped) {
    printf("\n");
    fflush(stdout);
  } else {
    pb.Finish();
  }
  fIn->Close();
}

#ifndef __CLING__
#include <cstdlib>
#include <iostream>
int main(int argc, char *argv[]) {
  static const char *const kUsage =
      "Usage: runPlotting -input file.root -config path [-outdir dir] "
      "[-closure true]\n"
      "                    [-calibration JEC|JER] [-tag name] [-sample "
      "true] [-flags \"...\"]\n"
      "       runPlotting args.config   # config file lines use: key = value\n"
      "  IMPORTANT: put -flags LAST. The vendored CommandLine parser "
      "greedily\n"
      "         absorbs every following bare word into whichever flag "
      "precedes it,\n"
      "         until it hits a token starting with '-'. A single-value "
      "flag placed\n"
      "         after -flags -- or any flag typo'd without its leading '-' "
      "-- is\n"
      "         silently swallowed into -flags's value instead of raising "
      "an error.\n"
      "  -flags: omit for the curated smart default, \"all\" for every plot\n"
      "         unconditionally, or a space-separated list of:\n"
      "         etasym methods finals normcomp adist roverlay alpha ptfit "
      "kinematics event response\n"
      "  -closure true: \"finals\" plots use a fixed 0.95-1.05 y-range with "
      "0.99/1.01\n"
      "         guide lines instead of the auto-scaled range, for checking a\n"
      "         closure pass's R_MC/R_data ~= 1. No effect on any other flag.\n"
      "  -calibration JEC|JER (default JEC): switches "
      "etasym/methods/finals/normcomp/\n"
      "         roverlay/alpha between the mean-derived JEC output and the\n"
      "         stddev-derived JER SF output. No effect on "
      "adist/kinematics/event/\n"
      "         ptfit/response.\n"
      "  -tag: optional, plain name (no '/'). Plots land in outdir/tag/ "
      "instead of\n"
      "         outdir/ directly, so separate runs against the same outdir "
      "don't\n"
      "         overwrite each other's PNGs.\n"
      "  -sample true (default false): stop after 1000 plots or 60s, "
      "whichever\n"
      "         comes first, instead of running the full -flags selection. "
      "For\n"
      "         previewing a cosmetic change against a few real plots.\n";

  CommandLine cl;
  if (!cl.parse(argc, argv))
    return 1;

  std::string input = cl.getValue<std::string>("input");
  std::string outDir = cl.getValue<std::string>("outdir", std::string(""));
  std::string flags = cl.getValue<std::string>("flags", std::string(""));
  bool isClosure = cl.getValue<bool>("closure", false);
  std::string calibration =
      cl.getValue<std::string>("calibration", std::string("JEC"));
  std::string tag = cl.getValue<std::string>("tag", std::string(""));
  bool sample = cl.getValue<bool>("sample", false);
  std::string config = cl.getValue<std::string>("config");

  if (!cl.check()) {
    std::cerr << kUsage;
    return 1;
  }

  setenv("L2RESIDUALS_CONFIG", config.c_str(), 1);

  if (calibration != "JEC" && calibration != "JER") {
    std::cerr << "Invalid -calibration \"" << calibration
              << "\" -- must be JEC or JER\n"
              << kUsage;
    return 1;
  }
  bool useJer = (calibration == "JER");

  runPlotting(input, outDir, flags, isClosure, useJer, tag, sample);
  return 0;
}
#endif
