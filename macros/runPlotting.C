// CMake:       ./build/bin/runPlotting -input residuals.root [-outdir dir] [-flags "..."] [-closure true] [-calibration JEC|JER] -config path
//              ./build/bin/runPlotting args.config  (with lines like: flags = finals etasym)
// Interpreted: export L2RESIDUALS_CONFIG=/path/to/cfg/2024ppRef.toml  (required -- no implicit default)
//              root -l -b -q 'macros/runPlotting.C("residuals.root")'
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

#include "TFile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

#include "Binning.h"
#include "AnalysisConfig.h"
#include "jetmet/CommandLine.h"
#include "ProgressBar.h"

#include "plotting/Style.h"
#include "plotting/Utilities.h"
#include "plotting/EventQA.h"
#include "plotting/Kinematics.h"
#include "plotting/AsymmetryDistributions.h"
#include "plotting/ROverlays.h"
#include "plotting/AlphaExtrapolations.h"
#include "plotting/PtExtrapolations.h"
#include "plotting/EtaSymmetry.h"
#include "plotting/MethodComparisons.h"
#include "plotting/FinalCorrections.h"
#include "plotting/NormComparisons.h"
#include "plotting/ResponsePlots.h"
#include "plotting/Counters.h"

#include <vector>
#include <iostream>

// ============================================================
// Entry point
//
// flags has three modes:
//   (empty, the default) — curated smart default, NOT every plot family:
//     event + kinematics (tag/probe only, no inclusive jets) for Step 1;
//     adist + roverlay + alpha for Step 2 (excludes etasym/methods/normcomp,
//     which need Step 2's own output but are comparison/QA plots, not
//     day-to-day output); ptfit + finals for Step 3. "finals" is included
//     here only when the file provides it via the Step 3 corrfinal-grid
//     fallback (see PlotFinals) — i.e. only on genuine Step 3 output, never
//     via Step 2's native intercepts.
//   "all" — literal firehose: every plot family below, unconditionally,
//     including inclusive-jet kinematics and "finals" from either data
//     source. Use this to regenerate everything a file supports.
//   explicit space-separated keywords — exactly the named flags run, each
//     behaving as documented below (e.g. "kinematics" always includes
//     inclusive jets; "finals" always tries both its data sources).
//
//   "etasym"    — full-eta vs |eta| reflected symmetry check (PlotEtaSym)
//   "methods"   — method comparison: gauss vs trunc90 vs trunc95 (PlotMethodComp)
//   "finals"    — final R_MC/R_data at alpha→0, all pT slices overlaid (PlotFinals)
//                 also reads Step 3's corrfinal TH2D grid as a fallback.
//                 CLOSURE=true fixes its y-range to 0.95-1.05 with 0.99/1.01
//                 guide lines, for checking a closure pass's R_MC/R_data ~= 1.
//   "normcomp"  — direct vs kFSR-norm correction factor overlay with ratio panel (PlotNormComp)
//   "adist"     — asymmetry distributions per bin with log-y and truncation lines
//   "roverlay"  — R_data and R_MC overlay with ratio panel per alpha/pT
//   "alpha"     — kFSR-normalized alpha fit plots: R_MC/R_data(α)/R_MC/R_data(0.30) vs alpha
//   "ptfit"     — Step 3 pT-dependence fit: correction factor vs pT_avg per eta bin (PlotPtFit)
//   "kinematics"— Step-1 inclusive/tag/probe jet kinematics from runAsymmetry output
//   "event"     — Step-1 event-level QA: vz, primary vertex filter, HLT trigger
//   "response"  — runResponse JES/JER: per-bin response distributions with a
//                 Gaussian guide fit, plus JES/JER vs eta_gen and vs pT_gen
//                 summary overlays (incl/tag/probe) (PlotResponse)
//
// CALIBRATION=JEC|JER (default JEC, orthogonal to FLAGS -- same shape as
// CLOSURE=true) switches which Step 2 output etasym/methods/finals/normcomp/
// roverlay/alpha read: the default JEC (mean-derived, what Step 3 turns into
// JEC text files) or JER SF (stddev-derived, "_jer"-suffixed objects --
// see CalibrationExtractor.cxx). adist reads the same underlying A
// distribution regardless, since that QA histogram isn't duplicated between
// the two. kinematics/event/ptfit/response don't have a JEC/JER axis at all
// and ignore CALIBRATION= entirely. finals' Step 3 corrfinal-grid fallback
// has no JER equivalent yet -- CALIBRATION=JER only ever finds Step 2's
// native intercept_jer_* objects there.
// ============================================================

// Smart-default "finals" inclusion. PlotFinals draws from Step 2's native
// intercept histograms when present, falling back to Step 3's corrfinal grid
// otherwise (FinalCorrections.h). The curated default only wants it when the
// fallback path is what's actually available (a genuine Step 3 file) — not
// when Step 2's own output happens to contain intercepts too. Checked once
// against the first cone; a single pipeline run is uniform across cones.
bool WantsFinalsByDefault( TFile* fIn, const TString& cone, const BinningConfig& bins, bool useJer ){
    if( bins.ptavgSlices.empty() ) return false;
    const TString etaMode = L2Name::EtaModeKey( false );
    const TString interceptName = L2Name::ObjectName( cone, CalibKind( "intercept", useJer ),
        { etaMode, L2Name::PtKey( bins.ptavgSlices[0] ) }, { kMethodKeys[0] } );
    if( HasHAny( fIn, { cone + "/" + interceptName } ) ) return false;

    // no JER SF equivalent of the corrfinal grid exists yet -- see PlotFinals
    if( useJer ) return false;

    const TString gridName = L2Name::ObjectName( cone, "corrfinal", { etaMode }, { kMethodKeys[0] } );
    return HasHAny( fIn, { cone + "/" + gridName + "_norm", gridName + "_norm",
                           cone + "/" + gridName, gridName } );
}

void runPlotting( TString residualsFile, TString outDir = "", TString flags = "", bool isClosure = false, bool useJer = false ){
    const AnalysisConfig& cfg = Config();
    PrintConfigSummary( cfg );

    SetupPlotStyle();

    TFile* fIn = TFile::Open( residualsFile, "read" );
    if( !fIn || fIn->IsZombie() ){
        std::cerr << "Cannot open " << residualsFile << "\n";
        return;
    }

    if( outDir.IsNull() ) outDir = MakePlotDir( "plots_residuals" );
    if( gSystem->mkdir( outDir, true ) < 0 && gSystem->AccessPathName( outDir ) ){
        std::cerr << "Cannot create output directory: " << outDir << "\n";
        return;
    }

    BinningConfig bins;

    const bool doAllLiteral   = ( flags == "all" );
    const bool doSmartDefault = flags.IsNull();

    const bool doEtaSym   = doAllLiteral || flags.Contains( "etasym" );
    const bool doMethods  = doAllLiteral || flags.Contains( "methods" );
    const bool doNormComp = doAllLiteral || flags.Contains( "normcomp" );
    const bool doAdist    = doAllLiteral || doSmartDefault || flags.Contains( "adist" );
    const bool doRover    = doAllLiteral || doSmartDefault || flags.Contains( "roverlay" );
    const bool doAlpha    = doAllLiteral || doSmartDefault || flags.Contains( "alpha" );
    const bool doPtFit    = doAllLiteral || doSmartDefault || flags.Contains( "ptfit" );
    const bool doKine     = doAllLiteral || doSmartDefault || flags.Contains( "kinematics" );
    const bool doEvent    = doAllLiteral || doSmartDefault || flags.Contains( "event" );
    const bool doResponse = doAllLiteral || doSmartDefault || flags.Contains( "response" );

    // "kinematics" named explicitly (or "all") pulls in inclusive jets too;
    // the smart default sticks to tag/probe only.
    const bool kineIncludeIncl = doAllLiteral || flags.Contains( "kinematics" );

    bool doFinals = doAllLiteral || flags.Contains( "finals" );
    if( doSmartDefault && !cfg.coneLabels.empty() )
        doFinals = WantsFinalsByDefault( fIn, cfg.coneLabels[0], bins, useJer );

    int totalPlots = 0;
    if( doEvent ) totalPlots += CountEventPlots( fIn );
    for( const TString& cone : cfg.coneLabels ){
        if( doEtaSym )   totalPlots += CountEtaSymPlots( fIn, cone, bins, useJer );
        if( doMethods )  totalPlots += CountMethodCompPlots( fIn, cone, bins, false, useJer );
        if( doMethods )  totalPlots += CountMethodCompPlots( fIn, cone, bins, true, useJer );
        if( doFinals )   totalPlots += CountFinalsPlots( fIn, cone, bins, useJer );
        if( doNormComp ) totalPlots += CountNormCompPlots( fIn, cone, bins, false, useJer );
        if( doNormComp ) totalPlots += CountNormCompPlots( fIn, cone, bins, true, useJer );
        if( doAdist )    totalPlots += CountAsymDistPlots( fIn, cone, bins );
        if( doRover )    totalPlots += CountROverlayPlots( fIn, cone, bins, useJer );
        if( doAlpha )    totalPlots += CountAlphaFitPlots( fIn, cone, bins, useJer );
        if( doPtFit )    totalPlots += CountPtFitPlots( fIn, cone );
        if( doKine )     totalPlots += CountKinematicsPlots( fIn, cone, kineIncludeIncl );
        if( doResponse ) totalPlots += CountResponsePlots( fIn, cone );
    }

    if( totalPlots == 0 ){
        std::cout << "No matching plots found for " << residualsFile << " with flags \"" << flags << "\"\n";
        fIn->Close();
        return;
    }

    ProgressBar pb( "Saving plots:", totalPlots );

    if( doEvent ) PlotEvent( fIn, outDir, pb );

    for( const TString& cone : cfg.coneLabels ){
        if( doEtaSym )   PlotEtaSym( fIn, outDir, cone, bins,        pb, useJer );
        if( doMethods )  PlotMethodComp( fIn, outDir, cone, bins, false, pb, useJer );
        if( doMethods )  PlotMethodComp( fIn, outDir, cone, bins, true,  pb, useJer );
        if( doFinals )   PlotFinals( fIn, outDir, cone, bins,        pb, isClosure, useJer );
        if( doNormComp ) PlotNormComp( fIn, outDir, cone, bins, false, pb, useJer );
        if( doNormComp ) PlotNormComp( fIn, outDir, cone, bins, true,  pb, useJer );
        if( doAdist )    PlotAsymDist( fIn, outDir, cone, bins,        pb );
        if( doRover )    PlotROverlay( fIn, outDir, cone, bins,        pb, useJer );
        if( doAlpha )    PlotAlphaFit( fIn, outDir, cone, bins,        pb, useJer );
        if( doPtFit )    PlotPtFit( fIn, outDir, cone,                   pb );
        if( doKine )     PlotKinematics( fIn, outDir, cone, kineIncludeIncl, pb );
        if( doResponse ) PlotResponse( fIn, outDir, cone, cfg.responseGausFitHalfWidth, pb );
    }

    pb.Finish();
    fIn->Close();
}

#ifndef __CLING__
#include <cstdlib>
#include <iostream>
int main( int argc, char* argv[] ){
    static const char* const kUsage =
        "Usage: runPlotting -input file.root [-outdir dir] [-flags \"...\"] [-closure true] [-calibration JEC|JER] -config path\n"
        "       runPlotting args.config   # config file lines use: key = value\n"
        "  -flags: omit for the curated smart default, \"all\" for every plot\n"
        "         unconditionally, or a space-separated list of:\n"
        "         etasym methods finals normcomp adist roverlay alpha ptfit kinematics event response\n"
        "  -closure true: \"finals\" plots use a fixed 0.95-1.05 y-range with 0.99/1.01\n"
        "         guide lines instead of the auto-scaled range, for checking a\n"
        "         closure pass's R_MC/R_data ~= 1. No effect on any other flag.\n"
        "  -calibration JEC|JER (default JEC): switches etasym/methods/finals/normcomp/\n"
        "         roverlay/alpha between the mean-derived JEC output and the\n"
        "         stddev-derived JER SF output. No effect on adist/kinematics/event/\n"
        "         ptfit/response.\n";

    CommandLine cl;
    if( !cl.parse( argc, argv ) ) return 1;

    std::string input       = cl.getValue<std::string>( "input" );
    std::string outDir      = cl.getValue<std::string>( "outdir", std::string( "" ) );
    std::string flags       = cl.getValue<std::string>( "flags", std::string( "" ) );
    bool isClosure          = cl.getValue<bool>( "closure", false );
    std::string calibration = cl.getValue<std::string>( "calibration", std::string( "JEC" ) );
    std::string config      = cl.getValue<std::string>( "config" );

    if( !cl.check() ){
        std::cerr << kUsage;
        return 1;
    }

    setenv( "L2RESIDUALS_CONFIG", config.c_str(), 1 );

    if( calibration != "JEC" && calibration != "JER" ){
        std::cerr << "Invalid -calibration \"" << calibration << "\" -- must be JEC or JER\n" << kUsage;
        return 1;
    }
    bool useJer = ( calibration == "JER" );

    runPlotting( input, outDir, flags, isClosure, useJer );
    return 0;
}
#endif
