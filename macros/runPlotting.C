// Compiled:    ./bin/runPlotting <residuals.root> [out_dir] [flags] [CONFIG=path]
// Interpreted: root -l -b -q 'macros/runPlotting.C("residuals.root")'
//              (for interpreted ROOT, run from the repo root or set L2RESIDUALS_HOME)

#ifdef __CLING__
R__ADD_INCLUDE_PATH(include)
R__ADD_INCLUDE_PATH(cfg)
#if defined(__APPLE__)
R__LOAD_LIBRARY(lib/libl2residuals.dylib)
#else
R__LOAD_LIBRARY(lib/libl2residuals.so)
#endif
#endif

#include "TFile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

#include "Binning.h"
#include "AnalysisConfig.h"
#include "ConfigCli.h"
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
#include "plotting/Counters.h"

#include <vector>
#include <iostream>

// ============================================================
// Entry point
//
// flags (space-separated keywords, default "all"):
//   "etasym"    — full-eta vs |eta| reflected symmetry check (PlotEtaSym)
//   "methods"   — method comparison: gauss vs trunc90 vs trunc95 (PlotMethodComp)
//   "finals"    — final R_MC/R_data at alpha→0, all pT slices overlaid (PlotFinals)
//                 also reads Step 3's corrfinal TH2D grid as a fallback
//   "normcomp"  — direct vs kFSR-norm correction factor overlay with ratio panel (PlotNormComp)
//   "adist"     — asymmetry distributions per bin with log-y and truncation lines
//   "roverlay"  — R_data and R_MC overlay with ratio panel per alpha/pT
//   "alpha"     — kFSR-normalized alpha fit plots: R_MC/R_data(α)/R_MC/R_data(0.30) vs alpha
//   "ptfit"     — Step 3 pT-dependence fit: correction factor vs pT_avg per eta bin (PlotPtFit)
//   "kinematics"— Step-1 inclusive/tag/probe jet kinematics from runAsymmetry output
//   "event"     — Step-1 event-level QA: vz, primary vertex filter, HLT trigger
//   "all"       — run all plots including kinematics and event (default)
// ============================================================

void runPlotting( TString residualsFile, TString outDir = "", TString flags = "all" ){
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

    const bool doAll      = flags.IsNull() || flags == "all";
    const bool doEtaSym   = doAll || flags.Contains( "etasym" );
    const bool doMethods  = doAll || flags.Contains( "methods" );
    const bool doFinals   = doAll || flags.Contains( "finals" );
    const bool doNormComp = doAll || flags.Contains( "normcomp" );
    const bool doAdist    = doAll || flags.Contains( "adist" );
    const bool doRover    = doAll || flags.Contains( "roverlay" );
    const bool doAlpha    = doAll || flags.Contains( "alpha" );
    const bool doPtFit    = doAll || flags.Contains( "ptfit" );
    const bool doKine     = doAll || flags.Contains( "kinematics" );
    const bool doEvent    = doAll || flags.Contains( "event" );

    int totalPlots = 0;
    if( doEvent ) totalPlots += CountEventPlots( fIn );
    for( const TString& cone : cfg.coneLabels ){
        if( doEtaSym )   totalPlots += CountEtaSymPlots( fIn, cone, bins );
        if( doMethods )  totalPlots += CountMethodCompPlots( fIn, cone, bins, false );
        if( doMethods )  totalPlots += CountMethodCompPlots( fIn, cone, bins, true );
        if( doFinals )   totalPlots += CountFinalsPlots( fIn, cone, bins );
        if( doNormComp ) totalPlots += CountNormCompPlots( fIn, cone, bins, false );
        if( doNormComp ) totalPlots += CountNormCompPlots( fIn, cone, bins, true );
        if( doAdist )    totalPlots += CountAsymDistPlots( fIn, cone, bins );
        if( doRover )    totalPlots += CountROverlayPlots( fIn, cone, bins );
        if( doAlpha )    totalPlots += CountAlphaFitPlots( fIn, cone, bins );
        if( doPtFit )    totalPlots += CountPtFitPlots( fIn, cone );
        if( doKine )     totalPlots += CountKinematicsPlots( fIn, cone );
    }

    if( totalPlots == 0 ){
        std::cout << "No matching plots found for " << residualsFile << " with flags \"" << flags << "\"\n";
        fIn->Close();
        return;
    }

    ProgressBar pb( "Saving plots:", totalPlots );

    if( doEvent ) PlotEvent( fIn, outDir, pb );

    for( const TString& cone : cfg.coneLabels ){
        if( doEtaSym )   PlotEtaSym( fIn, outDir, cone, bins,        pb );
        if( doMethods )  PlotMethodComp( fIn, outDir, cone, bins, false, pb );
        if( doMethods )  PlotMethodComp( fIn, outDir, cone, bins, true,  pb );
        if( doFinals )   PlotFinals( fIn, outDir, cone, bins,        pb );
        if( doNormComp ) PlotNormComp( fIn, outDir, cone, bins, false, pb );
        if( doNormComp ) PlotNormComp( fIn, outDir, cone, bins, true,  pb );
        if( doAdist )    PlotAsymDist( fIn, outDir, cone, bins,        pb );
        if( doRover )    PlotROverlay( fIn, outDir, cone, bins,        pb );
        if( doAlpha )    PlotAlphaFit( fIn, outDir, cone, bins,        pb );
        if( doPtFit )    PlotPtFit( fIn, outDir, cone,                   pb );
        if( doKine )     PlotKinematics( fIn, outDir, cone,              pb );
    }

    pb.Finish();
    fIn->Close();
}

#ifndef __CLING__
#include <iostream>
int main( int argc, char* argv[] ){
    L2ConfigCli::ApplyConfigArgument( argc, argv );
    std::vector<std::string> args = L2ConfigCli::PositionalArgs( argc, argv );
    if( args.size() < 1 ){
        std::cerr << "Usage: runPlotting <residuals.root> [out_dir] [flags]"
                  << L2ConfigCli::ConfigUsage() << "\n"
                  << "  flags: all etasym methods finals normcomp adist roverlay alpha ptfit kinematics event (space-separated)\n";
        return 1;
    }
    runPlotting( args[0], args.size() >= 2 ? args[1] : "", args.size() >= 3 ? args[2] : "all" );
    return 0;
}
#endif
