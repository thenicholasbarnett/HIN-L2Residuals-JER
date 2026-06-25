#ifndef CONFIG_2024PPREF_H
#define CONFIG_2024PPREF_H

// Run configuration for the 2024 pp reference run (387474-387721).
// Paths are for lxplus; override locally as needed.

#include "TString.h"
#include <vector>
#include <string>

// ---- jet energy corrections ----
// Applied to rawpt for every jet before any selection.
// L2Relative removes detector non-uniformities; L2Residuals is what we are measuring.
// For the asymmetry generator, include L2L3 through the current residuals iteration.
inline const std::vector<std::string> kJECFiles = {
    "/afs/cern.ch/user/n/nbarnett/public/txt_files/L2L3_ppReco_2024ppRef/2024ppRef_withPU_L2Relative_AK4PF.txt",
    "/afs/cern.ch/user/n/nbarnett/public/txt_files/L2L3_ppReco_2024ppRef/L2Residuals_2024ppRef_12_15_2025.txt",
};

// ---- veto map (DATA only) ----
inline const std::string kVetoMapPath =
    "/eos/cms/store/group/phys_heavyions/nbarnett/Winter25Prompt25_RunCDE.root";

// ---- golden JSON (DATA only) ----
inline const TString kJSONPath =
    "/eos/cms/store/group/phys_heavyions/nbarnett/JSON_files/Cert_Collisions2024_ppref_387474_387721_golden.json";

// ---- TTree paths ----
// DATA has three main trees + a separate HLT decision tree.
// MC has two trees (no filter tree, no trigger tree).
inline const TString kHiTreePath       = "hiEvtAnalyzer/HiTree";
inline const TString kSkimTreePath     = "skimanalysis/HltTree";   // ppvF filter lives here
inline const TString kJetTreePath      = "ak4PFJetAnalyzer/t";
inline const TString kTrigTreePath     = "hltanalysis/HltTree";    // HLT decisions (DATA only)
// TODO: confirm kTrigTreePath from the HiForest file — check with TFile::ls()

// ---- triggers (DATA only) ----
// Branch names in the HLT tree.  Version suffix (_v*) may change between
// runs — check the trigger menu for the 2024 pp reference period.
inline const TString kHLTJ40Branch     = "HLT_AK4PFJet40_v1";
inline const TString kHLTJ80Branch     = "HLT_AK4PFJet80_v1";
// TODO: confirm branch names from trigger menu / TTree::Print()

// Efficiency thresholds: the trigger is fully efficient above these pT values.
static constexpr float kHLTJ40Thresh  = 45.0f;   // GeV
static constexpr float kHLTJ80Thresh  = 95.0f;   // GeV

// ---- cone sizes ----
// Labels used to prefix histogram names; one ConeHistograms instance per entry.
// Extend this when adding AK3PF, AK6PF, etc.
inline const std::vector<TString> kConeLabels = { "ak4PF" };

#endif
