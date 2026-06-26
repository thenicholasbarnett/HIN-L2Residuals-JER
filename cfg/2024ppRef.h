#ifndef CONFIG_2024PPREF_H
#define CONFIG_2024PPREF_H

// Run configuration for the 2024 pp reference run (387474-387721).
// Local paths point to data/ — stage files there before running locally.
// EOS/AFS paths are noted in comments for lxplus/Condor use.

#include "TString.h"
#include <vector>
#include <string>

// ---- jet energy corrections ----
// Applied to rawpt for every jet before any selection.
// L2Relative removes detector non-uniformities; L2Residuals is what we are measuring.
// For the asymmetry generator, include L2L3 through the current residuals iteration.
// EOS: /afs/cern.ch/user/n/nbarnett/public/txt_files/L2L3_ppReco_2024ppRef/
inline const std::vector<std::string> kJECFiles = {
    "data/jec/2024ppRef_withPU_L2Relative_AK4PF.txt",
    "data/jec/L2Residuals_2024ppRef_12_15_2025.txt",
};

// ---- veto map (DATA only) ----
// EOS: /eos/cms/store/group/phys_heavyions/nbarnett/Winter25Prompt25_RunCDE.root
inline const std::string kVetoMapPath = "data/veto/Winter25Prompt25_RunCDE.root";

// ---- golden JSON (DATA only) ----
// EOS: /eos/cms/store/group/phys_heavyions/nbarnett/JSON_files/Cert_Collisions2024_ppref_387474_387721_golden.json
inline const TString kJSONPath = "data/json/Cert_Collisions2024_ppref_387474_387721_golden.json";

// ---- TTree paths ----
// DATA has three main trees + a separate HLT decision tree.
// MC has two trees (no filter tree, no trigger tree).
inline const TString kHiTreePath   = "hiEvtAnalyzer/HiTree";
inline const TString kSkimTreePath = "skimanalysis/HltTree";   // ppvF filter (DATA only)
inline const TString kTrigTreePath = "hltanalysis/HltTree";    // HLT decisions (HardProbes only)
// TODO: confirm kTrigTreePath from the HiForest file — check with TFile::ls()

// Jet tree paths indexed by cone size — parallel to kConeLabels
inline const std::vector<TString> kJetTreePaths = {
    "ak2PFJetAnalyzer/t","ak3PFJetAnalyzer/t","ak4PFJetAnalyzer/t","ak5PFJetAnalyzer/t","ak6PFJetAnalyzer/t",
};

// ---- trigger (HardProbes mode only) ----
// Branch name in the HLT tree.  Version suffix (_v*) may change between
// runs — check the trigger menu for the 2024 pp reference period.
inline const TString kHLTJ80Branch = "HLT_AK4PFJet80_v1";
// TODO: confirm branch name from trigger menu / TTree::Print()

// Leading jet pT must exceed this threshold for the trigger to be fully efficient.
static constexpr float kHLTJ80Thresh = 95.0f;   // GeV

// ---- cone sizes ----
// Labels used to prefix histogram names; one ConeHistograms instance per entry.
// Extend this when adding more cone sizes; kJetTreePaths must stay parallel.
inline const std::vector<TString> kConeLabels = { "ak4PF" };

#endif
