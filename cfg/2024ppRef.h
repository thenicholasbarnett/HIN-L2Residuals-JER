#ifndef CONFIG_2024PPREF_H
#define CONFIG_2024PPREF_H

#include "TString.h"
#include <vector>
#include <string>

// L2Relative
inline const std::vector<std::vector<std::string>> kJECFilesPerCone = {
    {"data/jec/Prompt24HIpp_V1_MC_L2Relative_AK2PF.txt"},
    {"data/jec/Prompt24HIpp_V1_MC_L2Relative_AK3PF.txt"},
    {"data/jec/Prompt24HIpp_V1_MC_L2Relative_AK4PF.txt"},
    {"data/jec/Prompt24HIpp_V1_MC_L2Relative_AK5PF.txt"},
    {"data/jec/Prompt24HIpp_V1_MC_L2Relative_AK6PF.txt"},
};

// jet veto map
inline const std::string kVetoMapPath = "data/veto/Winter25Prompt25_RunCDE.root";

// golden JSON
inline const TString kJSONPath = "data/json/Cert_Collisions2024_ppref_387474_387721_golden.json";

// TTrees
inline const TString kHiTreePath = "hiEvtAnalyzer/HiTree";
inline const TString kSkimTreePath = "skimanalysis/HltTree";
inline const TString kTrigTreePath = "hltanalysis/HltTree";
inline const std::vector<TString> kJetTreePaths = {
    "ak2PFJetAnalyzer/t", "ak3PFJetAnalyzer/t", "ak4PFJetAnalyzer/t", "ak5PFJetAnalyzer/t", "ak6PFJetAnalyzer/t",
};

// trigger
inline const TString  kHLTJ80Branch = "HLT_AK4PFJet80_v8";
static constexpr float kHLTJ80Thresh = 100.0f;

// input filelists (for Condor submission)
inline const TString kFilelistHP  = "data/txt/filelist_HiForest_2024ppref_DATA_HP0.txt";
inline const TString kFilelistZB  = "data/txt/filelist_HiForest_2024ppref_DATA_ZB0.txt";
inline const TString kFilelistMC  = "data/txt/filelist_HiForest_2024ppref_MC.txt";

// cone sizes
inline const std::vector<TString> kConeLabels = {"ak2PF","ak3PF","ak4PF","ak5PF","ak6PF"};
inline const TString kTrigCone = "ak4PF";

#endif
