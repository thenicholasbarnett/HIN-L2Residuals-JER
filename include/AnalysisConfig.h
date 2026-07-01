#ifndef ANALYSISCONFIG_H
#define ANALYSISCONFIG_H

#include "TString.h"

#include <string>
#include <vector>

struct AnalysisConfig {
    std::string configPath;
    std::string repoRoot;
    bool usedDefaultConfig = false;

    std::vector<std::vector<std::string>> jecFilesPerCone;
    std::vector<std::vector<std::string>> residualFilesPerCone;
    std::string vetoMapPath;
    std::string vetoMapHist;
    TString jsonPath;

    TString hiTreePath;
    TString skimTreePath;
    TString trigTreePath;
    std::vector<TString> jetTreePaths;
    TString filterBranch;

    TString hltJ80Branch;
    float hltJ80Thresh = 0.0f;

    std::vector<TString> coneLabels;
    TString trigCone;

    std::vector<float> ptavgEdges;

    float minJetPt = 0.0f;
    float dphiCut = 0.0f;
    float maxAbsA = 0.0f;
    int minEntriesPerBin = 0;
};

std::string DefaultConfigPath();
AnalysisConfig LoadAnalysisConfig( const std::string& path = "" );
const AnalysisConfig& Config();
void PrintConfigSummary( const AnalysisConfig& cfg );

#endif
