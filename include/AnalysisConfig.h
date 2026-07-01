#ifndef ANALYSISCONFIG_H
#define ANALYSISCONFIG_H

#include "TString.h"

#include <string>
#include <vector>

struct AnalysisConfig {
    std::vector<std::vector<std::string>> jecFilesPerCone;
    std::vector<std::vector<std::string>> residualFilesPerCone;
    std::string vetoMapPath;
    TString jsonPath;

    TString hiTreePath;
    TString skimTreePath;
    TString trigTreePath;
    std::vector<TString> jetTreePaths;

    TString hltJ80Branch;
    float hltJ80Thresh = 0.0f;

    std::vector<TString> coneLabels;
    TString trigCone;

    std::vector<float> ptavgEdges;

    float minPt = 0.0f;
};

std::string DefaultConfigPath();
AnalysisConfig LoadAnalysisConfig( const std::string& path = "" );
const AnalysisConfig& Config();

#endif
