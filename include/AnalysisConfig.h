#ifndef ANALYSISCONFIG_H
#define ANALYSISCONFIG_H

#include "TString.h"

#include <string>
#include <vector>

struct AnalysisConfig {
  std::string configPath;
  std::string repoRoot;

  std::vector<std::vector<std::string>> jecFilesPerCone;
  std::vector<std::vector<std::string>> residualFilesPerCone;
  // one JetResolutionObject text file per cone (no chaining, unlike JEC/
  // residuals above) -- only read when runAsymmetry gets -calibration jer
  // -closure true, to JER-smear MC jets for the JER SF closure check
  std::vector<std::string> jerResolutionFilesPerCone;
  std::vector<std::string> jerScaleFactorFilesPerCone;
  // "hybrid" (JME default) | "scaling" | "stochastic", JetSmearer.h
  std::string jerMethod = "hybrid";
  // JetSelector inputs: correctionlib-format jet ID + veto map JSON
  std::string jetIdPath;
  std::string vetoMapPath;
  std::string jetSystem;  // "pp" | "ion"
  // "analysis" (jetvetomap) | "calibration" (jetvetomap_all)
  std::string jetPurpose;
  // cone whose jets drive the JME Run 3 event veto in runAsymmetry, "" = off
  std::string jetEventVetoCone;
  TString jsonPath;

  TString hiTreePath;
  TString ggTreePath;
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
  double residualGausFitHalfWidth = 0.5;
  double residualAlphaFitHi = 0.31;
  double responseGausFitHalfWidth = 0.3;

  TString defaultMethod = "gauss";
  TString etaModeOutput = "both";
  // Step 3 pT-fit x values: "mean" (weighted <pT_avg> per slice) | "midpoint"
  TString ptCenter = "mean";
};

std::string DefaultConfigPath();
AnalysisConfig LoadAnalysisConfig(const std::string &path = "");
const AnalysisConfig &Config();
void PrintConfigSummary(const AnalysisConfig &cfg);

#endif
