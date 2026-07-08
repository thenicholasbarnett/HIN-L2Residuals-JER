#include "AnalysisConfig.h"

#include "toml.hpp"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {

namespace fs = std::filesystem;

bool IsAbsolutePath(const std::string &path) {
  return !path.empty() && fs::path(path).is_absolute();
}

std::string RepoRoot() {
  const char *envHome = std::getenv("L2RESIDUALS_HOME");
  if (envHome && envHome[0])
    return envHome;

#ifdef L2RESIDUALS_SOURCE_DIR
  return L2RESIDUALS_SOURCE_DIR;
#else
  return ".";
#endif
}

std::string ResolvePath(const std::string &path, const std::string &root) {
  if (path.empty() || IsAbsolutePath(path))
    return path;
  return (fs::path(root) / path).lexically_normal().string();
}

std::string ResolveConfigPath(const std::string &path) {
  if (path.empty())
    return path;
  return ResolvePath(path, RepoRoot());
}

std::string ConfigRoot(const std::string &configPath) {
  const char *envHome = std::getenv("L2RESIDUALS_HOME");
  if (envHome && envHome[0])
    return envHome;

  if (IsAbsolutePath(configPath)) {
    fs::path configDir = fs::path(configPath).parent_path();
    if (configDir.filename() == "cfg")
      return configDir.parent_path().string();
    return configDir.string();
  }

  return RepoRoot();
}

TString ResolvePath(const TString &path, const std::string &root) {
  return TString(ResolvePath(std::string(path.Data()), root).c_str());
}

} // namespace

// main calibration vs closure config w/ different residual_files load
// L2RESIDUALS_HOME is within TOML (see ConfigRoot() below)
std::string DefaultConfigPath() {
  const char *envConfig = std::getenv("L2RESIDUALS_CONFIG");
  if (envConfig && envConfig[0])
    return ResolveConfigPath(envConfig);

  throw std::runtime_error("No config: set L2RESIDUALS_CONFIG or "
                           "pass CONFIG=path to binaries or "
                           "set path to LoadAnalysisConfig()");
}

AnalysisConfig LoadAnalysisConfig(const std::string &path) {
  const std::string configPath =
      path.empty() ? DefaultConfigPath() : ResolveConfigPath(path);
  const auto doc = toml::parse_file(configPath);
  const std::string repoRoot = ConfigRoot(configPath);

  AnalysisConfig cfg;
  cfg.configPath = configPath;
  cfg.repoRoot = repoRoot;

  cfg.vetoMapPath =
      ResolvePath(doc["paths"]["veto_map"].value_or(std::string{}), repoRoot);
  cfg.jsonPath = ResolvePath(
      TString(doc["paths"]["golden_json"].value_or(std::string{}).c_str()),
      repoRoot);

  cfg.vetoMapHist = doc["jet_id"]["veto_map_histogram"].value_or(std::string{});

  cfg.hiTreePath = TString(doc["trees"]["hi"].value_or(std::string{}).c_str());
  cfg.ggTreePath = TString(doc["trees"]["gg"].value_or(std::string{}).c_str());
  cfg.skimTreePath =
      TString(doc["trees"]["skim"].value_or(std::string{}).c_str());
  cfg.trigTreePath =
      TString(doc["trees"]["trigger"].value_or(std::string{}).c_str());
  cfg.filterBranch =
      TString(doc["trees"]["filter"].value_or(std::string{}).c_str());

  for (const auto &v : *doc["trees"]["jets"].as_array())
    cfg.jetTreePaths.push_back(TString(v.value_or(std::string{}).c_str()));

  cfg.hltJ80Branch =
      TString(doc["trigger"]["branch"].value_or(std::string{}).c_str());
  cfg.hltJ80Thresh = doc["trigger"]["threshold"].value_or(0.0f);
  cfg.trigCone =
      TString(doc["trigger"]["cone"].value_or(std::string{}).c_str());

  for (const auto &v : *doc["cones"]["labels"].as_array())
    cfg.coneLabels.push_back(TString(v.value_or(std::string{}).c_str()));

  for (const auto &row : *doc["jec"]["files"].as_array()) {
    std::vector<std::string> files;
    for (const auto &v : *row.as_array())
      files.push_back(ResolvePath(v.value_or(std::string{}), repoRoot));
    cfg.jecFilesPerCone.push_back(std::move(files));
  }

  // Optional: absent or empty, both mean "no residuals for this cone yet"
  if (auto *residualArr = doc["jec"]["residual_files"].as_array()) {
    for (const auto &row : *residualArr) {
      std::vector<std::string> files;
      for (const auto &v : *row.as_array())
        files.push_back(ResolvePath(v.value_or(std::string{}), repoRoot));
      cfg.residualFilesPerCone.push_back(std::move(files));
    }
  }

  // Optional: absent or empty means "no JER closure smearing configured yet"
  if (auto *jerResArr = doc["jer_closure"]["resolution_files"].as_array()) {
    for (const auto &v : *jerResArr)
      cfg.jerResolutionFilesPerCone.push_back(
          ResolvePath(v.value_or(std::string{}), repoRoot));
  }
  if (auto *jerSfArr = doc["jer_closure"]["scale_factor_files"].as_array()) {
    for (const auto &v : *jerSfArr)
      cfg.jerScaleFactorFilesPerCone.push_back(
          ResolvePath(v.value_or(std::string{}), repoRoot));
  }

  for (const auto &v : *doc["binning"]["ptavg_edges"].as_array())
    cfg.ptavgEdges.push_back(v.value_or(0.0f));
    
  cfg.minJetPt = doc["cuts"]["min_jet_pt"].value_or(0.0f);
  cfg.dphiCut = doc["cuts"]["dphi"].value_or(0.0f);
  cfg.maxAbsA = doc["cuts"]["max_abs_a"].value_or(0.0f);
  cfg.minEntriesPerBin = doc["cuts"]["min_entries_per_bin"].value_or(0);
  cfg.residualGausFitHalfWidth =
      doc["cuts"]["residual_gaus_fit_half_width"].value_or(0.5);
  cfg.residualAlphaFitHi = doc["cuts"]["residual_alpha_fit_hi"].value_or(0.31);
  cfg.responseGausFitHalfWidth =
      doc["cuts"]["response_gaus_fit_half_width"].value_or(0.3);
  cfg.responseTruncFraction =
      doc["cuts"]["response_trunc_fraction"].value_or(0.95);

  // Step 3 output selection
  cfg.defaultMethod = TString(
      doc["step3"]["default_method"].value_or(std::string("gauss")).c_str());
  cfg.etaModeOutput =
      TString(doc["step3"]["eta_mode"].value_or(std::string("both")).c_str());

  if (cfg.jecFilesPerCone.empty())
    throw std::runtime_error("jec.files is empty");
  if (cfg.coneLabels.empty())
    throw std::runtime_error("cones.labels is empty");
  if (cfg.jecFilesPerCone.size() != cfg.coneLabels.size())
    throw std::runtime_error(
        "jec.files and cones.labels have different lengths");
  if (!cfg.residualFilesPerCone.empty() &&
      cfg.residualFilesPerCone.size() != cfg.coneLabels.size())
    throw std::runtime_error(
        "jec.residual_files and cones.labels have different lengths");
  if (!cfg.jerResolutionFilesPerCone.empty() &&
      cfg.jerResolutionFilesPerCone.size() != cfg.coneLabels.size())
    throw std::runtime_error(
        "jer_closure.resolution_files and cones.labels have different lengths");
  if (!cfg.jerScaleFactorFilesPerCone.empty() &&
      cfg.jerScaleFactorFilesPerCone.size() != cfg.coneLabels.size())
    throw std::runtime_error(
        "jer_closure.scale_factor_files and cones.labels have different "
        "lengths");
  if (cfg.jetTreePaths.size() != cfg.coneLabels.size())
    throw std::runtime_error(
        "trees.jets and cones.labels have different lengths");
  if (cfg.trigCone.IsNull())
    throw std::runtime_error("trigger.cone is missing");
  if (cfg.ptavgEdges.size() < 2)
    throw std::runtime_error("binning.ptavg_edges needs at least 2 edges");
  for (size_t i = 1; i < cfg.ptavgEdges.size(); i++)
    if (cfg.ptavgEdges[i] <= cfg.ptavgEdges[i - 1])
      throw std::runtime_error(
          "binning.ptavg_edges must be strictly ascending");
  if (cfg.etaModeOutput != "both" && cfg.etaModeOutput != "abseta" &&
      cfg.etaModeOutput != "eta")
    throw std::runtime_error(
        "step3.eta_mode must be \"both\", \"abseta\", or \"eta\"");

  return cfg;
}

const AnalysisConfig &Config() {
  static const AnalysisConfig cfg = LoadAnalysisConfig();
  return cfg;
}

void PrintConfigSummary(const AnalysisConfig &cfg) {
  std::cout << "Config: " << cfg.configPath << "\n";
}
