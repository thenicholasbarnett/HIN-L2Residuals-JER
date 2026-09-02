#ifndef JETSMEARER_H
#define JETSMEARER_H

// USAGE
// Instantiate JetSmearer per jet cone/collection
// Construct with JetResolutionObject, CMS JERC txt file format
// Call .SmearedPt(recoPt, eta, rho, genPt) per jet
//
//   JetSmearer smearer("Resolution_AK4PFchs.txt", "ScaleFactor_AK4PFchs.txt");
//   double smearedPt = smearer.SmearedPt(jet.pt, jet.eta, event.rho, genPt);

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <TFormula.h>
#include <TROOT.h>

// vendored from cmssw, private namespace (JME -> JetSmearerJME)
//   CondFormats/JetMETObjects/interface/JetResolutionObject.h
//   CondFormats/JetMETObjects/src/JetResolutionObject.cc
//   CondFormats/JetMETObjects/interface/Utilities.h
//   JetMETCorrections/Modules/interface/JetResolution.h
//   JetMETCorrections/Modules/src/JetResolution.cc

// BEGIN VENDORED CODE

namespace jer_detail {
inline std::vector<std::string> getTokens(const std::string &fLine) {
  std::vector<std::string> tokens;
  std::string currentToken;
  for (unsigned ipos = 0; ipos < fLine.length(); ++ipos) {
    char c = fLine[ipos];
    if (c == '#') {
      break;               // ignore comments
    } else if (c == ' ') { // flush current token if any
      if (!currentToken.empty()) {
        tokens.push_back(currentToken);
        currentToken.clear();
      }
    } else {
      currentToken += c;
    }
  }
  if (!currentToken.empty()) {
    tokens.push_back(currentToken); // flush end
  }
  return tokens;
}
} // namespace jer_detail

enum class Variation { NOMINAL = 0, DOWN = 1, UP = 2 };

template <typename T> T clip(const T &n, const T &lower, const T &upper) {
  return std::max(lower, std::min(n, upper));
}

namespace JetSmearerJME {

template <typename T, typename U> struct bimap {
  typedef std::unordered_map<T, U> left_type;
  typedef std::unordered_map<U, T> right_type;

  left_type left;
  right_type right;

  bimap(std::initializer_list<typename left_type::value_type> l) {
    for (auto &v : l) {
      left.insert(v);
      right.insert(typename right_type::value_type(v.second, v.first));
    }
  }

  bimap() {}
  bimap(bimap &&rhs) {
    left = std::move(rhs.left);
    right = std::move(rhs.right);
  }
};

enum class Binning {
  JetPt = 0,
  JetEta,
  JetAbsEta,
  JetE,
  JetArea,
  Mu,
  Rho,
  NPV,
};

} // namespace JetSmearerJME

namespace std {
template <> struct hash<JetSmearerJME::Binning> {
  typedef JetSmearerJME::Binning argument_type;
  typedef std::size_t result_type;
  hash<uint8_t> int_hash;
  result_type operator()(argument_type const &s) const {
    return int_hash(static_cast<uint8_t>(s));
  }
};
} // namespace std

namespace JetSmearerJME {

class JetParameters {
public:
  typedef std::unordered_map<Binning, float> value_type;

  JetParameters() = default;
  JetParameters(JetParameters &&rhs) { m_values = std::move(rhs.m_values); }
  JetParameters(std::initializer_list<typename value_type::value_type> init) {
    for (auto &i : init) {
      set(i.first, i.second);
    }
  }

  JetParameters &setJetPt(float pt) {
    m_values[Binning::JetPt] = pt;
    return *this;
  }
  JetParameters &setJetEta(float eta) {
    m_values[Binning::JetEta] = eta;
    m_values[Binning::JetAbsEta] = fabs(eta);
    return *this;
  }
  JetParameters &setJetE(float e) {
    m_values[Binning::JetE] = e;
    return *this;
  }
  JetParameters &setJetArea(float area) {
    m_values[Binning::JetArea] = area;
    return *this;
  }
  JetParameters &setMu(float mu) {
    m_values[Binning::Mu] = mu;
    return *this;
  }
  JetParameters &setRho(float rho) {
    m_values[Binning::Rho] = rho;
    return *this;
  }
  JetParameters &setNPV(float npv) {
    m_values[Binning::NPV] = npv;
    return *this;
  }
  JetParameters &set(const Binning &bin, float value) {
    m_values.emplace(bin, value);
    if (bin == Binning::JetEta) {
      m_values.emplace(Binning::JetAbsEta, fabs(value));
    }
    return *this;
  }
  JetParameters &set(const typename value_type::value_type &value) {
    set(value.first, value.second);
    return *this;
  }

  static const bimap<Binning, std::string> binning_to_string;

  std::vector<float> createVector(const std::vector<Binning> &binning) const {
    std::vector<float> values;
    for (const auto &bin : binning) {
      const auto &it = m_values.find(bin);
      if (it == m_values.cend()) {
        throw std::runtime_error(
            "JER parametrisation depends on '" +
            JetParameters::binning_to_string.left.at(bin) +
            "' but no value for this parameter has been specified. Please "
            "call the appropriate 'set' function of the "
            "JetSmearerJME::JetParameters "
            "object");
      }
      values.push_back(it->second);
    }
    return values;
  }

  std::vector<float>
  createVector(const std::vector<std::string> &binname) const {
    std::vector<float> values;
    for (const auto &name : binname) {
      Binning bi = binning_to_string.right.find(name)->second;
      const auto &it = m_values.find(bi);
      if (it == m_values.cend()) {
        std::cerr << "Bin name " << name << " not found!" << std::endl;
        throw std::runtime_error(
            "JER parametrisation depends on '" +
            JetParameters::binning_to_string.left.at(bi) +
            "' but no value for this parameter has been specified. Please "
            "call the appropriate 'set' function of the "
            "JetSmearerJME::JetParameters "
            "object");
      }
      values.push_back(it->second);
    }
    return values;
  }

private:
  value_type m_values;
};

inline const bimap<Binning, std::string> JetParameters::binning_to_string = {
    {Binning::JetPt, "JetPt"},
    {Binning::JetEta, "JetEta"},
    {Binning::JetAbsEta, "JetAbsEta"},
    {Binning::JetE, "JetE"},
    {Binning::JetArea, "JetArea"},
    {Binning::Mu, "Mu"},
    {Binning::Rho, "Rho"},
    {Binning::NPV, "NPV"}};

class JetResolutionObject {
public:
  struct Range {
    float min;
    float max;
    Range() {}
    Range(float mn, float mx) : min(mn), max(mx) {}
    bool is_inside(float value) const {
      return (value >= min) && (value < max);
    }
  };

  class Definition {
  public:
    Definition() {}
    Definition(const std::string &definition) {
      std::vector<std::string> tokens = jer_detail::getTokens(definition);
      if (tokens.size() < 3) {
        throw std::runtime_error(
            "Definition line needs at least three tokens. Please check file "
            "format.");
      }

      size_t n_bins = std::stoul(tokens[0]);
      if (tokens.size() < (n_bins + 2)) {
        throw std::runtime_error("Invalid file format. Please check.");
      }
      for (size_t i = 0; i < n_bins; i++) {
        m_bins_name.push_back(tokens[i + 1]);
      }

      size_t n_variables = std::stoul(tokens[n_bins + 1]);
      if (tokens.size() < (1 + n_bins + 1 + n_variables + 1)) {
        throw std::runtime_error("Invalid file format. Please check.");
      }
      for (size_t i = 0; i < n_variables; i++) {
        m_variables_name.push_back(tokens[n_bins + 2 + i]);
      }

      m_formula_str = tokens[n_bins + n_variables + 2];
      std::string formula_str_lower = m_formula_str;
      std::transform(formula_str_lower.begin(), formula_str_lower.end(),
                     formula_str_lower.begin(), ::tolower);

      if (formula_str_lower == "none") {
        m_formula_str = "";
        if ((tokens.size() > n_bins + n_variables + 3) &&
            (std::atoi(tokens[n_bins + n_variables + 3].c_str()))) {
          size_t n_parameters = std::stoul(tokens[n_bins + n_variables + 3]);
          if (tokens.size() <
              (1 + n_bins + 1 + n_variables + 1 + 1 + n_parameters)) {
            throw std::runtime_error("Invalid file format. Please check.");
          }
          for (size_t i = 0; i < n_parameters; i++) {
            m_formula_str += tokens[n_bins + n_variables + 4 + i] + " ";
          }
        }
      }
      init();
    }

    const std::vector<std::string> &getBinsName() const { return m_bins_name; }
    const std::vector<Binning> &getBins() const { return m_bins; }
    std::string getBinName(size_t bin) const { return m_bins_name[bin]; }
    size_t nBins() const { return m_bins_name.size(); }
    const std::vector<std::string> &getVariablesName() const {
      return m_variables_name;
    }
    const std::vector<Binning> &getVariables() const { return m_variables; }
    std::string getVariableName(size_t variable) const {
      return m_variables_name[variable];
    }
    size_t nVariables() const { return m_variables_name.size(); }
    const std::vector<std::string> &getParametersName() const {
      return m_parameters_name;
    }
    size_t nParameters() const { return m_parameters_name.size(); }
    std::string getFormulaString() const { return m_formula_str; }
    TFormula const *getFormula() const { return m_formula.get(); }

    void init() {
      if (!m_formula_str.empty()) {
        if (m_formula_str.find(' ') == std::string::npos) {
          // bug: vendored source uses a fixed literal name here
          // fix: per-instance unique name, doesn't change formula evaluation
          static std::atomic<unsigned long> sFormulaCounter{0};
          std::string uniqueName =
              "jet_resolution_formula_" + std::to_string(sFormulaCounter++);
          m_formula = std::make_shared<TFormula>(uniqueName.c_str(),
                                                 m_formula_str.c_str());
          if (gROOT) {
            // detach from ROOT global list, otherwise it still double-frees at teardown
            gROOT->GetListOfFunctions()->Remove(m_formula.get());
          }
          // end fix
        } else {
          m_parameters_name = jer_detail::getTokens(m_formula_str);
        }
      }
      for (const auto &bin : m_bins_name) {
        const auto &b = JetParameters::binning_to_string.right.find(bin);
        if (b == JetParameters::binning_to_string.right.cend()) {
          throw std::runtime_error("Bin name not supported: '" + bin + "'");
        }
        m_bins.push_back(b->second);
      }
      for (const auto &v : m_variables_name) {
        const auto &var = JetParameters::binning_to_string.right.find(v);
        if (var == JetParameters::binning_to_string.right.cend()) {
          throw std::runtime_error("Variable name not supported: '" + v + "'");
        }
        m_variables.push_back(var->second);
      }
    }

  private:
    std::vector<std::string> m_bins_name;
    std::vector<std::string> m_variables_name;
    std::string m_formula_str;
    std::shared_ptr<TFormula> m_formula;
    std::vector<Binning> m_bins;
    std::vector<Binning> m_variables;
    std::vector<std::string> m_parameters_name;
  };

  class Record {
  public:
    Record() {}
    Record(const std::string &line, const Definition &def) {
      std::vector<std::string> tokens = jer_detail::getTokens(line);
      if (tokens.size() < (def.nBins() * 2 + def.nVariables() * 2 + 1)) {
        throw std::runtime_error(
            "Invalid record. Please check file format. Record: " + line);
      }

      size_t pos = 0;
      for (size_t i = 0; i < def.nBins(); i++) {
        Range r(std::stof(tokens[pos]), std::stof(tokens[pos + 1]));
        pos += 2;
        m_bins_range.push_back(r);
      }

      size_t n_parameters = std::stoul(tokens[pos++]);
      if (tokens.size() < (def.nBins() * 2 + def.nVariables() * 2 + 1 +
                           (n_parameters - def.nVariables() * 2))) {
        throw std::runtime_error(
            "Invalid record. Please check file format. Record: " + line);
      }

      for (size_t i = 0; i < def.nVariables(); i++) {
        Range r(std::stof(tokens[pos]), std::stof(tokens[pos + 1]));
        pos += 2;
        m_variables_range.push_back(r);
        n_parameters -= 2;
      }
      for (size_t i = 0; i < n_parameters; i++) {
        m_parameters_values.push_back(std::stof(tokens[pos++]));
      }
    }

    const std::vector<Range> &getBinsRange() const { return m_bins_range; }
    const std::vector<Range> &getVariablesRange() const {
      return m_variables_range;
    }
    const std::vector<float> &getParametersValues() const {
      return m_parameters_values;
    }
    size_t nVariables() const { return m_variables_range.size(); }
    size_t nParameters() const { return m_parameters_values.size(); }

  private:
    std::vector<Range> m_bins_range;
    std::vector<Range> m_variables_range;
    std::vector<float> m_parameters_values;
  };

  JetResolutionObject(const std::string &filename) {
    std::ifstream f(filename);
    if (!f.good()) {
      throw std::runtime_error("Can't read input file '" + filename + "'");
    }

    for (std::string line; std::getline(f, line);) {
      if (line.empty() || line[0] == '#') {
        continue;
      }

      size_t first = line.find('{');
      size_t last = line.find('}');
      std::string definition =
          (first != std::string::npos && last != std::string::npos &&
           first < last)
              ? std::string(line, first + 1, last - first - 1)
              : "";

      if (!definition.empty()) {
        m_definition = Definition(definition);
      } else {
        m_records.push_back(Record(line, m_definition));
      }
    }
    m_valid = true;
  }

  JetResolutionObject(const JetResolutionObject &object) {
    m_definition = object.m_definition;
    m_records = object.m_records;
    m_valid = object.m_valid;
    m_definition.init();
  }

  JetResolutionObject() {}

  void dump() const {
    std::cout << "Definition: " << std::endl;
    std::cout << "    Number of binning variables: " << m_definition.nBins()
              << std::endl;
    std::cout << "        ";
    for (const auto &bin : m_definition.getBinsName()) {
      std::cout << bin << ", ";
    }
    std::cout << std::endl;
    std::cout << "    Number of variables: " << m_definition.nVariables()
              << std::endl;
    std::cout << "        ";
    for (const auto &bin : m_definition.getVariablesName()) {
      std::cout << bin << ", ";
    }
    std::cout << std::endl;
    std::cout << "    Formula: " << m_definition.getFormulaString()
              << std::endl;
    std::cout << std::endl << "Bin contents" << std::endl;
    for (const auto &record : m_records) {
      std::cout << "    Bins" << std::endl;
      size_t index = 0;
      for (const auto &bin : record.getBinsRange()) {
        std::cout << "        " << m_definition.getBinName(index) << " ["
                  << bin.min << " - " << bin.max << "]" << std::endl;
        index++;
      }
      std::cout << "    Variables" << std::endl;
      index = 0;
      for (const auto &r : record.getVariablesRange()) {
        std::cout << "        " << m_definition.getVariableName(index) << " ["
                  << r.min << " - " << r.max << "] " << std::endl;
        index++;
      }
      std::cout << "    Parameters" << std::endl;
      index = 0;
      for (const auto &par : record.getParametersValues()) {
        std::cout << "        Parameter #" << index << " = " << par
                  << std::endl;
        index++;
      }
    }
  }

  void saveToFile(const std::string &file) const {
    std::ofstream fout(file);
    fout.setf(std::ios::right);
    fout << "{" << m_definition.nBins();
    for (auto &bin : m_definition.getBinsName()) {
      fout << "    " << bin;
    }
    fout << "    " << m_definition.nVariables();
    for (auto &var : m_definition.getVariablesName()) {
      fout << "    " << var;
    }
    fout << "    "
         << (m_definition.getFormulaString().empty()
                 ? "None"
                 : m_definition.getFormulaString())
         << "    Resolution}" << std::endl;
    for (auto &record : m_records) {
      for (auto &r : record.getBinsRange()) {
        fout << std::left << std::setw(15) << r.min << std::setw(15) << r.max
             << std::setw(15);
      }
      fout << (record.nVariables() * 2 + record.nParameters()) << std::setw(15);
      for (auto &r : record.getVariablesRange()) {
        fout << r.min << std::setw(15) << r.max << std::setw(15);
      }
      for (auto &p : record.getParametersValues()) {
        fout << p << std::setw(15);
      }
      fout << std::endl << std::setw(0);
    }
  }

  const Record *getRecord(const JetParameters &bins_parameters) const {
    if (!m_valid) {
      return nullptr;
    }
    std::vector<float> bins =
        bins_parameters.createVector(m_definition.getBinsName());
    const Record *good_record = nullptr;
    for (const auto &record : m_records) {
      size_t valid_bins = 0;
      size_t current_bin = 0;
      for (const auto &bin : record.getBinsRange()) {
        if (bin.is_inside(bins[current_bin])) {
          valid_bins++;
        }
        current_bin++;
      }
      if (valid_bins == m_definition.nBins()) {
        good_record = &record;
        break;
      }
    }
    return good_record;
  }

  float evaluateFormula(const Record &record,
                        const JetParameters &variables_parameters) const {
    if (!m_valid) {
      return 1;
    }
    auto const *pFormula = m_definition.getFormula();
    if (!pFormula) {
      return 1;
    }
    auto formula = *pFormula;

    std::vector<float> variables =
        variables_parameters.createVector(m_definition.getVariablesName());
    double variables_[4] = {0};
    for (size_t index = 0; index < m_definition.nVariables(); index++) {
      variables_[index] =
          clip(variables[index], record.getVariablesRange()[index].min,
               record.getVariablesRange()[index].max);
    }

    const std::vector<float> &parameters = record.getParametersValues();
    for (size_t index = 0; index < parameters.size(); index++) {
      formula.SetParameter(index, parameters[index]);
    }
    return formula.EvalPar(variables_);
  }

  const std::vector<Record> &getRecords() const { return m_records; }
  const Definition &getDefinition() const { return m_definition; }

private:
  Definition m_definition;
  std::vector<Record> m_records;
  bool m_valid = false;
};

class JetResolution {
public:
  JetResolution(const std::string &filename) {
    m_object = std::make_shared<JetResolutionObject>(filename);
  }
  JetResolution(const JetResolutionObject &object) {
    m_object = std::make_shared<JetResolutionObject>(object);
  }
  JetResolution() {}

  float getResolution(const JetParameters &parameters) const {
    const JetResolutionObject::Record *record = m_object->getRecord(parameters);
    if (!record) {
      return 1;
    }
    return m_object->evaluateFormula(*record, parameters);
  }

  void dump() const { m_object->dump(); }
  const JetResolutionObject *getResolutionObject() const {
    return m_object.get();
  }

private:
  std::shared_ptr<JetResolutionObject> m_object;
};

class JetResolutionScaleFactor {
public:
  JetResolutionScaleFactor(const std::string &filename) {
    m_object = std::make_shared<JetResolutionObject>(filename);
  }
  JetResolutionScaleFactor(const JetResolutionObject &object) {
    m_object = std::make_shared<JetResolutionObject>(object);
  }
  JetResolutionScaleFactor() {}

  float getScaleFactor(const JetParameters &parameters,
                       Variation variation = Variation::NOMINAL,
                       std::string uncertaintySource = "") const {
    const JetResolutionObject::Record *record = m_object->getRecord(parameters);
    if (!record) {
      return 1;
    }

    const std::vector<float> &parameters_values = record->getParametersValues();
    const std::vector<std::string> &parameter_names =
        m_object->getDefinition().getParametersName();
    size_t parameter = static_cast<size_t>(variation);
    if (!uncertaintySource.empty()) {
      if (variation == Variation::DOWN) {
        parameter = std::distance(parameter_names.begin(),
                                  std::find(parameter_names.begin(),
                                            parameter_names.end(),
                                            uncertaintySource + "Down"));
      } else if (variation == Variation::UP) {
        parameter = std::distance(parameter_names.begin(),
                                  std::find(parameter_names.begin(),
                                            parameter_names.end(),
                                            uncertaintySource + "Up"));
      }
      if (parameter >= parameter_names.size()) {
        std::string s;
        for (const auto &piece : parameter_names) {
          s += piece + " ";
        }
        throw std::runtime_error(
            "Invalid value for 'uncertaintySource' parameter. Only " + s +
            " are supported.\n");
      }
    }
    return parameters_values[parameter];
  }

  void dump() const { m_object->dump(); }
  const JetResolutionObject *getResolutionObject() const {
    return m_object.get();
  }

private:
  std::shared_ptr<JetResolutionObject> m_object;
};

} // namespace JetSmearerJME

// END VENDORED CODE

// below implements same algorithm as
// cms-sw/cmssw/tree/master/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h

namespace JetSmearing {

struct Result {
  double smearFactor = 1.0;
  double resolution = 0.0; // sigma_JER used (getResolution() output)
  double scaleFactor = 1.0;
  bool matched = false;
};

// genPt < 0 is no matched gen jet
// this function requires user provides existing gen matched jet pt
inline Result
ComputeSmearFactor(double recoPt, double eta, double rho, double genPt,
                   const JetSmearerJME::JetResolution &resolution,
                   const JetSmearerJME::JetResolutionScaleFactor &resolutionSF,
                   std::mt19937 &rng, Variation variation = Variation::NOMINAL,
                   const std::string &uncertaintySource = "",
                   double dPtMaxFactor = 3.0) {
  Result r;
  r.resolution = resolution.getResolution(
      JetSmearerJME::JetParameters().setJetPt(recoPt).setJetEta(eta).setRho(
          rho));
  r.scaleFactor = resolutionSF.getScaleFactor(
      JetSmearerJME::JetParameters().setJetPt(recoPt).setJetEta(eta), variation,
      uncertaintySource);

  if (genPt >= 0 &&
      std::abs(recoPt - genPt) < dPtMaxFactor * r.resolution * recoPt) {
    // scaling method
    r.matched = true;
    r.smearFactor = 1.0 + (r.scaleFactor - 1.0) * (recoPt - genPt) / recoPt;
  } else if (r.scaleFactor > 1.0) {
    // stochastic method
    double sigma =
        r.resolution * std::sqrt(r.scaleFactor * r.scaleFactor - 1.0);
    std::normal_distribution<double> d(0.0, sigma);
    r.smearFactor = 1.0 + d(rng);
  }
  return r;
}

// pT floor mirrors SmearedJetProducerT's MIN_JET_ENERGY to avoid negative/flipped jet(s)
inline double SmearedPt(double recoPt, double smearFactor,
                        double minPt = 1e-2) {
  double smeared = recoPt * smearFactor;
  return (smeared < minPt) ? minPt : smeared;
}

} // namespace JetSmearing

class JetSmearer {
public:
  // Mirrors SmearedJetProducerT.h default seed
  static constexpr std::uint32_t kDefaultSeed = 37428479;

  JetSmearer(const std::string &resolutionFile,
             const std::string &scaleFactorFile,
             std::uint32_t seed = kDefaultSeed)
      : resolution_(resolutionFile), scaleFactor_(scaleFactorFile), rng_(seed) {
  }

  // smear factor, resolution/scale factor used, scaling or stochastic
  // see JetSmearing::Result above
  JetSmearing::Result Smear(double recoPt, double eta, double rho, double genPt,
                            Variation variation = Variation::NOMINAL,
                            const std::string &uncertaintySource = "",
                            double dPtMaxFactor = 3.0) {
    return JetSmearing::ComputeSmearFactor(recoPt, eta, rho, genPt, resolution_,
                                           scaleFactor_, rng_, variation,
                                           uncertaintySource, dPtMaxFactor);
  }

  // in: reco jet {pT, eta, rho}, gen matched jet pT
  // out: smeared reco jet pT
  double SmearedPt(double recoPt, double eta, double rho, double genPt,
                   Variation variation = Variation::NOMINAL,
                   const std::string &uncertaintySource = "",
                   double dPtMaxFactor = 3.0) {
    JetSmearing::Result r = Smear(recoPt, eta, rho, genPt, variation,
                                  uncertaintySource, dPtMaxFactor);
    return JetSmearing::SmearedPt(recoPt, r.smearFactor);
  }

private:
  JetSmearerJME::JetResolution resolution_;
  JetSmearerJME::JetResolutionScaleFactor scaleFactor_;
  std::mt19937 rng_;
};

#endif
