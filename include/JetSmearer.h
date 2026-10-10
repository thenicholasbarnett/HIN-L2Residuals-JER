#ifndef JETSMEARER_H
#define JETSMEARER_H

// JetSmearer v1.0
// Smear the width of jet energy responses
// Match two jet collections
// Sort a jet collection by pT
// Authored by Nicholas Shawn Barnett

// USAGE
//
// Instantiate JetSmearer object for each jet collection being smeared
//
//   JetSmearer smearer("Resolution_AK4PFchs.txt", "ScaleFactor_AK4PFchs.txt");
//
// Output smeared pT and up/down variations
//
//   double smearedPt = smearer.SmearedPt(jet.pt, jet.eta, event.rho, genPt, event.evt);
//   double smearedPtUp = smearer.SmearedPt(jet.pt, jet.eta, event.rho, genPt, event.evt, Variation::UP);
//   double smearedPtDown = smearer.SmearedPt(jet.pt, jet.eta, event.rho, genPt, event.evt, Variation::DOWN);
//
// or with setters and getters
//
//   smearer.SetJetPT(jet.pt);
//   smearer.SetJetEta(jet.eta);
//   smearer.SetRho(event.rho);
//   smearer.SetGenPT(genPt);
//   smearer.SetEventID(event.evt);
//   double smearedPt = smearer.GetSmearedPT();
//   double smearedPtUp = smearer.GetSmearedPT(Variation::UP);
//   double smearedPtDown = smearer.GetSmearedPT(Variation::DOWN);
//
// Scaling: 1 + (SF - 1)(pT - pT_gen)/pT
// Stochastic: 1 + sqrt(max(SF^2 - 1, 0)) sigma_JER N(0,1)
// N seeded from a hash of (pT, eta, rho, eventID)
// Hybrid: Scaling for matched jets within |pT - pT_gen| < 3 sigma_JER pT, else stochastic (Default & JME recommendation)
//
// smearer.Match returns array of indices with nref entries
// each entry maps first-collection jet index to matched second-collection jet index 
// entry is -999 when no match
//
// Nearest: closest pair of jets, each jet can be shared
// OneToOne: closest pairs first, no jet matched twice
// JME: closest pair, dR < R/2, and |pT - pT_gen| < 3 sigma_JER pT
//
//   std::vector<int> match = smearer.Match(nref, ptCorr, jteta, jtphi,
//                                          ngen, genpt, geneta, genphi,
//                                          0.4, rho, JetSmearing::Mode::JME);
// match[i] is now index of gen jet matched to reco jet i
//
// Reproduce HiForest ref jet collection with dRFraction=1.0 and OneToOne mode
// JetSmearing::Match(nref, jteta, jtphi, ngen, geneta, genphi, R, mode, dRFraction);
//
// JetSmearing::Order(n, pt) returns n long array,
// each entry indexes an entry in pt (input array) by descending value
//
//   std::vector<int> order = JetSmearing::Order(nref, ptSmeared);
// match[order[0]] is now leading reco jet's matched gen jet index

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
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

#ifndef JET_VARIATION_ENUM
#define JET_VARIATION_ENUM
enum class Variation { NOMINAL = 0, DOWN = 1, UP = 2 };
#endif

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
            // detach from ROOT global list, otherwise still double-frees at teardown
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

// below is JME's correctionlib JERSmear (jer_smear.json) and the hybrid method of
// cms-sw/cmssw/tree/master/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h

namespace JetSmearing {

enum class Method { Hybrid, JME, Scaling, Stochastic };

// "hybrid" | "jme" | "scaling" | "stochastic", throws otherwise
inline Method MethodFromString(const std::string &name) {
  if (name == "hybrid") {
    return Method::Hybrid;
  }
  if (name == "jme") {
    return Method::JME;
  }
  if (name == "scaling") {
    return Method::Scaling;
  }
  if (name == "stochastic") {
    return Method::Stochastic;
  }
  throw std::invalid_argument(
      "JetSmearing: unknown method \"" + name +
      "\", expected hybrid, jme, scaling or stochastic");
}

// correctionlib's hashprng, distribution "normal": XXH64 of the inputs' 64-bit
// patterns (seed 0) seeds pcg32_oneseq, then a Marsaglia polar draw
namespace detail {

constexpr std::uint64_t kP1 = 11400714785074694791ULL;
constexpr std::uint64_t kP2 = 14029467366897019727ULL;
constexpr std::uint64_t kP3 = 1609587929392839161ULL;
constexpr std::uint64_t kP4 = 9650029242287828579ULL;
constexpr std::uint64_t kP5 = 2870177450012600261ULL;

inline std::uint64_t Rotl64(std::uint64_t x, int r) {
  return (x << r) | (x >> (64 - r));
}

inline std::uint64_t Read64(const unsigned char *p) {
  std::uint64_t v;
  std::memcpy(&v, p, 8);
  return v;
}

inline std::uint64_t Round64(std::uint64_t acc, std::uint64_t in) {
  acc += in * kP2;
  return Rotl64(acc, 31) * kP1;
}

inline std::uint64_t Merge64(std::uint64_t acc, std::uint64_t val) {
  acc ^= Round64(0, val);
  return acc * kP1 + kP4;
}

// XXH64 (xxHash, Yann Collet, BSD 2-clause), little-endian
inline std::uint64_t XXH64(const void *input, std::size_t len,
                           std::uint64_t seed) {
  const unsigned char *p = (const unsigned char *)input;
  const unsigned char *end = p + len;
  std::uint64_t h;
  if (len >= 32) {
    const unsigned char *limit = end - 32;
    std::uint64_t v1 = seed + kP1 + kP2, v2 = seed + kP2, v3 = seed,
                  v4 = seed - kP1;
    do {
      v1 = Round64(v1, Read64(p));
      v2 = Round64(v2, Read64(p + 8));
      v3 = Round64(v3, Read64(p + 16));
      v4 = Round64(v4, Read64(p + 24));
      p += 32;
    } while (p <= limit);
    h = Rotl64(v1, 1) + Rotl64(v2, 7) + Rotl64(v3, 12) + Rotl64(v4, 18);
    h = Merge64(h, v1);
    h = Merge64(h, v2);
    h = Merge64(h, v3);
    h = Merge64(h, v4);
  } else {
    h = seed + kP5;
  }
  h += (std::uint64_t)len;
  while (p + 8 <= end) {
    h ^= Round64(0, Read64(p));
    h = Rotl64(h, 27) * kP1 + kP4;
    p += 8;
  }
  if (p + 4 <= end) {
    std::uint32_t v;
    std::memcpy(&v, p, 4);
    h ^= (std::uint64_t)v * kP1;
    h = Rotl64(h, 23) * kP2 + kP3;
    p += 4;
  }
  while (p < end) {
    h ^= (*p) * kP5;
    h = Rotl64(h, 11) * kP1;
    p++;
  }
  h ^= h >> 33;
  h *= kP2;
  h ^= h >> 29;
  h *= kP3;
  h ^= h >> 32;
  return h;
}

// pcg32_oneseq (pcg-cpp): 64-bit LCG state, XSH RR output of the old state
class Pcg32 {
public:
  explicit Pcg32(std::uint64_t seed) : state_(Bump(seed + kInc)) {}
  std::uint32_t operator()() {
    const std::uint64_t old = state_;
    state_ = Bump(state_);
    const std::uint32_t xsh = (std::uint32_t)(((old >> 18) ^ old) >> 27);
    const unsigned rot = (unsigned)(old >> 59);
    return (xsh >> rot) | (xsh << ((32 - rot) & 31));
  }

private:
  static constexpr std::uint64_t kMult = 6364136223846793005ULL;
  static constexpr std::uint64_t kInc = 1442695040888963407ULL;
  static std::uint64_t Bump(std::uint64_t s) { return s * kMult + kInc; }
  std::uint64_t state_;
};

// rounds a product before it's added, so no FMA fusing (bit-identity).
inline double Rounded(double x) {
  volatile double r = x;
  return r;
}

} // namespace detail

// N(0,1) from (pT, eta, rho, eventID), bit for bit as JERSmear's hashprng
// verified against correctionlib 2.9.0 on jer_smear.json, 200k jets
inline double HashNormal(double pt, double eta, double rho,
                         std::int64_t eventID) {
  std::uint64_t data[4];
  std::memcpy(&data[0], &pt, 8);
  std::memcpy(&data[1], &eta, 8);
  std::memcpy(&data[2], &rho, 8);
  data[3] = (std::uint64_t)eventID;
  detail::Pcg32 gen(detail::XXH64(data, sizeof(data), 0));
  double u, v, s;
  do {
    u = std::ldexp((double)gen(), -31) - 1;
    v = std::ldexp((double)gen(), -31) - 1;
    s = std::fma(u, u, detail::Rounded(v * v));
  } while (s >= 1.0 || s == 0.0);
  return u * std::sqrt(-2.0 * std::log(s) / s);
}

// JERSmear: genPt >= 0 scaling, genPt < 0 stochastic
inline double JERSmear(double pt, double eta, double genPt, double rho,
                       std::int64_t eventID, double jer, double jersf) {
  if (genPt >= 0) {
    return 1 + (jersf - 1) * (pt - genPt) / pt;
  }
  const double n = HashNormal(pt, eta, rho, eventID);
  const double m = detail::Rounded(jersf * jersf) - 1;
  const double t = std::sqrt(std::max(m, 0.0)) * jer;
  return 1 + detail::Rounded(t * n);
}

struct Result {
  double smearFactor = 1.0;
  double resolution = 0.0; // sigma_JER used (getResolution() output)
  double scaleFactor = 1.0;
  bool matched = false; // scaling method applied

  void Print() const {
    printf("smearFactor  %.6f\n", smearFactor);
    printf("sigma_JER    %.6f\n", resolution);
    printf("JER SF       %.6f\n", scaleFactor);
    printf("branch       %s\n",
           matched ? "scaling (gen matched)" : "stochastic, or left alone");
  }
};

// smear factor from sigma_JER and SF already looked up
inline Result ComputeSmearFactor(double recoPt, double eta, double rho,
                                 double genPt, std::int64_t eventID,
                                 double resolution, double scaleFactor,
                                 double dPtMaxFactor = 3.0,
                                 Method method = Method::Hybrid) {
  Result r;
  r.resolution = resolution;
  r.scaleFactor = scaleFactor;

  const bool wellMatched =
      genPt >= 0 &&
      std::abs(recoPt - genPt) < dPtMaxFactor * r.resolution * recoPt;
  double useGen = -1; // stochastic
  if (method == Method::JME) {
    useGen = genPt;
  } else if (method != Method::Stochastic && wellMatched) {
    useGen = genPt;
  } else if (method == Method::Scaling) {
    return r; // unmatched, left alone
  }
  r.matched = useGen >= 0;
  r.smearFactor =
      JERSmear(recoPt, eta, useGen, rho, eventID, r.resolution, r.scaleFactor);
  return r;
}

// same, looking sigma_JER and SF up in the text-file readers
inline Result
ComputeSmearFactor(double recoPt, double eta, double rho, double genPt,
                   std::int64_t eventID,
                   const JetSmearerJME::JetResolution &resolution,
                   const JetSmearerJME::JetResolutionScaleFactor &resolutionSF,
                   Variation variation = Variation::NOMINAL,
                   const std::string &uncertaintySource = "",
                   double dPtMaxFactor = 3.0, Method method = Method::Hybrid) {
  const double resolutionValue = resolution.getResolution(
      JetSmearerJME::JetParameters().setJetPt(recoPt).setJetEta(eta).setRho(
          rho));
  const double sf = resolutionSF.getScaleFactor(
      JetSmearerJME::JetParameters().setJetPt(recoPt).setJetEta(eta), variation,
      uncertaintySource);
  return ComputeSmearFactor(recoPt, eta, rho, genPt, eventID, resolutionValue,
                            sf, dPtMaxFactor, method);
}

inline double SmearedPt(double recoPt, double smearFactor) {
  return recoPt * smearFactor;
}

// unmatched entry
constexpr int kUnmatched = -999;

enum class Mode { OneToOne, Nearest, JME };

// JME matching conditions
constexpr double kJMEdRFraction = 0.5; // dR < R/2
constexpr double kJMENSigma = 3.0;     // |pT - pT_gen| < 3 sigma_JER pT

template <typename T> std::vector<int> Order(int n, const T *pt) {
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(),
                   [&](int a, int b) { return pt[a] > pt[b]; });
  return order;
}

inline double DeltaR(double eta1, double phi1, double eta2, double phi2) {
  const double dphi = std::remainder(phi1 - phi2, 2 * M_PI);
  const double deta = eta1 - eta2;
  return std::sqrt(deta * deta + dphi * dphi);
}

namespace detail {

// each reco jet matches closest gen jet within maxDR
// checks JME pT window if jtpt is provided
template <typename T>
std::vector<int> MatchNearest(int nref, const T *jteta, const T *jtphi,
                              int ngen, const T *geneta, const T *genphi,
                              double maxDR, const T *jtpt = nullptr,
                              const T *genpt = nullptr,
                              const double *sigmaJER = nullptr) {
  std::vector<int> match(nref, kUnmatched);
  for (int i = 0; i < nref; i++) {
    double best = maxDR;
    for (int g = 0; g < ngen; g++) {
      if (jtpt &&
          std::abs(jtpt[i] - genpt[g]) >= kJMENSigma * sigmaJER[i] * jtpt[i]) {
        continue;
      }
      const double dr = DeltaR(jteta[i], jtphi[i], geneta[g], genphi[g]);
      if (dr < best) {
        best = dr;
        match[i] = g;
      }
    }
  }
  return match;
}

// every pair within maxDR, smallest dR first, no jet used twice
template <typename T>
std::vector<int> MatchOneToOne(int nref, const T *jteta, const T *jtphi,
                               int ngen, const T *geneta, const T *genphi,
                               double maxDR) {
  std::vector<int> match(nref, kUnmatched);
  std::vector<std::tuple<double, int, int>> pairs;
  for (int i = 0; i < nref; i++) {
    for (int g = 0; g < ngen; g++) {
      const double dr = DeltaR(jteta[i], jtphi[i], geneta[g], genphi[g]);
      if (dr < maxDR) {
        pairs.emplace_back(dr, i, g);
      }
    }
  }
  std::sort(pairs.begin(), pairs.end());
  std::vector<bool> genUsed(ngen, false);
  for (const auto &p : pairs) {
    const int i = std::get<1>(p);
    const int g = std::get<2>(p);
    if (match[i] < 0 && !genUsed[g]) {
      match[i] = g;
      genUsed[g] = true;
    }
  }
  return match;
}

} // namespace detail

template <typename T>
std::vector<int> Match(int nref, const T *jteta, const T *jtphi, int ngen,
                       const T *geneta, const T *genphi, double coneR,
                       Mode mode = Mode::OneToOne, double dRFraction = 0.5) {
  if (mode == Mode::JME) {
    throw std::invalid_argument("JetSmearing::Match: Mode::JME needs "
                                "sigma_JER, use JetSmearer::Match");
  }
  if (mode == Mode::Nearest) {
    return detail::MatchNearest(nref, jteta, jtphi, ngen, geneta, genphi,
                                dRFraction * coneR);
  }
  return detail::MatchOneToOne(nref, jteta, jtphi, ngen, geneta, genphi,
                               dRFraction * coneR);
}

// JER SF in JME's 2024+ format: a JEC-style formula file, SF(pT) =
// sigma_data / sigma_MC per eta bin -- {1 JetEta 1 JetPt <formula> Correction
// ...}, each record: eta range, N, pT range, parameters; pT clamped to the
// record's range, evaluated in double. 1 outside every eta bin
class FormulaScaleFactor {
public:
  explicit FormulaScaleFactor(const std::string &file) {
    std::ifstream in(file);
    if (!in) {
      throw std::runtime_error("JetSmearer: can't open " + file);
    }
    std::string line;
    bool header = false;
    while (std::getline(in, line)) {
      std::vector<std::string> t = Tokens(line);
      if (t.empty()) {
        continue;
      }
      if (!header) {
        // {nBin binVar nDep depVar formula Correction label}
        if (t.size() < 5 || std::stoi(t[0]) != 1 || std::stoi(t[2]) != 1 ||
            (t[1] != "JetEta" && t[1] != "JetAbsEta") || t[3] != "JetPt") {
          throw std::runtime_error("JetSmearer: " + file +
                                   ": expected {1 JetEta 1 JetPt <formula> "
                                   "Correction ...}");
        }
        absEta_ = t[1] == "JetAbsEta";
        formula_ = std::make_shared<TFormula>(
            ("JetSmearerSF_" + std::to_string(Counter()++)).c_str(),
            t[4].c_str());
        if (gROOT) {
          gROOT->GetListOfFunctions()->Remove(formula_.get());
        }
        header = true;
        continue;
      }
      Record r;
      r.etaLo = std::stod(t[0]);
      r.etaHi = std::stod(t[1]);
      const int n = std::stoi(t[2]);
      if (n < 2 || (int)t.size() < 3 + n) {
        throw std::runtime_error("JetSmearer: " + file + ": bad record");
      }
      r.ptLo = std::stod(t[3]);
      r.ptHi = std::stod(t[4]);
      for (int i = 2; i < n; i++) {
        r.par.push_back(std::stod(t[3 + i]));
      }
      records_.push_back(r);
    }
    if (!header || records_.empty()) {
      throw std::runtime_error("JetSmearer: " + file + ": no records");
    }
  }

  double Evaluate(double pt, double eta) const {
    const double e = absEta_ ? std::abs(eta) : eta;
    for (const Record &r : records_) {
      if (e >= r.etaLo && e < r.etaHi) {
        formula_->SetParameters(r.par.data());
        return formula_->Eval(std::min(std::max(pt, r.ptLo), r.ptHi));
      }
    }
    return 1;
  }

  // "{... <formula> Correction ...}" rather than "{... ScaleFactor}"
  static bool IsFormulaFile(const std::string &file) {
    std::ifstream in(file);
    std::string line;
    while (std::getline(in, line)) {
      if (line.find_first_not_of(" \t\r") == std::string::npos) {
        continue;
      }
      return line.find("Correction") != std::string::npos &&
             line.find("ScaleFactor") == std::string::npos;
    }
    return false;
  }

  // whitespace tokens, braces dropped
  static std::vector<std::string> Tokens(std::string line) {
    for (char &c : line) {
      if (c == '{' || c == '}') {
        c = ' ';
      }
    }
    std::istringstream ss(line);
    std::vector<std::string> t;
    std::string w;
    while (ss >> w) {
      t.push_back(w);
    }
    return t;
  }

private:
  struct Record {
    double etaLo, etaHi, ptLo, ptHi;
    std::vector<double> par;
  };
  std::vector<Record> records_;
  std::shared_ptr<TFormula> formula_;
  bool absEta_ = false;

  static int &Counter() {
    static int n = 0;
    return n;
  }
};

// JER SF uncertainty in JME's 2024+ format: a JEC-style uncertainty file,
// {1 JetEta 1 JetPt "" Correction Uncertainty}, each record: eta range, N,
// then (pT, up, down) points. As JME's correctionlib JSON: the first value is
// the uncertainty, linear in pT between points, clamped at the ends, applied
// symmetric (SF * (1 +- unc)). 0 outside every eta bin
class ScaleFactorUncertainty {
public:
  explicit ScaleFactorUncertainty(const std::string &file) {
    std::ifstream in(file);
    if (!in) {
      throw std::runtime_error("JetSmearer: can't open " + file);
    }
    std::string line;
    bool header = false;
    while (std::getline(in, line)) {
      std::vector<std::string> t = FormulaScaleFactor::Tokens(line);
      if (t.empty()) {
        continue;
      }
      if (!header) {
        if (line.find("Uncertainty") == std::string::npos || t.size() < 4 ||
            t[1] != "JetEta" || t[3] != "JetPt") {
          throw std::runtime_error("JetSmearer: " + file +
                                   ": expected {1 JetEta 1 JetPt \"\" "
                                   "Correction Uncertainty}");
        }
        header = true;
        continue;
      }
      Record r;
      r.etaLo = std::stod(t[0]);
      r.etaHi = std::stod(t[1]);
      const int n = std::stoi(t[2]);
      for (int i = 0; i + 2 < n && 5 + i < (int)t.size(); i += 3) {
        r.pt.push_back(std::stod(t[3 + i]));
        r.unc.push_back(std::stod(t[4 + i]));
      }
      if (r.pt.empty()) {
        throw std::runtime_error("JetSmearer: " + file + ": bad record");
      }
      records_.push_back(r);
    }
    if (!header || records_.empty()) {
      throw std::runtime_error("JetSmearer: " + file + ": no records");
    }
  }

  double Evaluate(double pt, double eta) const {
    for (const Record &r : records_) {
      if (eta < r.etaLo || eta >= r.etaHi) {
        continue;
      }
      if (pt <= r.pt.front()) {
        return r.unc.front();
      }
      if (pt >= r.pt.back()) {
        return r.unc.back();
      }
      for (size_t i = 0; i + 1 < r.pt.size(); i++) {
        if (pt >= r.pt[i] && pt < r.pt[i + 1]) {
          return r.unc[i] + (r.unc[i + 1] - r.unc[i]) /
                                (r.pt[i + 1] - r.pt[i]) * (pt - r.pt[i]);
        }
      }
    }
    return 0;
  }

private:
  struct Record {
    double etaLo, etaHi;
    std::vector<double> pt, unc;
  };
  std::vector<Record> records_;
};

} // namespace JetSmearing

class JetSmearer {
public:
  // SF file: the table format ({... ScaleFactor}, down/up columns included)
  // or JME's 2024+ formula format ({... Correction ...}), told apart by the
  // header; the latter takes its variations from the SF uncertainty file
  JetSmearer(const std::string &resolutionFile,
             const std::string &scaleFactorFile,
             JetSmearing::Method method = JetSmearing::Method::Hybrid)
      : resolution_(resolutionFile), method_(method) {
    LoadScaleFactor(scaleFactorFile);
  }

  JetSmearer(const std::string &resolutionFile,
             const std::string &scaleFactorFile,
             const std::string &scaleFactorUncertaintyFile,
             JetSmearing::Method method = JetSmearing::Method::Hybrid)
      : resolution_(resolutionFile), method_(method) {
    LoadScaleFactor(scaleFactorFile);
    if (!formulaSF_) {
      throw std::runtime_error("JetSmearer: " + scaleFactorFile +
                               " is a table SF file, its down/up are in the "
                               "file; no SF uncertainty file with it");
    }
    sfUncertainty_ = std::make_shared<JetSmearing::ScaleFactorUncertainty>(
        scaleFactorUncertaintyFile);
  }

  void SetMethod(JetSmearing::Method method) { method_ = method; }
  JetSmearing::Method GetMethod() const { return method_; }

  // sigma_JER from the resolution file
  double Resolution(double pt, double eta, double rho) const {
    return resolution_.getResolution(
        JetSmearerJME::JetParameters().setJetPt(pt).setJetEta(eta).setRho(rho));
  }

  // JER SF, at the JEC-corrected pT
  double ScaleFactor(double pt, double eta,
                     Variation variation = Variation::NOMINAL,
                     const std::string &uncertaintySource = "") const {
    if (tableSF_) {
      return tableSF_->getScaleFactor(
          JetSmearerJME::JetParameters().setJetPt(pt).setJetEta(eta), variation,
          uncertaintySource);
    }
    if (!uncertaintySource.empty()) {
      throw std::runtime_error("JetSmearer: uncertainty sources need a table "
                               "SF file");
    }
    const double sf = formulaSF_->Evaluate(pt, eta);
    if (variation == Variation::NOMINAL) {
      return sf;
    }
    if (!sfUncertainty_) {
      throw std::runtime_error("JetSmearer: Variation::UP/DOWN with a formula "
                               "SF file needs the SF uncertainty file");
    }
    const double unc = sfUncertainty_->Evaluate(pt, eta);
    return variation == Variation::UP ? sf * (1 + unc) : sf * (1 - unc);
  }

  // smear factor, resolution/scale factor used
  JetSmearing::Result Smear(double recoPt, double eta, double rho, double genPt,
                            std::int64_t eventID,
                            Variation variation = Variation::NOMINAL,
                            const std::string &uncertaintySource = "",
                            double dPtMaxFactor = 3.0) const {
    if (tableSF_) {
      return JetSmearing::ComputeSmearFactor(
          recoPt, eta, rho, genPt, eventID, resolution_, *tableSF_, variation,
          uncertaintySource, dPtMaxFactor, method_);
    }
    return JetSmearing::ComputeSmearFactor(
        recoPt, eta, rho, genPt, eventID, Resolution(recoPt, eta, rho),
        ScaleFactor(recoPt, eta, variation, uncertaintySource), dPtMaxFactor,
        method_);
  }

  // in: reco jet {pT, eta, rho}, gen matched jet pT, event number
  // out: smeared reco jet pT
  double SmearedPt(double recoPt, double eta, double rho, double genPt,
                   std::int64_t eventID,
                   Variation variation = Variation::NOMINAL,
                   const std::string &uncertaintySource = "",
                   double dPtMaxFactor = 3.0) const {
    JetSmearing::Result r = Smear(recoPt, eta, rho, genPt, eventID, variation,
                                  uncertaintySource, dPtMaxFactor);
    return JetSmearing::SmearedPt(recoPt, r.smearFactor);
  }

  // setters and getters
  // throw if an input was never set
  void SetJetPT(double value) { setPt_ = value; }
  void SetJetEta(double value) { setEta_ = value; }
  void SetRho(double value) { setRho_ = value; }
  void SetGenPT(double value) { setGenPt_ = value; }
  void SetEventID(std::int64_t value) {
    setEventID_ = value;
    hasEventID_ = true;
  }

  JetSmearing::Result GetSmear(Variation variation = Variation::NOMINAL,
                               const std::string &uncertaintySource = "",
                               double dPtMaxFactor = 3.0) const {
    CheckSet(true);
    return Smear(setPt_, setEta_, setRho_, setGenPt_, setEventID_, variation,
                 uncertaintySource, dPtMaxFactor);
  }

  double GetSmearedPT(Variation variation = Variation::NOMINAL,
                      const std::string &uncertaintySource = "",
                      double dPtMaxFactor = 3.0) const {
    CheckSet(true);
    return SmearedPt(setPt_, setEta_, setRho_, setGenPt_, setEventID_,
                     variation, uncertaintySource, dPtMaxFactor);
  }

  double GetResolution() const {
    CheckSet(false);
    return Resolution(setPt_, setEta_, setRho_);
  }

  double GetScaleFactor(Variation variation = Variation::NOMINAL,
                        const std::string &uncertaintySource = "") const {
    CheckSet(false);
    return ScaleFactor(setPt_, setEta_, variation, uncertaintySource);
  }

  template <typename T>
  std::vector<int>
  Match(int nref, const T *jtpt, const T *jteta, const T *jtphi, int ngen,
        const T *genpt, const T *geneta, const T *genphi, double coneR,
        double rho, JetSmearing::Mode mode = JetSmearing::Mode::OneToOne,
        double dRFraction = 0.5) const {
    if (mode != JetSmearing::Mode::JME) {
      return JetSmearing::Match(nref, jteta, jtphi, ngen, geneta, genphi, coneR,
                                mode, dRFraction);
    }
    std::vector<double> sigma(nref);
    for (int i = 0; i < nref; i++) {
      sigma[i] = Resolution(jtpt[i], jteta[i], rho);
    }
    return JetSmearing::detail::MatchNearest(
        nref, jteta, jtphi, ngen, geneta, genphi,
        JetSmearing::kJMEdRFraction * coneR, jtpt, genpt, sigma.data());
  }

private:
  JetSmearerJME::JetResolution resolution_;
  // one of the two SF formats, and the formula format's uncertainty
  std::shared_ptr<JetSmearerJME::JetResolutionScaleFactor> tableSF_;
  std::shared_ptr<JetSmearing::FormulaScaleFactor> formulaSF_;
  std::shared_ptr<JetSmearing::ScaleFactorUncertainty> sfUncertainty_;
  JetSmearing::Method method_;

  void LoadScaleFactor(const std::string &file) {
    if (JetSmearing::FormulaScaleFactor::IsFormulaFile(file)) {
      formulaSF_ = std::make_shared<JetSmearing::FormulaScaleFactor>(file);
    } else {
      tableSF_ =
          std::make_shared<JetSmearerJME::JetResolutionScaleFactor>(file);
    }
  }

  // setter-style inputs, NaN until set
  double setPt_ = std::nan("");
  double setEta_ = std::nan("");
  double setRho_ = std::nan("");
  double setGenPt_ = std::nan("");
  std::int64_t setEventID_ = 0;
  bool hasEventID_ = false;

  void CheckSet(bool smearing) const {
    if (std::isnan(setPt_) || std::isnan(setEta_) || std::isnan(setRho_) ||
        (smearing && (std::isnan(setGenPt_) || !hasEventID_))) {
      throw std::logic_error(
          smearing ? "JetSmearer: SetJetPT, SetJetEta, SetRho, SetGenPT and "
                     "SetEventID before GetSmear / GetSmearedPT"
                   : "JetSmearer: SetJetPT, SetJetEta and SetRho before "
                     "GetResolution");
    }
  }
};

#endif
