#ifndef JETSELECTOR_H
#define JETSELECTOR_H

// JetSelector v1.0
// jet ID + jet veto maps from CMS correctionlib JSON with just this header
// Author: Nicholas Shawn Barnett

// USAGE
// Instantiate JetSelector once, with the collision system
// Construct with jet ID + jet veto map JSON files, CMS correctionlib format
// Call .JetSelection(...) per jet to drop bad jets, or .VetoEvent(...) once
// per event, with the forest arrays, to drop the whole event
//
//   JetSelector js(JetSelector::System::pp, "json/jetid.json",
//                  "json/jetvetomaps_Summer24Prompt24_RunBCDEFGHI_V1.json");
//   if (!js.JetSelection(jteta[i], jtphi[i], jtPfCHF[i], jtPfNHF[i],
//                        jtPfCEF[i], jtPfNEF[i], jtPfMUF[i], jtPfCHM[i],
//                        jtPfNHM[i], jtPfCEM[i], jtPfNEM[i], jtPfMUM[i])) {
//     continue;
//   }
//
//   if (js.VetoEvent(nref, ptCorr, jteta, jtphi, jtPfCHF, jtPfNHF, jtPfCEF,
//                    jtPfNEF, jtPfMUF, jtPfCHM, jtPfNHM, jtPfCEM, jtPfNEM,
//                    jtPfMUM)) {
//     continue;
//   }
//
// Calibration work uses the strictest veto map:
//
//   JetSelector js(JetSelector::System::pp, "json/jetid.json",
//                  "json/jetvetomaps_Summer24Prompt24_RunBCDEFGHI_V1.json",
//                  JetSelector::Purpose::Calibration);

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// minimal JSON reader, enough for correctionlib files
namespace JetSelectorJSON {

struct Value {
  enum class Type { Null, Bool, Number, String, Array, Object };
  Type type = Type::Null;
  bool boolean = false;
  double number = 0.0;
  std::string string;
  std::vector<Value> array;
  std::vector<std::string> keys;
  std::vector<Value> values;

  bool Has(const std::string &key) const {
    return std::find(keys.begin(), keys.end(), key) != keys.end();
  }

  const Value &At(const std::string &key) const {
    auto it = std::find(keys.begin(), keys.end(), key);
    if (it == keys.end()) {
      throw std::runtime_error("JetSelector: missing JSON key '" + key + "'");
    }
    return values[it - keys.begin()];
  }
};

class Parser {
public:
  explicit Parser(const std::string &text) : s_(text) {}

  Value Parse() {
    Value v = ParseValue();
    SkipSpace();
    if (pos_ != s_.size()) {
      Fail("trailing characters");
    }
    return v;
  }

private:
  const std::string &s_;
  size_t pos_ = 0;

  [[noreturn]] void Fail(const std::string &what) const {
    throw std::runtime_error("JetSelector: JSON parse error at byte " +
                             std::to_string(pos_) + ": " + what);
  }

  void SkipSpace() {
    while (pos_ < s_.size() && (s_[pos_] == ' ' || s_[pos_] == '\n' ||
                                s_[pos_] == '\r' || s_[pos_] == '\t')) {
      pos_++;
    }
  }

  void Expect(char c) {
    SkipSpace();
    if (pos_ >= s_.size() || s_[pos_] != c) {
      Fail(std::string("expected '") + c + "'");
    }
    pos_++;
  }

  bool Consume(const char *word) {
    size_t n = std::char_traits<char>::length(word);
    if (s_.compare(pos_, n, word) == 0) {
      pos_ += n;
      return true;
    }
    return false;
  }

  Value ParseValue() {
    SkipSpace();
    if (pos_ >= s_.size()) {
      Fail("unexpected end of input");
    }
    Value v;
    char c = s_[pos_];
    if (c == '{') {
      v.type = Value::Type::Object;
      pos_++;
      SkipSpace();
      if (pos_ < s_.size() && s_[pos_] == '}') {
        pos_++;
        return v;
      }
      while (true) {
        SkipSpace();
        v.keys.push_back(ParseString());
        Expect(':');
        v.values.push_back(ParseValue());
        SkipSpace();
        if (pos_ < s_.size() && s_[pos_] == ',') {
          pos_++;
          continue;
        }
        Expect('}');
        return v;
      }
    }
    if (c == '[') {
      v.type = Value::Type::Array;
      pos_++;
      SkipSpace();
      if (pos_ < s_.size() && s_[pos_] == ']') {
        pos_++;
        return v;
      }
      while (true) {
        v.array.push_back(ParseValue());
        SkipSpace();
        if (pos_ < s_.size() && s_[pos_] == ',') {
          pos_++;
          continue;
        }
        Expect(']');
        return v;
      }
    }
    if (c == '"') {
      v.type = Value::Type::String;
      v.string = ParseString();
      return v;
    }
    if (Consume("true")) {
      v.type = Value::Type::Bool;
      v.boolean = true;
      return v;
    }
    if (Consume("false")) {
      v.type = Value::Type::Bool;
      return v;
    }
    if (Consume("null")) {
      return v;
    }
    const char *start = s_.c_str() + pos_;
    char *end = nullptr;
    v.number = std::strtod(start, &end);
    if (end == start) {
      Fail("unexpected character");
    }
    v.type = Value::Type::Number;
    pos_ += end - start;
    return v;
  }

  std::string ParseString() {
    if (pos_ >= s_.size() || s_[pos_] != '"') {
      Fail("expected string");
    }
    pos_++;
    std::string out;
    while (pos_ < s_.size() && s_[pos_] != '"') {
      char c = s_[pos_++];
      if (c != '\\') {
        out += c;
        continue;
      }
      if (pos_ >= s_.size()) {
        Fail("unterminated escape");
      }
      char e = s_[pos_++];
      switch (e) {
      case '"':
      case '\\':
      case '/':
        out += e;
        break;
      case 'b':
        out += '\b';
        break;
      case 'f':
        out += '\f';
        break;
      case 'n':
        out += '\n';
        break;
      case 'r':
        out += '\r';
        break;
      case 't':
        out += '\t';
        break;
      default:
        Fail("unsupported string escape");
      }
    }
    if (pos_ >= s_.size()) {
      Fail("unterminated string");
    }
    pos_++;
    return out;
  }
};

inline Value ParseFile(const std::string &path) {
  std::ifstream f(path);
  if (!f.is_open()) {
    throw std::runtime_error("JetSelector: could not open " + path);
  }
  std::stringstream ss;
  ss << f.rdbuf();
  std::string text = ss.str();
  return Parser(text).Parse();
}

} // namespace JetSelectorJSON

// evaluates the correctionlib node types jet ID + veto maps use:
// binning, multibinning, category, transform (abs only), plain numbers
namespace JetSelectorCorrection {

struct Input {
  double number = 0.0;
  std::string string;
};

struct Node {
  enum class Kind { Leaf, Binning, MultiBinning, Category, Transform };
  enum class Flow { Clamp, Error, Default };

  Kind kind = Kind::Leaf;
  double value = 0.0;

  // binning (1 input), multibinning (N inputs), transform, category
  std::vector<int> inputs;
  std::vector<std::vector<double>> edges;
  std::vector<Node> content;
  Flow flow = Flow::Error;
  std::vector<Node> fallback; // flow/default node, 0 or 1 entries

  // category
  std::vector<std::string> keys;

  // binning cut ignored, always follows content[passChild]
  bool skip = false;
  int passChild = -1;
};

inline bool IsFail(const Node &n) {
  return n.kind == Node::Kind::Leaf && n.value == 0.0;
}

// ignore every binning cut on these inputs, e.g. hadronic ID criteria
inline void SkipInputs(Node &n, const std::vector<int> &skipped) {
  for (Node &c : n.content) {
    SkipInputs(c, skipped);
  }
  for (Node &c : n.fallback) {
    SkipInputs(c, skipped);
  }
  bool hit = false;
  for (int in : n.inputs) {
    hit = hit || std::find(skipped.begin(), skipped.end(), in) !=
                     skipped.end();
  }
  if (!hit) {
    return;
  }
  if (n.kind != Node::Kind::Binning) {
    throw std::runtime_error("JetSelector: can only skip plain binning cuts");
  }
  for (size_t i = 0; i < n.content.size(); i++) {
    if (IsFail(n.content[i])) {
      continue;
    }
    if (n.passChild >= 0) {
      throw std::runtime_error("JetSelector: skipped cut has more than one "
                               "passing branch");
    }
    n.passChild = (int)i;
  }
  if (n.passChild < 0) {
    throw std::runtime_error("JetSelector: skipped cut has no passing branch");
  }
  n.skip = true;
}

inline int InputIndex(const std::vector<std::string> &names,
                      const std::string &name) {
  auto it = std::find(names.begin(), names.end(), name);
  if (it == names.end()) {
    throw std::runtime_error("JetSelector: node refers to unknown input '" +
                             name + "'");
  }
  return (int)(it - names.begin());
}

inline Node Build(const JetSelectorJSON::Value &v,
                  const std::vector<std::string> &names) {
  using JetSelectorJSON::Value;
  Node n;
  if (v.type == Value::Type::Number) {
    n.value = v.number;
    return n;
  }
  if (v.type != Value::Type::Object) {
    throw std::runtime_error("JetSelector: unexpected JSON node");
  }

  const std::string type = v.At("nodetype").string;
  auto setFlow = [&](const Value &flow) {
    if (flow.type == Value::Type::String && flow.string == "clamp") {
      n.flow = Node::Flow::Clamp;
    } else if (flow.type == Value::Type::String && flow.string == "error") {
      n.flow = Node::Flow::Error;
    } else {
      n.flow = Node::Flow::Default;
      n.fallback.push_back(Build(flow, names));
    }
  };

  if (type == "binning") {
    n.kind = Node::Kind::Binning;
    n.inputs.push_back(InputIndex(names, v.At("input").string));
    std::vector<double> e;
    for (const auto &x : v.At("edges").array) {
      e.push_back(x.number);
    }
    n.edges.push_back(e);
    for (const auto &c : v.At("content").array) {
      n.content.push_back(Build(c, names));
    }
    if (n.content.size() + 1 != e.size()) {
      throw std::runtime_error("JetSelector: binning edges/content mismatch");
    }
    setFlow(v.At("flow"));
  } else if (type == "multibinning") {
    n.kind = Node::Kind::MultiBinning;
    size_t total = 1;
    for (const auto &in : v.At("inputs").array) {
      n.inputs.push_back(InputIndex(names, in.string));
    }
    for (const auto &axis : v.At("edges").array) {
      std::vector<double> e;
      for (const auto &x : axis.array) {
        e.push_back(x.number);
      }
      total *= e.size() - 1;
      n.edges.push_back(e);
    }
    for (const auto &c : v.At("content").array) {
      n.content.push_back(Build(c, names));
    }
    if (n.content.size() != total) {
      throw std::runtime_error(
          "JetSelector: multibinning edges/content mismatch");
    }
    setFlow(v.At("flow"));
  } else if (type == "category") {
    n.kind = Node::Kind::Category;
    n.inputs.push_back(InputIndex(names, v.At("input").string));
    for (const auto &item : v.At("content").array) {
      const Value &key = item.At("key");
      if (key.type != Value::Type::String) {
        throw std::runtime_error("JetSelector: only string category keys "
                                 "are supported");
      }
      n.keys.push_back(key.string);
      n.content.push_back(Build(item.At("value"), names));
    }
    if (v.Has("default") && v.At("default").type != Value::Type::Null) {
      n.fallback.push_back(Build(v.At("default"), names));
    }
  } else if (type == "transform") {
    n.kind = Node::Kind::Transform;
    n.inputs.push_back(InputIndex(names, v.At("input").string));
    const Value &rule = v.At("rule");
    if (rule.At("nodetype").string != "formula" ||
        rule.At("expression").string != "abs(x)" ||
        rule.At("variables").array.size() != 1 ||
        rule.At("variables").array[0].string != v.At("input").string) {
      throw std::runtime_error("JetSelector: only abs(x) of the transformed "
                               "input is supported");
    }
    n.content.push_back(Build(v.At("content"), names));
  } else {
    throw std::runtime_error("JetSelector: unsupported node type '" + type +
                             "'");
  }
  return n;
}

// correctionlib convention: bins are [lo, hi), x >= last edge is overflow
// returns -1 underflow, nBins overflow
inline int FindBin(const std::vector<double> &e, double x) {
  auto it = std::upper_bound(e.begin(), e.end(), x);
  if (it == e.begin()) {
    return -1;
  }
  return (int)(it - e.begin()) - 1;
}

inline double Evaluate(const Node &n, const std::vector<Input> &in) {
  switch (n.kind) {
  case Node::Kind::Leaf:
    return n.value;
  case Node::Kind::Transform: {
    std::vector<Input> out = in;
    out[n.inputs[0]].number = std::abs(out[n.inputs[0]].number);
    return Evaluate(n.content[0], out);
  }
  case Node::Kind::Category: {
    const std::string &key = in[n.inputs[0]].string;
    for (size_t i = 0; i < n.keys.size(); i++) {
      if (n.keys[i] == key) {
        return Evaluate(n.content[i], in);
      }
    }
    if (!n.fallback.empty()) {
      return Evaluate(n.fallback[0], in);
    }
    throw std::runtime_error("JetSelector: category key '" + key +
                             "' not found");
  }
  case Node::Kind::Binning:
  case Node::Kind::MultiBinning: {
    if (n.skip) {
      return Evaluate(n.content[n.passChild], in);
    }
    size_t index = 0;
    for (size_t d = 0; d < n.inputs.size(); d++) {
      const std::vector<double> &e = n.edges[d];
      int nBins = (int)e.size() - 1;
      int b = FindBin(e, in[n.inputs[d]].number);
      if (b < 0 || b >= nBins) {
        if (n.flow == Node::Flow::Default) {
          return Evaluate(n.fallback[0], in);
        }
        if (n.flow == Node::Flow::Error) {
          throw std::runtime_error("JetSelector: input out of range, value " +
                                   std::to_string(in[n.inputs[d]].number));
        }
        b = (b < 0) ? 0 : nBins - 1;
      }
      // row-major, last axis fastest
      index = index * nBins + b;
    }
    return Evaluate(n.content[index], in);
  }
  }
  return 0.0;
}

struct Correction {
  std::string name;
  std::vector<std::string> inputNames;
  Node root;

  Correction() = default;
  explicit Correction(const JetSelectorJSON::Value &v) {
    name = v.At("name").string;
    for (const auto &in : v.At("inputs").array) {
      inputNames.push_back(in.At("name").string);
    }
    root = Build(v.At("data"), inputNames);
  }

  int Index(const std::string &inputName) const {
    return InputIndex(inputNames, inputName);
  }

  double Evaluate(const std::vector<Input> &in) const {
    if (in.size() != inputNames.size()) {
      throw std::runtime_error("JetSelector: " + name + " expects " +
                               std::to_string(inputNames.size()) + " inputs");
    }
    return JetSelectorCorrection::Evaluate(root, in);
  }
};

} // namespace JetSelectorCorrection

class JetSelector {
public:
  enum class System { pp, Ion };
  enum class Purpose { Analysis, Calibration };

  // Run 3 AK4 CHS recommendations, cms-jme-jmar/jerc docs
  static constexpr const char *kIDKey = "AK4CHS_TightLeptonVeto";
  static constexpr const char *kAnalysisMap = "jetvetomap";
  static constexpr const char *kCalibrationMap = "jetvetomap_all";

  // Run 3 event veto minimal jet selection
  static constexpr double kEventVetoMinPt = 15.0;
  static constexpr double kEventVetoMaxEMF = 0.9;

  JetSelector(System system, const std::string &jetIDFile,
              const std::string &vetoMapFile,
              Purpose purpose = Purpose::Analysis)
      : system_(system), purpose_(purpose),
        vetoMapType_(purpose == Purpose::Calibration ? kCalibrationMap
                                                     : kAnalysisMap) {
    LoadJetID(jetIDFile);
    LoadVetoMap(vetoMapFile);
  }

  // multiplicities in forest convention, summed like pat::Jet
  bool PassesID(double eta, double CHF, double NHF, double CEF, double NEF,
                double MUF, int CHM, int NHM, int CEM, int NEM,
                int MUM) const {
    int chMult = CHM + CEM + MUM;
    int neMult = NHM + NEM;
    int mult = chMult + neMult;

    std::vector<JetSelectorCorrection::Input> in(id_.inputNames.size());
    in[idEta_].number = eta;
    in[idCHF_].number = CHF;
    in[idNHF_].number = NHF;
    in[idCEF_].number = CEF;
    in[idNEF_].number = NEF;
    in[idMUF_].number = MUF;
    in[idChMult_].number = chMult;
    in[idNeMult_].number = neMult;
    in[idMult_].number = mult;
    return id_.Evaluate(in) != 0.0;
  }

  bool IsVeto(double eta, double phi) const {
    std::vector<JetSelectorCorrection::Input> in(veto_.inputNames.size());
    in[vetoType_].string = vetoMapType_;
    in[vetoEta_].number = eta;
    in[vetoPhi_].number = phi;
    return veto_.Evaluate(in) != 0.0;
  }

  // per jet: passes ID and outside vetoed regions
  bool JetSelection(double eta, double phi, double CHF, double NHF, double CEF,
                    double NEF, double MUF, int CHM, int NHM, int CEM, int NEM,
                    int MUM) const {
    return PassesID(eta, CHF, NHF, CEF, NEF, MUF, CHM, NHM, CEM, NEM, MUM) &&
           !IsVeto(eta, phi);
  }

  // per event: true if this jet should veto the whole event, i.e. all of
  //   corrected pT > kEventVetoMinPt (15 GeV)
  //   passes the jet ID (leptonic-only in Ion mode)
  //   CEF + NEF < kEventVetoMaxEMF (0.9), not an electron/photon
  //   inside a vetoed region of the map
  // Run 3 recommendation, call for every jet, reject event on any true
  bool VetoEvent(double pt, double eta, double phi, double CHF, double NHF,
                 double CEF, double NEF, double MUF, int CHM, int NHM, int CEM,
                 int NEM, int MUM) const {
    if (pt <= kEventVetoMinPt) {
      return false;
    }
    if (CEF + NEF >= kEventVetoMaxEMF) {
      return false;
    }
    if (!PassesID(eta, CHF, NHF, CEF, NEF, MUF, CHM, NHM, CEM, NEM, MUM)) {
      return false;
    }
    return IsVeto(eta, phi);
  }

  // whole event: forest-style arrays of n jets, true if any jet vetoes it
  template <typename P, typename F, typename M>
  bool VetoEvent(int n, const P *pt, const F *eta, const F *phi, const F *CHF,
                 const F *NHF, const F *CEF, const F *NEF, const F *MUF,
                 const M *CHM, const M *NHM, const M *CEM, const M *NEM,
                 const M *MUM) const {
    for (int j = 0; j < n; j++) {
      if (VetoEvent(pt[j], eta[j], phi[j], CHF[j], NHF[j], CEF[j], NEF[j],
                    MUF[j], CHM[j], NHM[j], CEM[j], NEM[j], MUM[j])) {
        return true;
      }
    }
    return false;
  }

  System GetSystem() const { return system_; }
  Purpose GetPurpose() const { return purpose_; }
  const std::string &VetoMapTag() const { return veto_.name; }
  const std::string &VetoMapType() const { return vetoMapType_; }

private:
  System system_;
  Purpose purpose_;
  std::string vetoMapType_;

  JetSelectorCorrection::Correction id_;
  JetSelectorCorrection::Correction veto_;

  int idEta_ = 0, idCHF_ = 0, idNHF_ = 0, idCEF_ = 0, idNEF_ = 0, idMUF_ = 0;
  int idChMult_ = 0, idNeMult_ = 0, idMult_ = 0;
  int vetoType_ = 0, vetoEta_ = 0, vetoPhi_ = 0;

  void LoadJetID(const std::string &path) {
    JetSelectorJSON::Value file = JetSelectorJSON::ParseFile(path);
    for (const auto &c : file.At("corrections").array) {
      if (c.At("name").string == kIDKey) {
        id_ = JetSelectorCorrection::Correction(c);
        idEta_ = id_.Index("eta");
        idCHF_ = id_.Index("chHEF");
        idNHF_ = id_.Index("neHEF");
        idCEF_ = id_.Index("chEmEF");
        idNEF_ = id_.Index("neEmEF");
        idMUF_ = id_.Index("muEF");
        idChMult_ = id_.Index("chMultiplicity");
        idNeMult_ = id_.Index("neMultiplicity");
        idMult_ = id_.Index("multiplicity");

        // ion: leptonic/EM criteria only, hadronic cuts skipped
        if (system_ == System::Ion) {
          JetSelectorCorrection::SkipInputs(
              id_.root, {idCHF_, idNHF_, idChMult_, idNeMult_, idMult_});
        }
        return;
      }
    }
    throw std::runtime_error("JetSelector: " + std::string(kIDKey) +
                             " not found in " + path);
  }

  void LoadVetoMap(const std::string &path) {
    JetSelectorJSON::Value file = JetSelectorJSON::ParseFile(path);
    const auto &corrections = file.At("corrections").array;
    if (corrections.size() != 1) {
      throw std::runtime_error("JetSelector: expected one veto map "
                               "correction in " +
                               path);
    }
    veto_ = JetSelectorCorrection::Correction(corrections[0]);
    vetoType_ = veto_.Index("type");
    vetoEta_ = veto_.Index("eta");
    vetoPhi_ = veto_.Index("phi");

    // fail at construction, not on the first jet
    const auto &keys = veto_.root.keys;
    if (veto_.root.kind != JetSelectorCorrection::Node::Kind::Category ||
        std::find(keys.begin(), keys.end(), vetoMapType_) == keys.end()) {
      throw std::runtime_error("JetSelector: veto map type '" + vetoMapType_ +
                               "' not found in " + path);
    }
  }
};

#endif
