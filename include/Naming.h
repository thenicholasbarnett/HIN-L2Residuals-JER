#ifndef NAMING_H
#define NAMING_H

#include "TString.h"
#include "Binning.h"

#include <cmath>
#include <vector>

namespace L2Name {

inline TString CleanKey(const TString &s) {
  return s.BeginsWith("_") ? s(1, s.Length() - 1) : s;
}

inline TString EdgeKey(double x) {
  const double ax = std::abs(x);
  TString s = Form("%.3f", ax);
  while (s.EndsWith("0"))
    s.Chop();
  if (s.EndsWith("."))
    s += "0";
  s.ReplaceAll(".", "p");
  return x < 0.0 ? "m" + s : s;
}

inline TString EtaKey(double lo, double hi) {
  return Form("eta_%s_%s", EdgeKey(lo).Data(), EdgeKey(hi).Data());
}

inline TString EtaKey(int ieta, bool fullEta) {
  const auto &edges = fullEta ? kEtaEdges : kAbsEtaEdges;
  if (ieta < 0 || ieta + 1 >= (int)edges.size())
    return Form("eta_invalid_%d", ieta);
  return EtaKey(edges[ieta], edges[ieta + 1]);
}

inline TString PtKey(const RangeBin &ptSlice) {
  return CleanKey(ptSlice.shortName);
}

inline TString AlphaKey(const RangeBin &alphaSlice) {
  return CleanKey(alphaSlice.shortName);
}

// For a raw numeric pT_gen bin (runResponse's JES/JER-vs-pT_gen extraction,
// which bins directly off the response THnSparse's own uniform pT axis
// rather than a curated RangeBin slice list like ptavgSlices).
inline TString PtGenKey(double lo, double hi) {
  return Form("ptgen_%d_%d", (int)lo, (int)hi);
}

inline TString EtaModeKey(bool fullEta) {
  return fullEta ? "fulleta" : "abseta";
}

inline TString Join(const std::vector<TString> &parts,
                    const TString &sep = "_") {
  TString out;
  for (const auto &part : parts) {
    if (part.IsNull())
      continue;
    if (!out.IsNull())
      out += sep;
    out += part;
  }
  return out;
}

inline TString ObjectName(const TString &cone, const TString &objectKind,
                          const std::vector<TString> &orderedKeys = {},
                          const std::vector<TString> &detailKeys = {}) {
  std::vector<TString> parts = {cone, objectKind};
  parts.insert(parts.end(), orderedKeys.begin(), orderedKeys.end());
  parts.insert(parts.end(), detailKeys.begin(), detailKeys.end());
  return Join(parts);
}

} // namespace L2Name

#endif
