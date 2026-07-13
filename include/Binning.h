#ifndef BINNING_H
#define BINNING_H

#include "TAxis.h"
#include "TString.h"
#include "THnSparse.h"
#include "Rtypes.h"

#include "Colors.h"
#include "AnalysisConfig.h"

#include <vector>

struct AxisBins {
  Int_t nBins;
  Float_t lo;
  Float_t hi;
  TString title = "";
};

struct RangeBin {
  Float_t lo;
  Float_t hi;
  TString title;
  TString shortName;
  Color_t color = kBlack;
};

inline const std::vector<Double_t> kAbsEtaEdges = {
    0,     0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 2.172,
    2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};

inline const std::vector<Double_t> kEtaEdges = {
    -5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.650, -2.500,
    -2.322, -2.172, -1.930, -1.653, -1.479, -1.305, -1.044, -0.783,
    -0.522, -0.261, 0,      0.261,  0.522,  0.783,  1.044,  1.305,
    1.479,  1.653,  1.930,  2.172,  2.322,  2.500,  2.650,  2.853,
    2.964,  3.139,  3.489,  3.839,  5.191};

// Coarse pileup-density (rho) bins for the response THnSparse. 1-GeV bins
// across the populated 0-6 region (~98% of pp-ref jets), one wide catch-all
// so no jet overflows the axis -- inclusive-in-rho projections stay complete.
inline const std::vector<Double_t> kRhoEdges = {0.0, 1.0, 2.0, 3.0,
                                                4.0, 6.0, 25.0};

// Coarse |eta| detector regions for the JES/JER-vs-pT_gen overlay plots
// (barrel / endcap / transition / forward-HF). Edges are real kAbsEtaEdges
// members, so each slice sums whole fine bins -- kept purpose-built and
// separate from the fine residual binning so the overlay legend stays legible.
inline const std::vector<Double_t> kEtaRangeEdges = {0.0, 1.305, 2.5, 2.964,
                                                     5.191};

// Derive AxisBins from an edge vector
inline AxisBins AxisBinsFromEdges(const std::vector<Double_t> &edges,
                                  const TString &title = "") {
  return {(Int_t)edges.size() - 1, (Float_t)edges.front(),
          (Float_t)edges.back(), title};
}

inline void SetAbsEtaBins(THnSparse *h, Int_t axis) {
  h->GetAxis(axis)->Set((Int_t)kAbsEtaEdges.size() - 1, kAbsEtaEdges.data());
}

inline void SetEtaBins(THnSparse *h, Int_t axis) {
  h->GetAxis(axis)->Set((Int_t)kEtaEdges.size() - 1, kEtaEdges.data());
}

inline void SetRhoBins(THnSparse *h, Int_t axis) {
  h->GetAxis(axis)->Set((Int_t)kRhoEdges.size() - 1, kRhoEdges.data());
}

// Build pT_avg RangeBins from cfg.ptavgEdges (cfg/2024ppRef.toml, [binning] ptavg_edges)
inline std::vector<RangeBin> BuildPtAvgSlices(const std::vector<float> &edges) {
  static Color_t (*const kColorCycle[])() = {
      HiroshigeYellow, HiroshigeIceBlue,  HiroshigeLightBlue,
      HiroshigeBlue,   HiroshigeGrayBlue, HiroshigeNightBlue};
  static constexpr int kNColors = sizeof(kColorCycle) / sizeof(kColorCycle[0]);

  std::vector<RangeBin> slices;
  const int n = (int)edges.size() - 1;
  for (int i = 0; i < n; i++) {
    RangeBin sl;
    sl.lo = edges[i];
    sl.hi = edges[i + 1];
    sl.title = Form("%g < p_{T,avg} < %g GeV", sl.lo, sl.hi);
    sl.shortName = Form("_ptavg_%d_%d", (int)sl.lo, (int)sl.hi);
    sl.color = kColorCycle[i % kNColors]();
    slices.push_back(sl);
  }
  return slices;
}

// Coarse |eta| slices for the JES/JER detector-region overlay: title + color
// per slice, ordered cool (central) -> warm (forward). Object keys are built
// from lo/hi via L2Name::EtaKey at use time, so shortName isn't needed here.
inline std::vector<RangeBin> BuildEtaRangeSlices() {
  static Color_t (*const kColorCycle[])() = {
      HiroshigeNightBlue, HiroshigeLightBlue, HiroshigeOrange,
      HiroshigeLightRed};
  static constexpr int kNColors = sizeof(kColorCycle) / sizeof(kColorCycle[0]);

  std::vector<RangeBin> slices;
  const int n = (int)kEtaRangeEdges.size() - 1;
  for (int i = 0; i < n; i++) {
    RangeBin sl;
    sl.lo = kEtaRangeEdges[i];
    sl.hi = kEtaRangeEdges[i + 1];
    sl.title = Form("%.1f #leq |#eta| < %.1f", sl.lo, sl.hi);
    sl.shortName = "";
    sl.color = kColorCycle[i % kNColors]();
    slices.push_back(sl);
  }
  return slices;
}

// Cumulative (nested-from-zero) variant of BuildEtaRangeSlices, same edge
// values: [0, edge_i] for each edge_i > 0 in kEtaRangeEdges, e.g. "just
// barrel", "barrel+endcap", ..., up to the full range -- lets a comparison
// show how a quantity changes as progressively more of the endcap/forward
// region gets folded in, alongside (not instead of) the disjoint detector
// regions BuildEtaRangeSlices already covers.
inline std::vector<RangeBin> BuildEtaRangeSlicesCumulative() {
  static Color_t (*const kColorCycle[])() = {
      HiroshigeNightBlue, HiroshigeLightBlue, HiroshigeOrange,
      HiroshigeLightRed};
  static constexpr int kNColors = sizeof(kColorCycle) / sizeof(kColorCycle[0]);

  std::vector<RangeBin> slices;
  const int n = (int)kEtaRangeEdges.size() - 1;
  for (int i = 0; i < n; i++) {
    RangeBin sl;
    sl.lo = kEtaRangeEdges[0];
    sl.hi = kEtaRangeEdges[i + 1];
    sl.title = Form("|#eta| < %.1f", sl.hi);
    sl.shortName = "";
    sl.color = kColorCycle[i % kNColors]();
    slices.push_back(sl);
  }
  return slices;
}

struct BinningConfig {

  // 4D THnSparse axes
  AxisBins eta = AxisBinsFromEdges(kEtaEdges, "#eta_{probe}");
  AxisBins ptavg = {990, 10.0, 1000.0, "p_{T,avg} [GeV/c]"};
  AxisBins alpha = {50, 0.0, 0.5, "#alpha = p_{T,3} / p_{T,avg}"};
  AxisBins asymmetry = {100, -1.0, 1.0, "A"};

  // |eta| axis
  AxisBins abseta = AxisBinsFromEdges(kAbsEtaEdges, "|#eta_{probe}|");

  // control histogram axes
  AxisBins vz = {40, -20.0, 20.0, "v_{z} (cm)"};
  AxisBins pt = {100, 0.0, 1000.0, "p_{T} (GeV)"};
  AxisBins phi = {64, -3.2, 3.2, "#phi (rad)"};
  AxisBins trig = {2, 0, 2, "trigger decision"};

  AxisBins response = {200, 0.0, 2.0, "p_{T}^{reco}/p_{T}^{gen}"};
  AxisBins rho = AxisBinsFromEdges(kRhoEdges, "#rho [GeV]");
  std::vector<RangeBin> ptavgSlices = BuildPtAvgSlices(Config().ptavgEdges);
  std::vector<RangeBin> etaRangeSlices = BuildEtaRangeSlices();

  // cumulative alpha ranges
  std::vector<RangeBin> alphaSlices = {
      {0.0, 0.05, "#alpha < 0.05", "_alpha_0p05", HiroshigeLightRed()},
      {0.0, 0.10, "#alpha < 0.10", "_alpha_0p10", HiroshigeOrange()},
      {0.0, 0.15, "#alpha < 0.15", "_alpha_0p15", HiroshigeLightOrange()},
      {0.0, 0.20, "#alpha < 0.20", "_alpha_0p20", HiroshigeYellow()},
      {0.0, 0.25, "#alpha < 0.25", "_alpha_0p25", HiroshigeIceBlue()},
      {0.0, 0.30, "#alpha < 0.30", "_alpha_0p30", HiroshigeLightBlue()},
      {0.0, 0.35, "#alpha < 0.35", "_alpha_0p35", HiroshigeBlue()},
      {0.0, 0.40, "#alpha < 0.40", "_alpha_0p40", HiroshigeGrayBlue()},
      {0.0, 0.45, "#alpha < 0.45", "_alpha_0p45", HiroshigeNightBlue()},
  };
};

#endif
