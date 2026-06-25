#ifndef BINNING_H
#define BINNING_H

#include "TAxis.h"
#include "TString.h"
#include "THnSparse.h"
#include "Rtypes.h"

#include "Colors.h"

#include <vector>

// AxisBins and RangeBin are used throughout — AxisBins defines a uniform axis
// for THnSparse creation, RangeBin defines a named slice for projection

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
    0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930,
    2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191
};

inline const std::vector<Double_t> kEtaEdges = {
    -5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.650, -2.500, -2.322,
    -2.172, -1.930, -1.653, -1.479, -1.305, -1.044, -0.783, -0.522, -0.261, 0,
     0.261,  0.522,  0.783,  1.044,  1.305,  1.479,  1.653,  1.930,  2.172,
     2.322,  2.500,  2.650,  2.853,  2.964,  3.139,  3.489,  3.839,  5.191
};

// Derive AxisBins from an edge vector — nBins, lo, and hi stay in sync with the edges

inline AxisBins AxisBinsFromEdges(const std::vector<Double_t>& edges, const TString& title = ""){
    return {(Int_t)edges.size() - 1, (Float_t)edges.front(), (Float_t)edges.back(), title};
}

// After creating a THnSparse with kAbsEta/kEta AxisBins, call these to apply
// the correct variable-width CMS JEC bin edges to those axes

inline void SetAbsEtaBins(THnSparse* h, Int_t axis){
    h->GetAxis(axis)->Set((Int_t)kAbsEtaEdges.size() - 1, kAbsEtaEdges.data());
}

inline void SetEtaBins(THnSparse* h, Int_t axis){
    h->GetAxis(axis)->Set((Int_t)kEtaEdges.size() - 1, kEtaEdges.data());
}

// The main physics histogram is a 4D THnSparse filled once per event:
//   axis 0: eta_probe  — full eta, variable-width CMS JEC bins, call SetEtaBins(h, 0) after construction
//   axis 1: pT_avg     — 1 GeV granularity; coarse slices applied at extraction via SetRangeUser
//   axis 2: alpha      — filled at the event's actual alpha = pT_third / pT_avg (0 for 2-jet events)
//   axis 3: asymmetry  — A = (pT_probe - pT_tag) / (pT_probe + pT_tag)
//
// alphaSlices all have lo = 0 — they are cumulative ranges for SetRangeUser at extraction,
// not differential bins. SetRangeUser(alphaAxis, 0, threshold) gives all events with alpha < threshold.

struct BinningConfig {

    // 4D THnSparse axes — listed in axis order (0 through 3)
    AxisBins eta       = AxisBinsFromEdges(kEtaEdges,    "#eta_{probe}");
    AxisBins ptavg     = {990, 10.0, 1000.0, "p_{T,avg} (GeV)"};
    AxisBins alpha     = {50, 0.0, 0.5, "#alpha = p_{T,3} / p_{T,avg}"};
    AxisBins asymmetry = {100, -1.0, 1.0, "A"};

    // |eta| axis — used only at extraction after FoldEtaAxis, not stored in the THnSparse
    AxisBins abseta = AxisBinsFromEdges(kAbsEtaEdges, "|#eta_{probe}|");

    // control histogram axes
    AxisBins vz  = {40, -20.0, 20.0, "v_{z} (cm)"};
    AxisBins pt  = {500, 0.0, 1000.0, "p_{T} (GeV)"};
    AxisBins phi = {64, -3.2, 3.2, "#phi (rad)"};
    AxisBins trig = {2, 0, 2, "trigger decision"};

    // named pT_avg slices for extraction — applied via SetRangeUser(pT_avg axis, lo, hi)
    std::vector<RangeBin> ptavgSlices = {
        { 40,  90, "40 < p_{T,avg} < 90 GeV",    "_ptavg_40_90",   HiroshigeNightBlue  },
        { 90, 120, "90 < p_{T,avg} < 120 GeV",   "_ptavg_90_120",  HiroshigeGrayBlue   },
        {120, 190, "120 < p_{T,avg} < 190 GeV",  "_ptavg_120_190", HiroshigeBlue       },
        {190, 260, "190 < p_{T,avg} < 260 GeV",  "_ptavg_190_260", HiroshigeLightBlue  },
        {260, 1000,"260 < p_{T,avg} < 1000 GeV", "_ptavg_260_1000",HiroshigeIceBlue    },
    };

    // cumulative alpha ranges for extraction — lo is always 0
    // use SetRangeUser(alpha axis, slice.lo, slice.hi) to get all events with alpha < threshold
    std::vector<RangeBin> alphaSlices = {
        {0.0, 0.05, "#alpha < 0.05", "_alpha_0p05", HiroshigeLightRed    },
        {0.0, 0.10, "#alpha < 0.10", "_alpha_0p10", HiroshigeOrange      },
        {0.0, 0.15, "#alpha < 0.15", "_alpha_0p15", HiroshigeLightOrange },
        {0.0, 0.20, "#alpha < 0.20", "_alpha_0p20", HiroshigeYellow      },
        {0.0, 0.25, "#alpha < 0.25", "_alpha_0p25", HiroshigeIceBlue     },
        {0.0, 0.30, "#alpha < 0.30", "_alpha_0p30", HiroshigeLightBlue   },
        {0.0, 0.35, "#alpha < 0.35", "_alpha_0p35", HiroshigeBlue        },
        {0.0, 0.40, "#alpha < 0.40", "_alpha_0p40", HiroshigeGrayBlue    },
        {0.0, 0.45, "#alpha < 0.45", "_alpha_0p45", HiroshigeNightBlue   },
    };
};

#endif
