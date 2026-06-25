#ifndef BINNINGCONFIG_H
#define BINNINGCONFIG_H

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

// CMS JEC standard eta bin edges (identical to DijetResiduals binning)
// 18 |eta| bins, 36 full-eta bins — non-uniform widths

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

// After creating a THnSparse with kAbsEta/kEta AxisBins, call these to apply
// the correct variable-width CMS JEC bin edges to those axes

inline void SetAbsEtaBins(THnSparse* h, Int_t axis){
    h->GetAxis(axis)->Set((Int_t)kAbsEtaEdges.size() - 1, kAbsEtaEdges.data());
}

inline void SetEtaBins(THnSparse* h, Int_t axis){
    h->GetAxis(axis)->Set((Int_t)kEtaEdges.size() - 1, kEtaEdges.data());
}

struct BinningConfig {

    // --- THnSparse axes ---

    // pT_avg: 1 GeV bins, fine enough to use any coarse slice via SetRangeUser
    AxisBins ptavg = {990, 10.0, 1000.0, "p_{T,avg} (GeV)"};

    // alpha = pT_third / pT_avg, filled as a continuous variable
    // 2-jet events fill at alpha = 0 and naturally enter every threshold slice
    AxisBins alpha = {50, 0.0, 0.5, "#alpha = p_{T,3} / p_{T,avg}"};

    // asymmetry: the measured quantity
    AxisBins asymmetry = {100, -1.0, 1.0, "A"};

    // |eta| and full eta: created with these AxisBins, then call SetAbsEtaBins /
    // SetEtaBins to apply the correct CMS JEC variable-width edges
    AxisBins abseta = {18, 0.0, 5.191, "|#eta_{probe}|"};
    AxisBins eta    = {36, -5.191, 5.191, "#eta_{probe}"};

    // --- Control histogram axes ---

    AxisBins vz  = {40, -20.0, 20.0, "v_{z} (cm)"};
    AxisBins pt  = {500, 0.0, 1000.0, "p_{T} (GeV)"};
    AxisBins phi = {64, -3.2, 3.2, "#phi (rad)"};
    AxisBins trig = {2, 0, 2, "trigger decision"};

    // --- pT_avg slices (coarse bins for residuals extraction) ---
    // lo/hi match the ptlow/pthigh used in DijetResiduals

    std::vector<RangeBin> ptavgSlices = {
        { 40,  90, "40 < p_{T,avg} < 90 GeV",    "_ptavg_40_90",   HiroshigeNightBlue  },
        { 90, 120, "90 < p_{T,avg} < 120 GeV",   "_ptavg_90_120",  HiroshigeGrayBlue   },
        {120, 190, "120 < p_{T,avg} < 190 GeV",  "_ptavg_120_190", HiroshigeBlue       },
        {190, 260, "190 < p_{T,avg} < 260 GeV",  "_ptavg_190_260", HiroshigeLightBlue  },
        {260, 1000,"260 < p_{T,avg} < 1000 GeV", "_ptavg_260_1000",HiroshigeIceBlue    },
    };

    // --- alpha threshold slices (project as SetRangeUser(0, hi)) ---

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
