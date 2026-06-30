#ifndef BINNING_H
#define BINNING_H

#include "TAxis.h"
#include "TString.h"
#include "THnSparse.h"
#include "Rtypes.h"

#include "Colors.h"

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
    0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930,
    2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191
};

inline const std::vector<Double_t> kEtaEdges = {
    -5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.650, -2.500, -2.322,
    -2.172, -1.930, -1.653, -1.479, -1.305, -1.044, -0.783, -0.522, -0.261, 0,
     0.261,  0.522,  0.783,  1.044,  1.305,  1.479,  1.653,  1.930,  2.172,
     2.322,  2.500,  2.650,  2.853,  2.964,  3.139,  3.489,  3.839,  5.191
};

// Derive AxisBins from an edge vector
inline AxisBins AxisBinsFromEdges( const std::vector<Double_t>& edges, const TString& title = "" ){
    return { ( Int_t )edges.size() - 1, ( Float_t )edges.front(), ( Float_t )edges.back(), title };
}

inline void SetAbsEtaBins( THnSparse* h, Int_t axis ){
    h->GetAxis( axis )->Set( ( Int_t )kAbsEtaEdges.size() - 1, kAbsEtaEdges.data() );
}

inline void SetEtaBins( THnSparse* h, Int_t axis ){
    h->GetAxis( axis )->Set( ( Int_t )kEtaEdges.size() - 1, kEtaEdges.data() );
}

struct BinningConfig {

    // 4D THnSparse axes
    AxisBins eta = AxisBinsFromEdges( kEtaEdges, "#eta_{probe}" );
    AxisBins ptavg = { 990, 10.0, 1000.0, "p_{T,avg} [GeV/c]" };
    AxisBins alpha = { 50, 0.0, 0.5, "#alpha = p_{T,3} / p_{T,avg}" };
    AxisBins asymmetry = { 100, -1.0, 1.0, "A" };

    // |eta| axis
    AxisBins abseta = AxisBinsFromEdges( kAbsEtaEdges, "|#eta_{probe}|" );

    // control histogram axes
    AxisBins vz = { 40, -20.0, 20.0, "v_{z} (cm)" };
    AxisBins pt = { 100, 0.0, 1000.0, "p_{T} (GeV)" };
    AxisBins phi = { 64, -3.2, 3.2, "#phi (rad)" };
    AxisBins trig = { 2, 0, 2, "trigger decision" };

    // pT_avg slicing
    std::vector<RangeBin> ptavgSlices = {
        { 30, 70, "30 < p_{T,avg} < 70 GeV", "_ptavg_30_70", HiroshigeYellow() },
        { 70, 100, "70 < p_{T,avg} < 100 GeV", "_ptavg_70_100", HiroshigeIceBlue() },
        { 100, 175, "100 < p_{T,avg} < 175 GeV", "_ptavg_100_175", HiroshigeLightBlue() },
        { 175, 250, "175 < p_{T,avg} < 250 GeV", "_ptavg_175_250", HiroshigeBlue() },
        { 250, 500, "250 < p_{T,avg} < 500 GeV", "_ptavg_250_500", HiroshigeGrayBlue() },
        { 500, 1000, "500 < p_{T,avg} < 1000 GeV", "_ptavg_500_1000", HiroshigeNightBlue() },
    };

    // cumulative alpha ranges
    std::vector<RangeBin> alphaSlices = {
        {0.0, 0.05, "#alpha < 0.05", "_alpha_0p05", HiroshigeLightRed() },
        {0.0, 0.10, "#alpha < 0.10", "_alpha_0p10", HiroshigeOrange() },
        {0.0, 0.15, "#alpha < 0.15", "_alpha_0p15", HiroshigeLightOrange() },
        {0.0, 0.20, "#alpha < 0.20", "_alpha_0p20", HiroshigeYellow() },
        {0.0, 0.25, "#alpha < 0.25", "_alpha_0p25", HiroshigeIceBlue() },
        {0.0, 0.30, "#alpha < 0.30", "_alpha_0p30", HiroshigeLightBlue() },
        {0.0, 0.35, "#alpha < 0.35", "_alpha_0p35", HiroshigeBlue() },
        {0.0, 0.40, "#alpha < 0.40", "_alpha_0p40", HiroshigeGrayBlue() },
        {0.0, 0.45, "#alpha < 0.45", "_alpha_0p45", HiroshigeNightBlue() },
    };
};

#endif
