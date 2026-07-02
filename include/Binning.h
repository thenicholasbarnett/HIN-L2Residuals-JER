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

// Build pT_avg RangeBins from cfg.ptavgEdges (cfg/2024ppRef.toml, [binning] ptavg_edges).
// title/shortName are derived from the edges so they can't drift out of sync
// with a TOML edit; shortName matches the pre-TOML hardcoded format exactly
// ("_ptavg_<lo>_<hi>") so existing Step 2 output object names stay readable.
// Colors cycle through the Hiroshige palette by slice index — plot styling
// is deliberately kept out of the physics-binning config.
inline std::vector<RangeBin> BuildPtAvgSlices( const std::vector<float>& edges ){
    static Color_t( * const kColorCycle[] )() = {
        HiroshigeYellow, HiroshigeIceBlue, HiroshigeLightBlue,
        HiroshigeBlue, HiroshigeGrayBlue, HiroshigeNightBlue
    };
    static constexpr int kNColors = sizeof( kColorCycle ) / sizeof( kColorCycle[0] );

    std::vector<RangeBin> slices;
    const int n = ( int )edges.size() - 1;
    for( int i = 0; i < n; i++ ){
        RangeBin sl;
        sl.lo = edges[i];
        sl.hi = edges[i + 1];
        sl.title = Form( "%g < p_{T,avg} < %g GeV", sl.lo, sl.hi );
        sl.shortName = Form( "_ptavg_%d_%d", ( int )sl.lo, ( int )sl.hi );
        sl.color = kColorCycle[i % kNColors]();
        slices.push_back( sl );
    }
    return slices;
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

    // MC-only: p_{T}^{reco}/p_{T}^{gen} response axis, added to the incl/tag/probe
    // kinematics THnSparses in MC mode (JES/JER inputs — see ConeHistograms).
    AxisBins response = { 200, 0.0, 2.0, "p_{T}^{reco}/p_{T}^{gen}" };

    // pT_avg slicing — edges come from cfg/2024ppRef.toml ([binning] ptavg_edges),
    // so rebinning Step 2/3's pT slices needs only a TOML edit + rerun, no recompile.
    std::vector<RangeBin> ptavgSlices = BuildPtAvgSlices( Config().ptavgEdges );

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
