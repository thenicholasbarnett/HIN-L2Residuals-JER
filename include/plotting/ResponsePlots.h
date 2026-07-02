#ifndef L2RESIDUALS_PLOTTING_RESPONSE_H
#define L2RESIDUALS_PLOTTING_RESPONSE_H

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"

#include "plotting/Style.h"
#include "plotting/Utilities.h"
#include "plotting/AsymmetryDistributions.h"   // kMinEntriesPlot, StyleFit
#include "Binning.h"
#include "Naming.h"
#include "ProgressBar.h"

#include <algorithm>
#include <vector>

// ============================================================
// Plot type: runResponse (JES/JER extraction) QA + summary
//
// Response objects come in three variants riding the same 5D sparse
// (eta_reco, pT_gen, corrPt/pT_gen, jtPt/pT_gen, rawPt/pT_gen) -- "corr"
// (this framework's L2Relative+L2Residual correction), "reco" (the ntuple's
// own baked-in "jtpt" correction), "raw" (uncorrected). See
// ResponseExtractor.h/.cxx.
//
// Per (cone, collection = incl/tag/probe, variant = corr/reco/raw):
//   response_dist    — per pT_gen bin: the raw response distribution with a
//                       Gaussian guide fit redone at plot time. The
//                       extraction step (ResponseExtractor.cxx) writes only
//                       the raw histogram, never an embedded fit -- same
//                       split as Step 2's QA_data/QA_mc +
//                       AsymmetryDistributions.h's FitGaussianGuide.
//   response_summary — JES and JER vs pT_gen, one canvas per (quantity,
//                       variant) with all three collections overlaid (same
//                       incl/tag/probe 3-color scheme as Kinematics.h's
//                       overview overlay).
// Plus, per (cone, collection): a variant-comparison overlay -- corr/reco/raw
// JES (and separately JER) vs pT_gen on one canvas, to see directly how much
// response this framework's correction buys over the ntuple's own
// correction and over no correction at all.
// eta_reco-binned plots are deliberately not produced yet (extraction
// doesn't slice on that axis -- see ResponseExtractor.h).
// Expects a runResponse output file; gracefully does nothing if the input
// has no response objects (e.g. any other step's output).
// ============================================================

static const char* const kResponseCollections[] = { "incl", "tag", "probe" };
static constexpr int kNResponseCollections = 3;

static const char* const kResponseVariants[] = { "corr", "reco", "raw" };
static const char* const kResponseVariantLabels[] = {
    "p_{T}^{corr}/p_{T}^{gen}", "p_{T}^{reco}/p_{T}^{gen}", "p_{T}^{raw}/p_{T}^{gen}" };
static constexpr int kNResponseVariants = 3;

inline TF1* FitResponseGuide( TH1D* h, const TString& name, Color_t col, double halfWidth ){
    if( !h || h->GetEntries() < kMinEntriesPlot ) return nullptr;
    TF1* fit = new TF1( name, "gaus", 1.0 - halfWidth, 1.0 + halfWidth );
    fit->SetParameter( 0, h->GetMaximum() );
    fit->SetParameter( 1, h->GetMean() );
    fit->SetParameter( 2, std::max( h->GetRMS(), 1e-3 ) );
    StyleFit( fit, col );
    h->Fit( fit, "NQSR" );
    return fit;
}

inline void DrawResponseDist( TFile* fIn, const TString& outDir, const TString& cone,
                             const TString& dirName, const TString& objName,
                             const TString& collection, const TString& label,
                             const std::vector<TString>& plotKeys, const TString& xTitle,
                             double halfWidth, ProgressBar& pb ){
    TDirectory* d = ( TDirectory* )fIn->Get( cone + "/" + dirName );
    TH1D* h = GetHAny( d, { objName } );
    if( !h || h->GetEntries() < kMinEntriesPlot ){
        if( h ) delete h;
        return;
    }

    StyleH( h, HiroshigeNightBlue(), 20, 1.5f );
    h->SetTitle( "" );
    h->GetXaxis()->SetTitle( xTitle );
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->SetTitle( "Jets" );
    h->GetYaxis()->CenterTitle();
    h->GetXaxis()->SetTitleOffset( 1.25 );

    TString cvName = "response_" + objName;
    TCanvas* c = new TCanvas( cvName, "", 800, 600 );
    RealAspectRatio( c );
    c->SetLeftMargin( 0.13 );

    h->Draw( "E1" );

    TF1* fit = FitResponseGuide( h, cvName + "_fit", HiroshigeLightRed(), halfWidth );
    if( fit ) fit->Draw( "same" );

    DrawCMSInternalHeader( 0.13, 0.95 );
    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize( 0.033 );
    tex->SetTextFont( 42 );
    tex->SetTextAlign( 31 );
    tex->DrawLatex( 0.93, 0.84, Form( "%s  |  %s", cone.Data(), collection.Data() ) );
    tex->DrawLatex( 0.93, 0.79, label );
    DrawEntriesLabel( ( Long64_t )h->GetEntries(), 0.93, 0.74 );

    SavePlot( c, outDir, cone, "response_dist", plotKeys, cvName );
    pb.Update();

    delete fit;
    delete h;
    delete c;
}

inline void DrawResponseSummary( TFile* fIn, const TString& outDir, const TString& cone,
                                const TString& quantity, // "JES" or "JER"
                                const std::vector<TString>& orderedKeys,
                                const TString& xTitle, const TString& tag, ProgressBar& pb ){
    TH1D* h[kNResponseCollections] = {};
    bool any = false;
    for( int ic = 0; ic < kNResponseCollections; ic++ ){
        TString name = L2Name::ObjectName( cone, quantity, orderedKeys, { kResponseCollections[ic] } );
        h[ic] = GetHAny( fIn, { cone + "/" + name } );
        if( h[ic] ) any = true;
    }
    if( !any ) return;

    static const Color_t cols[kNResponseCollections] = { HiroshigeNightBlue(), HiroshigeBlue(), KlimtRed() };
    static const int markers[kNResponseCollections] = { 20, 21, 22 };

    TH1D* frame = nullptr;
    for( int ic = 0; ic < kNResponseCollections; ic++ ) if( h[ic] ){ frame = h[ic]; break; }

    auto [ylo, yhi] = YRange( { h[0], h[1], h[2] } );
    const double xMin = frame->GetXaxis()->GetXmin();
    const double xMax = frame->GetXaxis()->GetXmax();

    TCanvas* c = new TCanvas( "response_summary_" + cone + "_" + tag, "", 800, 600 );
    RealAspectRatio( c );
    c->SetLeftMargin( 0.14 ); c->SetGridx(); c->SetGridy();

    frame->SetTitle( "" );
    frame->GetXaxis()->SetTitle( xTitle );
    frame->GetYaxis()->SetTitle( quantity );
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
    frame->GetXaxis()->SetTitleOffset( 1.25 );
    frame->GetYaxis()->SetRangeUser( ylo, yhi );

    bool first = true;
    TLegend* leg = new TLegend( 0.60, 0.72, 0.88, 0.88 );
    leg->SetBorderSize( 0 ); leg->SetFillStyle( 0 ); leg->SetTextSize( 0.034 );
    for( int ic = 0; ic < kNResponseCollections; ic++ ){
        if( !h[ic] ) continue;
        StyleH( h[ic], cols[ic], markers[ic], 1.5f );
        h[ic]->Draw( first ? "E1" : "E1 same" );
        first = false;
        leg->AddEntry( h[ic], kResponseCollections[ic], "lp" );
    }
    leg->Draw();

    if( quantity == "JES" ){
        TLine* rl = new TLine( xMin, 1.0, xMax, 1.0 );
        rl->SetLineStyle( 2 );
        rl->SetLineColor( kGray + 2 );
        rl->SetLineWidth( 1 );
        rl->Draw();
    }

    DrawCMSInternalHeader( 0.14, 0.90 );
    TLatex* lab = new TLatex( 0.17, 0.855, cone.Data() );
    lab->SetNDC(); lab->SetTextFont( 42 ); lab->SetTextSize( 0.035 ); lab->Draw();

    SavePlot( c, outDir, cone, "response_summary", { tag }, c->GetName() );
    pb.Update();

    delete c;   // cascade-deletes h[*], leg, lab (and rl if drawn)
}

// Overlays the three response variants (corr/reco/raw) for one fixed
// collection on a single canvas -- directly answers "how much response does
// this framework's correction buy over the ntuple's own correction, and over
// no correction at all" for JES (and separately JER) vs pT_gen.
inline void DrawVariantComparison( TFile* fIn, const TString& outDir, const TString& cone,
                                  const TString& collection, const TString& quantity, // "JES" or "JER"
                                  ProgressBar& pb ){
    TH1D* h[kNResponseVariants] = {};
    bool any = false;
    for( int iv = 0; iv < kNResponseVariants; iv++ ){
        TString name = L2Name::ObjectName( cone, quantity, { kResponseVariants[iv], "vs_ptgen" }, { collection } );
        h[iv] = GetHAny( fIn, { cone + "/" + name } );
        if( h[iv] ) any = true;
    }
    if( !any ) return;

    static const Color_t cols[kNResponseVariants] = { HiroshigeNightBlue(), HiroshigeBlue(), KlimtRed() };
    static const int markers[kNResponseVariants] = { 20, 21, 22 };

    TH1D* frame = nullptr;
    for( int iv = 0; iv < kNResponseVariants; iv++ ) if( h[iv] ){ frame = h[iv]; break; }

    auto [ylo, yhi] = YRange( { h[0], h[1], h[2] } );
    const double xMin = frame->GetXaxis()->GetXmin();
    const double xMax = frame->GetXaxis()->GetXmax();

    TCanvas* c = new TCanvas( "response_variants_" + cone + "_" + collection + "_" + quantity, "", 800, 600 );
    RealAspectRatio( c );
    c->SetLeftMargin( 0.14 ); c->SetGridx(); c->SetGridy();

    frame->SetTitle( "" );
    frame->GetXaxis()->SetTitle( "p_{T}^{gen} [GeV]" );
    frame->GetYaxis()->SetTitle( quantity );
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
    frame->GetXaxis()->SetTitleOffset( 1.25 );
    frame->GetYaxis()->SetRangeUser( ylo, yhi );

    bool first = true;
    TLegend* leg = new TLegend( 0.60, 0.72, 0.88, 0.88 );
    leg->SetBorderSize( 0 ); leg->SetFillStyle( 0 ); leg->SetTextSize( 0.034 );
    for( int iv = 0; iv < kNResponseVariants; iv++ ){
        if( !h[iv] ) continue;
        StyleH( h[iv], cols[iv], markers[iv], 1.5f );
        h[iv]->Draw( first ? "E1" : "E1 same" );
        first = false;
        leg->AddEntry( h[iv], kResponseVariants[iv], "lp" );
    }
    leg->Draw();

    if( quantity == "JES" ){
        TLine* rl = new TLine( xMin, 1.0, xMax, 1.0 );
        rl->SetLineStyle( 2 );
        rl->SetLineColor( kGray + 2 );
        rl->SetLineWidth( 1 );
        rl->Draw();
    }

    DrawCMSInternalHeader( 0.14, 0.90 );
    TLatex* lab = new TLatex( 0.17, 0.855, Form( "%s  |  %s", cone.Data(), collection.Data() ) );
    lab->SetNDC(); lab->SetTextFont( 42 ); lab->SetTextSize( 0.035 ); lab->Draw();

    SavePlot( c, outDir, cone, "response_variants", { collection, quantity }, c->GetName() );
    pb.Update();

    delete c;   // cascade-deletes h[*], leg, lab (and rl if drawn)
}

inline void PlotResponse( TFile* fIn, const TString& outDir, const TString& cone,
                         double halfWidth, ProgressBar& pb ){
    // sentinel: does this cone have runResponse output at all?
    bool any = false;
    for( int ic = 0; ic < kNResponseCollections && !any; ic++ ){
        TString name = L2Name::ObjectName( cone, "JES", { "corr", "vs_ptgen" }, { kResponseCollections[ic] } );
        if( HasHAny( fIn, { cone + "/" + name } ) ) any = true;
    }
    if( !any ) return;

    for( int ic = 0; ic < kNResponseCollections; ic++ ){
        const TString collection = kResponseCollections[ic];

        // ---- QA distributions vs pt_gen, one per response variant -- bin
        // edges read back from each variant's own JES-vs-ptgen output
        // histogram axis rather than re-deriving them, so this stays correct
        // even if BinningConfig.pt ever changes ----
        for( int iv = 0; iv < kNResponseVariants; iv++ ){
            const TString variant = kResponseVariants[iv];
            TString jesName = L2Name::ObjectName( cone, "JES", { variant, "vs_ptgen" }, { collection } );
            TH1D* hJes = GetHAny( fIn, { cone + "/" + jesName } );
            if( hJes ){
                for( int ip = 1; ip <= hJes->GetNbinsX(); ip++ ){
                    const double lo = hJes->GetXaxis()->GetBinLowEdge( ip );
                    const double hi = hJes->GetXaxis()->GetBinUpEdge( ip );
                    TString objName = L2Name::ObjectName( cone, "response_" + variant,
                        { L2Name::PtGenKey( lo, hi ) }, { collection } );
                    TString label = Form( "p_{T}^{gen} bin: [%.0f, %.0f) GeV  (%s)", lo, hi, variant.Data() );
                    DrawResponseDist( fIn, outDir, cone, "QA_response_ptgen", objName, collection, label,
                        { collection, variant, "vs_ptgen" }, kResponseVariantLabels[iv], halfWidth, pb );
                }
                delete hJes;
            }
        }
    }

    // ---- per-variant summary: 3 collections overlaid per canvas ----
    for( int iv = 0; iv < kNResponseVariants; iv++ ){
        const TString variant = kResponseVariants[iv];
        DrawResponseSummary( fIn, outDir, cone, "JES", { variant, "vs_ptgen" }, "p_{T}^{gen} [GeV]", "JES_" + variant + "_vs_ptgen", pb );
        DrawResponseSummary( fIn, outDir, cone, "JER", { variant, "vs_ptgen" }, "p_{T}^{gen} [GeV]", "JER_" + variant + "_vs_ptgen", pb );
    }

    // ---- variant comparison: corr/reco/raw overlaid per canvas, one per collection ----
    for( int ic = 0; ic < kNResponseCollections; ic++ ){
        DrawVariantComparison( fIn, outDir, cone, kResponseCollections[ic], "JES", pb );
        DrawVariantComparison( fIn, outDir, cone, kResponseCollections[ic], "JER", pb );
    }
}

#endif
