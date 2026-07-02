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
// Per (cone, collection = incl/tag/probe):
//   response_dist    — per eta_gen bin (both |eta_gen| and full eta_gen) and
//                       per pT_gen bin: the raw response distribution with a
//                       Gaussian guide fit redone at plot time. The
//                       extraction step (ResponseExtractor.cxx) writes only
//                       the raw histogram, never an embedded fit -- same
//                       split as Step 2's QA_data/QA_mc +
//                       AsymmetryDistributions.h's FitGaussianGuide.
//   response_summary — JES and JER vs eta_gen (both eta modes) and vs
//                       pT_gen, one canvas per quantity with all three
//                       collections overlaid (same incl/tag/probe 3-color
//                       scheme as Kinematics.h's overview overlay).
// Expects a runResponse output file; gracefully does nothing if the input
// has no response objects (e.g. any other step's output).
// ============================================================

static const char* const kResponseCollections[] = { "incl", "tag", "probe" };
static constexpr int kNResponseCollections = 3;

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
                             const std::vector<TString>& plotKeys,
                             double halfWidth, ProgressBar& pb ){
    TDirectory* d = ( TDirectory* )fIn->Get( cone + "/" + dirName );
    TH1D* h = GetHAny( d, { objName } );
    if( !h || h->GetEntries() < kMinEntriesPlot ){
        if( h ) delete h;
        return;
    }

    StyleH( h, HiroshigeNightBlue(), 20, 1.5f );
    h->SetTitle( "" );
    h->GetXaxis()->SetTitle( "p_{T}^{reco}/p_{T}^{gen}" );
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

inline void PlotResponse( TFile* fIn, const TString& outDir, const TString& cone,
                         double halfWidth, ProgressBar& pb ){
    // sentinel: does this cone have runResponse output at all?
    bool any = false;
    for( int ic = 0; ic < kNResponseCollections && !any; ic++ ){
        TString name = L2Name::ObjectName( cone, "JES", { "vs_ptgen" }, { kResponseCollections[ic] } );
        if( HasHAny( fIn, { cone + "/" + name } ) ) any = true;
    }
    if( !any ) return;

    for( int ic = 0; ic < kNResponseCollections; ic++ ){
        const TString collection = kResponseCollections[ic];

        // ---- QA distributions vs eta_gen (both eta modes) ----
        for( int mode = 0; mode < 2; mode++ ){
            const bool fullEta = ( mode == 1 );
            const TString etaMode = L2Name::EtaModeKey( fullEta );
            const TString dirName = fullEta ? "QA_response_fulleta" : "QA_response_abseta";
            const int nEta = fullEta ? ( int )kEtaEdges.size() - 1 : ( int )kAbsEtaEdges.size() - 1;
            for( int ie = 0; ie < nEta; ie++ ){
                TString etaKey = L2Name::EtaKey( ie, fullEta );
                TString objName = L2Name::ObjectName( cone, "response", { etaMode, etaKey }, { collection } );
                TString label = Form( "%s bin: %s", fullEta ? "eta_{gen}" : "|eta_{gen}|", etaKey.Data() );
                DrawResponseDist( fIn, outDir, cone, dirName, objName, collection, label,
                    { collection, "vs_etagen_" + etaMode }, halfWidth, pb );
            }
        }

        // ---- QA distributions vs pt_gen -- bin edges read back from the
        // JES-vs-ptgen output histogram's own axis rather than re-deriving
        // them, so this stays correct even if BinningConfig.pt ever changes ----
        {
            TString jesName = L2Name::ObjectName( cone, "JES", { "vs_ptgen" }, { collection } );
            TH1D* hJes = GetHAny( fIn, { cone + "/" + jesName } );
            if( hJes ){
                for( int ip = 1; ip <= hJes->GetNbinsX(); ip++ ){
                    const double lo = hJes->GetXaxis()->GetBinLowEdge( ip );
                    const double hi = hJes->GetXaxis()->GetBinUpEdge( ip );
                    TString objName = L2Name::ObjectName( cone, "response",
                        { L2Name::PtGenKey( lo, hi ) }, { collection } );
                    TString label = Form( "p_{T}^{gen} bin: [%.0f, %.0f) GeV", lo, hi );
                    DrawResponseDist( fIn, outDir, cone, "QA_response_ptgen", objName, collection, label,
                        { collection, "vs_ptgen" }, halfWidth, pb );
                }
                delete hJes;
            }
        }
    }

    // ---- summary plots (all 3 collections overlaid per canvas) ----
    DrawResponseSummary( fIn, outDir, cone, "JES", { "abseta", "vs_etagen" }, "|#eta_{gen}|", "JES_abseta_vs_etagen", pb );
    DrawResponseSummary( fIn, outDir, cone, "JER", { "abseta", "vs_etagen" }, "|#eta_{gen}|", "JER_abseta_vs_etagen", pb );
    DrawResponseSummary( fIn, outDir, cone, "JES", { "fulleta", "vs_etagen" }, "#eta_{gen}", "JES_fulleta_vs_etagen", pb );
    DrawResponseSummary( fIn, outDir, cone, "JER", { "fulleta", "vs_etagen" }, "#eta_{gen}", "JER_fulleta_vs_etagen", pb );
    DrawResponseSummary( fIn, outDir, cone, "JES", { "vs_ptgen" }, "p_{T}^{gen} [GeV]", "JES_vs_ptgen", pb );
    DrawResponseSummary( fIn, outDir, cone, "JER", { "vs_ptgen" }, "p_{T}^{gen} [GeV]", "JER_vs_ptgen", pb );
}

#endif
