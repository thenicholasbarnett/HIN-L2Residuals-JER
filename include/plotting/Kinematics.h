#ifndef L2RESIDUALS_PLOTTING_KINEMATICS_H
#define L2RESIDUALS_PLOTTING_KINEMATICS_H

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TPad.h"

#include "plotting/Style.h"
#include "plotting/Utilities.h"
#include "Colors.h"
#include "ProgressBar.h"
#include "jetmet/Variables.h"   // JME-standard axis titles for reco jet pT/eta

#include <algorithm>
#include <cctype>

// ============================================================
// Plot type 7: Step-1 jet kinematics
//
// Reads the TH3D(eta, phi, pT) control histograms written by runAsymmetry:
//   {cone}_incl, {cone}_tag, {cone}_probe
//
// For each cone and collection, writes pT, eta, phi projections plus eta-phi
// maps above a few pT thresholds. This mode expects a Step-1 runAsymmetry file,
// not a Step-2 residuals file.
//
// includeIncl: when false, skips the inclusive-jet collection (both its own
// per-collection plots and its contribution to the incl/probe/tag overview
// overlays) — used by runPlotting's curated smart-default flag set, which
// only wants tag+probe by default; explicit "kinematics" or "all" still get
// all three collections.
// ============================================================

static constexpr int kNKinematicsCollections = 3;
struct KinematicsCollection {
    const char* inputKey;
    const char* plotKey;
};
static const KinematicsCollection kKinematicsCollections[] = {
    {"incl", "inclusive"},
    {"probe", "probe"},
    {"tag", "tag"}
};

static constexpr int kNKinematicsPtMins = 3;
static const double kKinematicsPtMins[] = { 40.0, 100.0, 200.0 };

// incl/tag/probe here are all genuinely reconstructed jets (Step 1's real
// data/MC jet collections, not gen-matched), so JetMET's jtpt/jteta ("Reco")
// convention is the honest equivalent for these two axes. phi has no JetMET
// Variables.hh entry at all, so it and the eta-phi map (which pairs eta with
// phi and would look inconsistent with only one axis relabeled) stay local.
static const TString kRecoPtAxisTitle  = VARIABLES::getVariableAxisTitleString( VARIABLES::jtpt, true );
static const TString kRecoEtaAxisTitle = VARIABLES::getVariableAxisTitleString( VARIABLES::jteta, false );

inline TH1D* ProjectTH3D1D( TH3D* h3, const char* axis, const TString& name ){
    TDirectory::TContext nodir( nullptr );
    TH1D* h = ( TH1D* )h3->Project3D( axis );
    h->SetName( name );
    h->SetDirectory( 0 );
    return h;
}

inline TH2D* ProjectEtaPhi( TH3D* h3, double ptMin, const TString& name ){
    h3->GetZaxis()->SetRangeUser( ptMin, h3->GetZaxis()->GetXmax() );
    TDirectory::TContext nodir( nullptr );
    TH2D* h = ( TH2D* )h3->Project3D( "yx" );
    h->SetName( name );
    h->SetDirectory( 0 );
    h3->GetZaxis()->SetRange( 0, 0 );
    return h;
}

inline void DrawKinematics1D( TH1D* h, const TString& xTitle, const TString& yTitle, bool logy ){
    if( h->Integral() > 0 ) h->Scale( 1.0 / h->Integral() );
    h->SetTitle( "" );
    h->GetXaxis()->SetTitle( xTitle );
    h->GetYaxis()->SetTitle( yTitle );
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->GetXaxis()->SetTitleOffset( 1.25 );
    StyleH( h, HiroshigeNightBlue(), 20 );
    h->Draw( "hist" );
    if( logy && h->GetMaximum() > 0 ){
        h->SetMinimum( std::max( 1e-6, h->GetMinimum( 0.0 ) * 0.5 ) );
        gPad->SetLogy();
    }
}

inline TString PtMinKey( double ptMin ){
    return Form( "ptmin_%g", ptMin );
}

inline void SaveKinematicsPlot( TCanvas* c, const TString& outDir,
                               const TString& cone, const TString& coll,
                               const TString& fileBase ){
    SavePlot( c, outDir, "", "kinematics", {cone, coll}, fileBase );
}

inline void PlotKinematics( TFile* fIn, const TString& outDir,
                           const TString& cone, bool includeIncl, ProgressBar& pb ){
    {
        // "_incl" is used purely as a marker that this is a Step-1 runAsymmetry
        // file, independent of whether inclusive-jet plots are actually drawn below.
        TH3D* hTest = ( TH3D* )fIn->Get( cone + "/" + cone + "_incl" );
        if( !hTest ){
            return;
        }
    }

    TH1D* hPtAll[kNKinematicsCollections]  = {};
    TH1D* hEtaAll[kNKinematicsCollections] = {};
    TH1D* hPhiAll[kNKinematicsCollections] = {};

    auto drawConeLabel = [&]( const TString& collPlot, double xLeft ){
        TLatex* lab = new TLatex( xLeft + 0.03, 0.855,
                                 Form( "%s  |  %s", cone.Data(), collPlot.Data() ) );
        lab->SetNDC(); lab->SetTextFont( 42 ); lab->SetTextSize( 0.035 ); lab->Draw();
    };

    for( int ic = 0; ic < kNKinematicsCollections; ic++ ){
        const TString coll = kKinematicsCollections[ic].inputKey;
        const TString collPlot = kKinematicsCollections[ic].plotKey;
        if( !includeIncl && coll == "incl" ) continue;
        TH3D* h3 = ( TH3D* )fIn->Get( cone + "/" + cone + "_" + coll );

        if( !h3 ){
            continue;
        }

        h3->GetXaxis()->SetRange( 0, 0 );
        h3->GetYaxis()->SetRange( 0, 0 );
        h3->GetZaxis()->SetRange( 0, 0 );

        {
            TString cvName = Form( "kinematics_%s_%s_pt", cone.Data(), collPlot.Data() );
            TH1D* h = ProjectTH3D1D( h3, "z", cvName + "_h" );
            TCanvas* c = new TCanvas( cvName, "", 800, 800 );
            RealAspectRatio( c );
            c->SetLeftMargin( 0.14 ); c->SetGridx(); c->SetGridy();
            DrawKinematics1D( h, kRecoPtAxisTitle, "1/N  dN/dp_{T}", true );
            hPtAll[ic] = ( TH1D* )h->Clone( Form( "hPtAll_%d", ic ) ); hPtAll[ic]->SetDirectory( 0 );
            DrawCMSInternalHeader( 0.14, 0.90 );
            drawConeLabel( collPlot, 0.14 );
            SaveKinematicsPlot( c, outDir, cone, collPlot, cvName );
            delete c;
            delete h;
            pb.Update();
        }

        {
            TString cvName = Form( "kinematics_%s_%s_eta", cone.Data(), collPlot.Data() );
            TH1D* h = ProjectTH3D1D( h3, "x", cvName + "_h" );
            TCanvas* c = new TCanvas( cvName, "", 800, 800 );
            RealAspectRatio( c );
            c->SetLeftMargin( 0.14 ); c->SetGridx(); c->SetGridy();
            DrawKinematics1D( h, kRecoEtaAxisTitle, "1/N  dN/d#eta", false );
            hEtaAll[ic] = ( TH1D* )h->Clone( Form( "hEtaAll_%d", ic ) ); hEtaAll[ic]->SetDirectory( 0 );
            DrawCMSInternalHeader( 0.14, 0.90 );
            drawConeLabel( collPlot, 0.14 );
            SaveKinematicsPlot( c, outDir, cone, collPlot, cvName );
            delete c;
            delete h;
            pb.Update();
        }

        {
            TString cvName = Form( "kinematics_%s_%s_phi", cone.Data(), collPlot.Data() );
            TH1D* h = ProjectTH3D1D( h3, "y", cvName + "_h" );
            TCanvas* c = new TCanvas( cvName, "", 800, 800 );
            RealAspectRatio( c );
            c->SetLeftMargin( 0.14 ); c->SetGridx(); c->SetGridy();
            DrawKinematics1D( h, "#phi", "1/N  dN/d#phi", false );
            hPhiAll[ic] = ( TH1D* )h->Clone( Form( "hPhiAll_%d", ic ) ); hPhiAll[ic]->SetDirectory( 0 );
            DrawCMSInternalHeader( 0.14, 0.90 );
            drawConeLabel( collPlot, 0.14 );
            SaveKinematicsPlot( c, outDir, cone, collPlot, cvName );
            delete c;
            delete h;
            pb.Update();
        }

        for( int ip = 0; ip < kNKinematicsPtMins; ip++ ){
            const double ptMin = kKinematicsPtMins[ip];
            const TString ptKey = PtMinKey( ptMin );
            TString cvName = Form( "kinematics_%s_%s_eta_phi_%s",
                cone.Data(), collPlot.Data(), ptKey.Data() );

            TH2D* h = ProjectEtaPhi( h3, ptMin, cvName + "_h" );
            if( h->Integral( "width" ) > 0 ) h->Scale( 1.0 / h->Integral( "width" ) );

            TCanvas* c = new TCanvas( cvName, "", 800, 800 );
            RealAspectRatio( c );
            c->SetLeftMargin( 0.12 );
            c->SetRightMargin( 0.16 );
            c->SetTopMargin( 0.14 );
            c->SetBottomMargin( 0.13 );

            h->SetTitle( "" );
            h->GetXaxis()->SetTitle( "#eta" );
            h->GetYaxis()->SetTitle( "#phi" );
            h->GetZaxis()->SetTitle( "1/N d^{2}N/d#etad#phi" );
            h->GetXaxis()->CenterTitle();
            h->GetYaxis()->CenterTitle();
            h->GetXaxis()->SetTitleOffset( 1.25 );
            h->Draw( "colz" );

            // top: CMS Internal (left) and 2024 pp lumi (right), above frame
            DrawCMSInternalHeader( 0.12, 0.84, 0.935 );

            // collection title centered in top margin
            TString collTitle = collPlot;
            if( !collTitle.IsNull() ) collTitle[0] = toupper( collTitle[0] );
            TLatex* tColl = new TLatex( 0.48, 0.895, Form( "%s jets", collTitle.Data() ) );
            tColl->SetNDC(); tColl->SetTextFont( 42 ); tColl->SetTextSize( 0.038 );
            tColl->SetTextAlign( 22 ); tColl->Draw();

            // bottom: pT cut (left) and cone label (right), below frame
            TLatex* tPt = new TLatex( 0.12, 0.045, Form( "p_{T} > %.0f GeV/c", ptMin ) );
            tPt->SetNDC(); tPt->SetTextFont( 42 ); tPt->SetTextSize( 0.035 );
            tPt->SetTextAlign( 11 ); tPt->Draw();

            TLatex* tCone = new TLatex( 0.84, 0.045, cone.Data() );
            tCone->SetNDC(); tCone->SetTextFont( 42 ); tCone->SetTextSize( 0.035 );
            tCone->SetTextAlign( 31 ); tCone->Draw();

            SaveKinematicsPlot( c, outDir, cone, collPlot, cvName );
            delete c;
            delete h;
            pb.Update();
        }
    }

    // Overview: all three collections overlaid on a single canvas per variable
    auto DrawOverview1D = [&]( TH1D* hArr[], const TString& xTitle, const TString& yTitle,
                              const TString& varName, bool logy ){
        bool anyValid = false;
        for( int ic = 0; ic < kNKinematicsCollections; ic++ ) if( hArr[ic] ) anyValid = true;
        if( !anyValid ) return;

        // normalize
        for( int ic = 0; ic < kNKinematicsCollections; ic++ ){
            if( hArr[ic] && hArr[ic]->Integral() > 0 ) hArr[ic]->Scale( 1.0 / hArr[ic]->Integral() );
        }

        Color_t cols[kNKinematicsCollections]    = {HiroshigeNightBlue(), HiroshigeBlue(), KlimtRed()};
        int     markers[kNKinematicsCollections] = {20, 21, 22};
        for( int ic = 0; ic < kNKinematicsCollections; ic++ ){
            if( hArr[ic] ) StyleH( hArr[ic], cols[ic], markers[ic] );
        }

        TH1D* hFrame = nullptr;
        for( int ic = 0; ic < kNKinematicsCollections; ic++ ) if( hArr[ic] ){ hFrame = hArr[ic]; break; }

        double ymax = 0;
        for( int ic = 0; ic < kNKinematicsCollections; ic++ ) if( hArr[ic] ) ymax = std::max( ymax, hArr[ic]->GetMaximum() );

        hFrame->SetTitle( "" );
        hFrame->GetXaxis()->SetTitle( xTitle );
        hFrame->GetYaxis()->SetTitle( yTitle );
        hFrame->GetXaxis()->CenterTitle();
        hFrame->GetYaxis()->CenterTitle();
        hFrame->GetXaxis()->SetTitleOffset( 1.25 );
        hFrame->SetMaximum( ymax * 1.3 );
        if( logy ) hFrame->SetMinimum( 1e-6 );

        TString cvName = Form( "kinematics_%s_overview_%s", cone.Data(), varName.Data() );
        TCanvas* c = new TCanvas( cvName, "", 800, 800 );
        RealAspectRatio( c );
        c->SetLeftMargin( 0.14 ); c->SetGridx(); c->SetGridy();
        hFrame->Draw( "hist" );
        for( int ic = 0; ic < kNKinematicsCollections; ic++ ){
            if( hArr[ic] && hArr[ic] != hFrame ) hArr[ic]->Draw( "hist same" );
        }
        if( logy ) gPad->SetLogy();

        TLegend* leg = new TLegend( 0.55, 0.72, 0.88, 0.87 );
        leg->SetBorderSize( 0 ); leg->SetFillStyle( 0 ); leg->SetTextSize( 0.034 );
        for( int ic = 0; ic < kNKinematicsCollections; ic++ ){
            if( hArr[ic] ) leg->AddEntry( hArr[ic], kKinematicsCollections[ic].plotKey, "l" );
        }
        leg->Draw();

        DrawCMSInternalHeader( 0.14, 0.90 );
        TLatex* lab = new TLatex( 0.17, 0.855, cone.Data() );
        lab->SetNDC(); lab->SetTextFont( 42 ); lab->SetTextSize( 0.035 ); lab->Draw();

        SaveKinematicsPlot( c, outDir, cone, "overview", cvName );
        delete c;
        pb.Update();
    };

    DrawOverview1D( hPtAll,  kRecoPtAxisTitle,  "1/N  dN/dp_{T}", "pt",  true );
    DrawOverview1D( hEtaAll, kRecoEtaAxisTitle, "1/N  dN/d#eta",  "eta", false );
    DrawOverview1D( hPhiAll, "#phi",            "1/N  dN/d#phi",  "phi", false );

    for( int ic = 0; ic < kNKinematicsCollections; ic++ ){
        delete hPtAll[ic];
        delete hEtaAll[ic];
        delete hPhiAll[ic];
    }
}

#endif
