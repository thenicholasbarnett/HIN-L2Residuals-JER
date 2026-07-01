#ifndef L2RESIDUALS_PLOTTING_EVENT_QA_H
#define L2RESIDUALS_PLOTTING_EVENT_QA_H

#include "TFile.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"

#include "plotting/Style.h"
#include "plotting/Utilities.h"
#include "ProgressBar.h"

#include <algorithm>

// ============================================================
// Plot type: Event-level QA
//
// Reads event-level control histograms written by runAsymmetry:
//   hvz_all   — vz before any cuts
//   hvz       — vz after vz + vertex filter cuts
//   hfilt     — pprimaryVertexFilter (DATA only)
//   h_hlt_j80 — HLT_AK4PFJet80 bits (hard-probe DATA only)
//
// Expects a Step-1 runAsymmetry output file; gracefully skips
// histograms that are absent (e.g. hfilt/h_hlt_j80 on MC).
// ============================================================

inline void PlotEvent( TFile* fIn, const TString& outDir, ProgressBar& pb ){

    // ---- vz: all events vs after cuts ----
    {
        TH1D* hall = ( TH1D* )fIn->Get( "hvz_all" );
        TH1D* hvz  = ( TH1D* )fIn->Get( "hvz" );
        if( hall && hvz ){
            TH1D* hallc = ( TH1D* )hall->Clone( "hvz_all_c" ); hallc->SetDirectory( 0 );
            TH1D* hvzc  = ( TH1D* )hvz ->Clone( "hvz_c" );     hvzc ->SetDirectory( 0 );

            hallc->SetLineColor( kGray + 2 );
            hallc->SetFillColor( kGray + 1 );
            hallc->SetFillStyle( 1001 );
            hallc->SetLineWidth( 1 );
            StyleH( hvzc, HiroshigeNightBlue(), 20, 2.0f );

            double ymax = std::max( hallc->GetMaximum(), hvzc->GetMaximum() ) * 1.2;
            hallc->SetMaximum( ymax );
            hallc->SetTitle( "" );
            hallc->GetXaxis()->SetTitle( "v_{z} [cm]" );
            hallc->GetYaxis()->SetTitle( "Events" );
            hallc->GetXaxis()->CenterTitle();
            hallc->GetYaxis()->CenterTitle();
            hallc->GetXaxis()->SetTitleOffset( 1.25 );

            TCanvas* c = new TCanvas( "event_vz", "", 800, 600 );
            RealAspectRatio( c );
            c->SetLeftMargin( 0.13 ); c->SetGridx(); c->SetGridy();
            hallc->Draw( "hist" );
            hvzc->Draw( "E1 same" );

            TLegend* leg = new TLegend( 0.55, 0.72, 0.89, 0.87 );
            leg->SetBorderSize( 0 ); leg->SetFillStyle( 0 ); leg->SetTextSize( 0.034 );
            leg->AddEntry( hallc, Form( "All events (%s)", FormatEntriesText( ( Long64_t )hall->GetEntries() ).Data() ), "f" );
            leg->AddEntry( hvzc, Form( "After cuts  (%s)", FormatEntriesText( ( Long64_t )hvz->GetEntries() ).Data() ), "lp" );
            leg->Draw();
            DrawCMSInternalHeader( 0.13, 0.90 );

            SavePlot( c, outDir, "", "event", {}, "event_vz" );
            delete c;
            pb.Update();
        }
    }

    // ---- primary vertex filter ----
    {
        TH1I* hfilt = ( TH1I* )fIn->Get( "hfilt" );
        if( hfilt ){
            TCanvas* c = new TCanvas( "event_ppvF", "", 800, 600 );
            RealAspectRatio( c );
            c->SetLeftMargin( 0.15 );
            c->SetLogy();
            hfilt->SetTitle( "" );
            hfilt->GetXaxis()->SetTitle( "pprimaryVertexFilter" );
            hfilt->GetYaxis()->SetTitle( "Events" );
            hfilt->GetXaxis()->CenterTitle();
            hfilt->GetYaxis()->CenterTitle();
            hfilt->GetXaxis()->SetTitleOffset( 1.25 );
            hfilt->SetLineColor( HiroshigeNightBlue() );
            hfilt->SetFillColor( HiroshigeLightBlue() );
            hfilt->SetFillStyle( 1001 );
            hfilt->SetMinimum( 0.5 );
            hfilt->SetMaximum( std::max( 1.0, hfilt->GetMaximum() ) * 10.0 );
            hfilt->Draw( "hist" );
            DrawCMSInternalHeader( 0.15, 0.90 );
            DrawEntriesLabel( ( Long64_t )hfilt->GetEntries() );
            SavePlot( c, outDir, "", "event", {}, "event_ppvF" );
            delete c;
            pb.Update();
        }
    }

    // ---- HLT trigger ----
    {
        TH1I* htrig = ( TH1I* )fIn->Get( "h_hlt_j80" );
        if( htrig ){
            TCanvas* c = new TCanvas( "event_hlt", "", 800, 600 );
            RealAspectRatio( c );
            c->SetLeftMargin( 0.15 );
            c->SetLogy();
            htrig->SetTitle( "" );
            htrig->GetXaxis()->SetTitle( "HLT_AK4PFJet80" );
            htrig->GetYaxis()->SetTitle( "Events" );
            htrig->GetXaxis()->CenterTitle();
            htrig->GetYaxis()->CenterTitle();
            htrig->GetXaxis()->SetTitleOffset( 1.25 );
            htrig->SetLineColor( HiroshigeNightBlue() );
            htrig->SetFillColor( HiroshigeLightBlue() );
            htrig->SetFillStyle( 1001 );
            htrig->SetMinimum( 0.5 );
            htrig->SetMaximum( std::max( 1.0, htrig->GetMaximum() ) * 10.0 );
            htrig->Draw( "hist" );
            DrawCMSInternalHeader( 0.15, 0.90 );
            DrawEntriesLabel( ( Long64_t )htrig->GetEntries() );
            SavePlot( c, outDir, "", "event", {}, "event_hlt_j80" );
            delete c;
            pb.Update();
        }
    }
}

#endif
