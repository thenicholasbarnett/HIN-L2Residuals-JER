#ifndef L2RESIDUALS_PLOTTING_ALPHA_EXTRAPOLATIONS_H
#define L2RESIDUALS_PLOTTING_ALPHA_EXTRAPOLATIONS_H

#include "TFile.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"

#include "plotting/Style.h"
#include "plotting/Utilities.h"
#include "Binning.h"
#include "Naming.h"
#include "ProgressBar.h"

// ============================================================
// Plot type 5: Alpha fit plots
//
// For each (cone, method, pT slice, eta bin): one canvas showing
//   all 9 alpha threshold points with fit line drawn only through [0, 0.31],
//   points at alpha > 0.30 are shown but outside the fit line.
// ============================================================

inline void PlotAlphaFit( TFile* fIn, const TString& outDir,
                         const TString& cone, const BinningConfig& bins,
                         ProgressBar& pb ){
    TDirectory* dGraphs = ( TDirectory* )fIn->Get( cone + "/graphs" );

    const int nPt  = ( int )bins.ptavgSlices.size();
    const int nEta = ( int )kAbsEtaEdges.size() - 1;

    for( int m = 0; m < kNMethods; m++ ){
        for( int ip = 0; ip < nPt; ip++ ){
            const auto& ptSl = bins.ptavgSlices[ip];
            for( int ie = 0; ie < nEta; ie++ ){
                TString etaKey = L2Name::EtaKey( ie, false );
                TString ptKey = L2Name::PtKey( ptSl );
                TString gname = L2Name::ObjectName( cone, "R",
                    {L2Name::EtaModeKey( false ), etaKey, ptKey}, {kMethodKeys[m]} );
                TString gnorm = gname + "_norm";

                TGraphErrors* gr = GetGraphAny( dGraphs, {gnorm, gname} );

                if( !gr || gr->GetN() < 2 ){
                    continue;
                }
                TGraphErrors* gc = ( TGraphErrors* )gr->Clone( gnorm + "_c" );

                const TString cvName = Form( "alphafit_%s_%s_%s_%s",
                    cone.Data(), kMethodKeys[m], etaKey.Data(), ptKey.Data() );
                TCanvas* c = new TCanvas( cvName, "", 800, 600 );
            RealAspectRatio( c );
                c->SetLeftMargin( 0.13 );
                c->SetGridx(); c->SetGridy();

                gc->SetMarkerStyle( 20 );
                gc->SetMarkerColor( ptSl.color );
                gc->SetLineColor( ptSl.color );
                gc->SetMarkerSize( 0.9 );

                gc->GetXaxis()->SetTitle( "#alpha threshold" );
                gc->GetYaxis()->SetTitle( "R_{MC}/R_{data}|_{#alpha} / R_{MC}/R_{data}|_{#alpha=0.30}" );
                gc->GetXaxis()->CenterTitle();
                gc->GetYaxis()->CenterTitle();
                gc->GetXaxis()->SetLimits( 0.0, 0.50 );
                gc->SetTitle( "" );

                gc->Draw( "AP" );   // embedded fit function draws automatically

                // vertical reference at x=0.30 to mark the fit boundary
                double ylo = gc->GetHistogram()->GetMinimum();
                double yhi = gc->GetHistogram()->GetMaximum();
                TLine* vl = new TLine( 0.30, ylo, 0.30, yhi );
                vl->SetLineStyle( 2 ); vl->SetLineColor( kGray+2 ); vl->SetLineWidth( 1 );
                vl->Draw();

                // horizontal reference at y=1
                TLine* hl = new TLine( 0.0, 1.0, 0.50, 1.0 );
                hl->SetLineStyle( 2 ); hl->SetLineColor( kGray+2 ); hl->SetLineWidth( 1 );
                hl->Draw();

                TLatex* tex = new TLatex();
                tex->SetNDC(); tex->SetTextSize( 0.042 ); tex->SetTextFont( 62 );
                tex->DrawLatex( 0.14, 0.92, Form( "%s  |  %s  |  %s  |  |#eta| bin %d",
                    cone.Data(), kMethodLabels[m], ptSl.title.Data(), ie ) );

                SavePlot( c, outDir, cone, "alpha", {etaKey, ptKey}, cvName );
                pb.Update();

                delete gc;
                delete c;
            }
        }
    }
}

#endif
