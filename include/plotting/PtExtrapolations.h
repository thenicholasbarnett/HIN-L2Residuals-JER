#ifndef L2RESIDUALS_PLOTTING_PT_EXTRAPOLATIONS_H
#define L2RESIDUALS_PLOTTING_PT_EXTRAPOLATIONS_H

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
// Plot type: pT-dependence fit plots (Step 3)
//
// For each (cone, method, eta mode, eta bin): one canvas showing the merged
// HP/ZB correction points vs pT_avg with the 3-parameter fit curve drawn
// (the fit function is embedded in the graph by TextFileWriter.cxx, so it
// draws automatically — never re-fit here).
// ============================================================

inline void PlotPtFit( TFile* fIn, const TString& outDir,
                      const TString& cone, ProgressBar& pb ){
    TDirectory* dGraphs = ( TDirectory* )fIn->Get( cone + "/graphs" );

    for( int m = 0; m < kNMethods; m++ ){
        for( int ieta = 0; ieta < 2; ieta++ ){   // 0 = |eta|, 1 = full eta
            const bool fullEta = ( ieta == 1 );
            const int nEta = fullEta ? ( int )kEtaEdges.size() - 1 : ( int )kAbsEtaEdges.size() - 1;
            const TString etaMode = L2Name::EtaModeKey( fullEta );

            for( int ie = 0; ie < nEta; ie++ ){
                TString etaKey = L2Name::EtaKey( ie, fullEta );
                TString gname = L2Name::ObjectName( cone, "ptcorr", {etaMode, etaKey}, {kMethodKeys[m]} );
                TString gnorm = gname + "_norm";

                // Step 3 stores only whichever of direct/kFSR-norm it was run with —
                // try norm first (matches PlotAlphaFit's convention), then direct.
                TGraphErrors* gr = GetGraphAny( dGraphs, {gnorm, gname} );
                if( !gr || gr->GetN() < 2 ){
                    continue;
                }
                TGraphErrors* gc = ( TGraphErrors* )gr->Clone( gname + "_c" );

                const TString cvName = Form( "ptfit_%s_%s_%s_%s",
                    cone.Data(), kMethodKeys[m], etaMode.Data(), etaKey.Data() );
                TCanvas* c = new TCanvas( cvName, "", 800, 600 );
                RealAspectRatio( c );
                c->SetLeftMargin( 0.13 );
                c->SetGridx(); c->SetGridy();

                gc->SetMarkerStyle( 20 );
                gc->SetMarkerColor( kMethodColors[m] );
                gc->SetLineColor( kMethodColors[m] );
                gc->SetMarkerSize( 0.9 );

                gc->GetXaxis()->SetTitle( "p_{T,avg} [GeV]" );
                gc->GetYaxis()->SetTitle( "Correction factor" );
                gc->GetXaxis()->CenterTitle();
                gc->GetYaxis()->CenterTitle();
                gc->SetTitle( "" );

                gc->Draw( "AP" );   // embedded fit function draws automatically

                // horizontal reference at y=1
                double xlo = gc->GetXaxis()->GetXmin();
                double xhi = gc->GetXaxis()->GetXmax();
                TLine* hl = new TLine( xlo, 1.0, xhi, 1.0 );
                hl->SetLineStyle( 2 ); hl->SetLineColor( kGray+2 ); hl->SetLineWidth( 1 );
                hl->Draw();

                TLatex* tex = new TLatex();
                tex->SetNDC(); tex->SetTextSize( 0.042 ); tex->SetTextFont( 62 );
                tex->DrawLatex( 0.14, 0.92, Form( "%s  |  %s  |  %s  |  eta bin %d",
                    cone.Data(), kMethodLabels[m], etaMode.Data(), ie ) );

                SavePlot( c, outDir, cone, "ptfit", {etaMode, etaKey}, cvName );
                pb.Update();

                delete gc;
                delete c;
            }
        }
    }
}

#endif
