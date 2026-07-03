#ifndef L2RESIDUALS_PLOTTING_METHOD_COMPARISONS_H
#define L2RESIDUALS_PLOTTING_METHOD_COMPARISONS_H

#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"

#include "plotting/Style.h"
#include "plotting/Utilities.h"
#include "Binning.h"
#include "Naming.h"
#include "ProgressBar.h"

#include <vector>

// ============================================================
// Plot type 2: method comparison (gauss vs trunc90 vs trunc95)
//
// For each (cone, ptavg slice, eta type): one canvas with
//   top panel  — all three methods overlaid
//   bottom panel — trunc90/gauss and trunc95/gauss
// ============================================================

inline void PlotMethodComp( TFile* fIn, const TString& outDir,
                           const TString& cone, const BinningConfig& bins,
                           bool fullEta, ProgressBar& pb, bool useJer = false ){
    const TString etaMode  = L2Name::EtaModeKey( fullEta );
    const TString etaLabel = fullEta ? "Full #eta" : "|#eta|";
    const double  xMin     = fullEta ? kEtaEdges.front()    :( double )kAbsEtaEdges.front();
    const double  xMax     = fullEta ? kEtaEdges.back()     :( double )kAbsEtaEdges.back();
    const TString xTitle   = fullEta ? "#eta" : "|#eta|";
    const TString calibKey = useJer ? "jer" : "jec";

    for( int ip = 0; ip <( int )bins.ptavgSlices.size(); ip++ ){
        const auto& sl = bins.ptavgSlices[ip];

        std::vector<TH1D*> hists( kNMethods, nullptr );
        for( int m = 0; m < kNMethods; m++ ){
            TString name = L2Name::ObjectName( cone, CalibKind( "intercept", useJer ),
                {etaMode, L2Name::PtKey( sl )}, {kMethodKeys[m]} );
            hists[m] = GetHAny( fIn, {cone + "/" + name} );
        }

        if( !hists[0] ){
            for( auto* h : hists ) delete h;
            continue;
        }

        // ratios to gauss for non-null trunc methods
        std::vector<TH1D*> ratios;
        std::vector<int>   ratioIdx;
        for( int m = 1; m < kNMethods; m++ ){
            if( !hists[m] ) continue;
            TString rname = Form( "%s_mcomp%s%s_r%d",
                cone.Data(), sl.shortName.Data(), etaMode.Data(), m );
            ratios.push_back( RatioH( hists[m], hists[0], rname ) );
            ratioIdx.push_back( m );
        }

        TString ptKey = L2Name::PtKey( sl );
        const TString cvName = Form( "methods_%s_%s_%s_%s",
            cone.Data(), calibKey.Data(), etaMode.Data(), ptKey.Data() );
        TwoPad cv = MakeTwoPad( cvName );

        // ---- main pad ----
        cv.main->cd();
        cv.main->SetGridx();
        cv.main->SetGridy();

        auto [ylo, yhi] = YRange( {hists[0], hists[1], hists[2]} );
        bool first = true;

        TLegend* leg = new TLegend( 0.16, 0.14, 0.58, 0.14 + 0.065 * kNMethods );
        leg->SetBorderSize( 0 );
        leg->SetFillStyle( 0 );
        leg->SetTextSize( 0.046 );

        for( int m = 0; m < kNMethods; m++ ){
            if( !hists[m] ) continue;
            StyleH( hists[m], kMethodColors[m], kMethodStyles[m], 2.f );
            hists[m]->GetYaxis()->SetRangeUser( ylo, yhi );
            hists[m]->GetYaxis()->SetTitle( CalibYTitle( useJer ) );
            hists[m]->GetYaxis()->SetTitleSize( 0.055 );
            hists[m]->GetYaxis()->SetTitleOffset( 1.10 );
            hists[m]->GetYaxis()->SetLabelSize( 0.050 );
            hists[m]->GetXaxis()->SetLabelSize( 0.0 );
            hists[m]->GetXaxis()->SetTitle( "" );
            hists[m]->SetTitle( "" );
            hists[m]->Draw( first ? "E1" : "E1 same" );
            first = false;
            leg->AddEntry( hists[m], kMethodLabels[m], "lp" );
        }

        RefLine( cv.main, xMin, xMax, 1.0 );
        leg->Draw();

        TLatex* tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize( 0.051 );
        tex->SetTextFont( 62 );
        tex->DrawLatex( 0.16, 0.91, Form( "%s   |   %s   |   %s   |   %s",
            cone.Data(), etaLabel.Data(), sl.title.Data(), CalibTag( useJer ).Data() ) );

        // ---- ratio pad ----
        cv.ratio->cd();
        cv.ratio->SetGridx();
        cv.ratio->SetGridy();

        TLegend* rleg = new TLegend( 0.16, 0.62, 0.58, 0.95 );
        rleg->SetBorderSize( 0 );
        rleg->SetFillStyle( 0 );
        rleg->SetTextSize( 0.115 );

        bool firstR = true;
        for( int k = 0; k <( int )ratios.size(); k++ ){
            const int m = ratioIdx[k];
            StyleH( ratios[k], kMethodColors[m], kMethodStyles[m], 1.5f );
            TuneRatio( ratios[k], xTitle, "/ Gauss", 0.993, 1.007 );
            ratios[k]->Draw( firstR ? "E1" : "E1 same" );
            firstR = false;
            rleg->AddEntry( ratios[k], kMethodLabels[m], "lp" );
        }
        if( !firstR ){
            RefLine( cv.ratio, xMin, xMax, 1.0 );
            rleg->Draw();
        }

        cv.c->cd();
        SavePlot( cv.c, outDir, cone, "methods", {calibKey, etaMode, ptKey}, cvName );
        pb.Update();

        // ratios and hists were drawn on pads — canvas cascade deletes them
        // rleg, leg, tex also drawn — also cascade-deleted
        delete cv.c;
    }
}

#endif
