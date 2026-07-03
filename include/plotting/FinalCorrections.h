#ifndef L2RESIDUALS_PLOTTING_FINAL_CORRECTIONS_H
#define L2RESIDUALS_PLOTTING_FINAL_CORRECTIONS_H

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"

#include "plotting/Style.h"
#include "plotting/Utilities.h"
#include "Binning.h"
#include "Naming.h"
#include "ProgressBar.h"
#include "jetmet/RootStyle.h"
#include "jetmet/Style.h"

#include <vector>

// ============================================================
// Plot type 6: Final extrapolated values — all pT slices overlaid
//
// For each (cone, method): two canvases
//   finals_{cone}_{method}_abseta  — R_data/R_MC at alpha→0 vs |eta|, all pT bins overlaid
//   finals_{cone}_{method}_fulleta — same vs full eta
// PlotEtaSym does the |eta|-vs-fulleta symmetry check per pT slice; this overlays pT slices.
//
// Style pilot (2026-07-03): this is the second of two plot families piloting
// the vendored JetMET setTDRStyle()/CMS_lumi() (external/jetmet/RootStyle.h,
// Style.h) instead of this repo's own SetupPlotStyle() -- deliberately picked
// because, unlike every other plot family, PlotFinals never called
// DrawCMSInternalHeader() at all (it substitutes its own descriptive
// "cone | method | xTitle | JEC/JER" title instead) -- a real test of whether
// CMS_lumi() composes cleanly alongside a plot that has no prior CMS-header
// call site to swap out. Adds CMS_lumi() as a genuinely new label (closing
// that gap) and nudges the existing custom title down slightly to avoid the
// two overlapping -- see the exact y-position note at that call site.
// setTDRStyle() replaces the *global* gStyle, so it's saved/restored around
// this function the same way as EventQA.h's pilot, to avoid leaking into
// every other (non-piloted) plot family drawn later in the same
// runPlotting() call.
// ============================================================

inline void PlotFinals( TFile* fIn, const TString& outDir,
                       const TString& cone, const BinningConfig& bins,
                       ProgressBar& pb, bool isClosure = false, bool useJer = false ){
    TStyle* prevStyle = gStyle;
    setTDRStyle();
    SetupPlotStyle();

    const TString calibKey = useJer ? "jer" : "jec";
    for( int m = 0; m < kNMethods; m++ ){
        for( int ieta = 0; ieta < 2; ieta++ ){   // 0 = |eta|, 1 = full eta
            const bool   fullEta  = ( ieta == 1 );
            const TString xTitle  = fullEta ? "#eta" : "|#eta|";
            const double  xMin    = fullEta ? kEtaEdges.front()    :( double )kAbsEtaEdges.front();
            const double  xMax    = fullEta ? kEtaEdges.back()     :( double )kAbsEtaEdges.back();
            const TString etaMode = L2Name::EtaModeKey( fullEta );

            std::vector<TH1D*> hists;
            for( const auto& ptSl : bins.ptavgSlices ){
                TString name = L2Name::ObjectName( cone, CalibKind( "intercept", useJer ),
                    {etaMode, L2Name::PtKey( ptSl )}, {kMethodKeys[m]} );
                hists.push_back( GetHAny( fIn, {cone + "/" + name} ) );
            }

            bool anyValid = false;
            for( auto* h : hists ) if( h ){ anyValid = true; break; }

            // Step 3 output stores one TH2D grid (eta vs pT_avg) per (cone, etaMode,
            // method) instead of per-slice histograms — project each pT slice's row
            // as a fallback. Y-bin (ip+1) maps 1:1 to bins.ptavgSlices[ip] since the
            // grid's Y edges are built from exactly those slices in order. No JER SF
            // equivalent exists yet (runTextFile doesn't write one) -- this fallback
            // is JEC-only regardless of useJer, so a JER request against a Step 3
            // file simply finds nothing, same as any other unmet flag.
            if( !anyValid && !useJer ){
                // Step 3 stores only whichever of direct/kFSR-norm it was run with —
                // try norm first (matches PlotAlphaFit's convention), then direct.
                TString gridName = L2Name::ObjectName( cone, "corrfinal", {etaMode}, {kMethodKeys[m]} );
                TString gridNormName = gridName + "_norm";
                TH2D* h2 = GetH2Any( fIn, {cone + "/" + gridNormName, gridNormName,
                                           cone + "/" + gridName, gridName} );
                if( h2 ){
                    for( int ip = 0; ip <( int )bins.ptavgSlices.size(); ip++ ){
                        TH1D* px;
                        {
                            TDirectory::TContext nodir( nullptr );
                            px = h2->ProjectionX( Form( "%s_px%d", gridName.Data(), ip ), ip + 1, ip + 1 );
                        }
                        px->SetDirectory( 0 );
                        hists[ip] = px;
                        anyValid = true;
                    }
                    delete h2;
                }
            }

            if( !anyValid ){
                for( auto* h : hists ) delete h;
                continue;
            }

            const TString cvName = Form( "finals_%s_%s_%s_%s",
                cone.Data(), calibKey.Data(), kMethodKeys[m], etaMode.Data() );
            TCanvas* c = new TCanvas( cvName, "", 800, 600 );
            RealAspectRatio( c );
            c->SetLeftMargin( 0.13 );
            c->SetGridx();
            c->SetGridy();

            auto [ylo, yhi] = YRange( hists );
            // Closure passes are checking R_MC/R_data ~= 1 to a tight tolerance --
            // fixed 0.95-1.05 range with 0.99/1.01 guide lines reads that much
            // more clearly than the auto-scaled range used for the correction
            // derivation itself.
            if( isClosure ){ ylo = 0.95; yhi = 1.05; }

            TLegend* leg = new TLegend( 0.60, 0.68, 0.93, 0.88 );
            leg->SetBorderSize( 0 );
            leg->SetFillStyle( 0 );
            leg->SetTextSize( 0.038 );

            bool first = true;
            for( int ip = 0; ip <( int )bins.ptavgSlices.size(); ip++ ){
                if( !hists[ip] ) continue;
                const auto& ptSl = bins.ptavgSlices[ip];
                StyleH( hists[ip], ptSl.color, kMethodStyles[m], 1.5f );
                hists[ip]->GetYaxis()->SetRangeUser( ylo, yhi );
                hists[ip]->GetYaxis()->SetTitle( CalibYTitle( useJer ) );
                hists[ip]->GetYaxis()->SetTitleSize( 0.052 );
                hists[ip]->GetYaxis()->SetTitleOffset( 1.15 );
                hists[ip]->GetYaxis()->SetLabelSize( 0.048 );
                hists[ip]->GetXaxis()->SetTitle( xTitle );
                hists[ip]->GetXaxis()->SetTitleSize( 0.052 );
                hists[ip]->GetXaxis()->SetLabelSize( 0.048 );
                hists[ip]->GetXaxis()->CenterTitle();
                hists[ip]->GetYaxis()->CenterTitle();
                hists[ip]->SetTitle( "" );
                hists[ip]->Draw( first ? "E1" : "E1 same" );
                first = false;
                leg->AddEntry( hists[ip], ptSl.title, "lp" );
            }

            TLine* rl = new TLine( xMin, 1.0, xMax, 1.0 );
            rl->SetLineStyle( 2 );
            rl->SetLineColor( kGray + 2 );
            rl->SetLineWidth( 1 );
            rl->Draw();

            if( isClosure ){
                // Distinct accent color + heavier weight than the plain gray
                // dotted gridlines (gStyle's own grid is also style 3 at every
                // 0.01 here) -- same gray would make these guide lines
                // indistinguishable from the automatic grid.
                TLine* rl99 = new TLine( xMin, 0.99, xMax, 0.99 );
                rl99->SetLineStyle( 3 );
                rl99->SetLineColor( HiroshigeLightRed() );
                rl99->SetLineWidth( 2 );
                rl99->Draw();

                TLine* rl101 = new TLine( xMin, 1.01, xMax, 1.01 );
                rl101->SetLineStyle( 3 );
                rl101->SetLineColor( HiroshigeLightRed() );
                rl101->SetLineWidth( 2 );
                rl101->Draw();
            }

            leg->Draw();

            // CMS_lumi's two-line "CMS"/"Internal" block runs from about
            // y=0.92 ("CMS" baseline) down to y=0.875 ("Internal" baseline)
            // under tdrStyle's margins -- verified visually, an initial
            // y=0.86 for the line below still collided with "Internal"'s
            // own font extent. y=0.78 clears it with real margin.
            CMS_lumi( c, 16, 11 );

            TLatex* tex = new TLatex();
            tex->SetNDC();
            tex->SetTextSize( 0.048 );
            tex->SetTextFont( 62 );
            tex->DrawLatex( 0.14, 0.78, Form( "%s  |  %s  |  %s  |  %s",
                cone.Data(), kMethodLabels[m], xTitle.Data(), CalibTag( useJer ).Data() ) );

            SavePlot( c, outDir, cone, "finals", {calibKey, etaMode}, cvName );
            pb.Update();

            delete c;   // cascade-deletes hists, leg, tex, rl (and rl99/rl101 if drawn)
        }
    }

    gStyle = prevStyle;
}

#endif
