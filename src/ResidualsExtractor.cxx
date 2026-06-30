#include "ResidualsExtractor.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "THnSparse.h"
#include "TDirectory.h"
#include "TString.h"
#include "TMath.h"

#include "Binning.h"
#include "Naming.h"
#include "Utilities.h"
#include "ProgressBar.h"
#include "2024ppRef.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

// global constants
static constexpr int kMinEntries = 100;
static constexpr int kNAlphaFit = 6;     // alpha thresholds 0.05–0.30 used for linear fit
static constexpr double kGausFitHW = 0.5;
static constexpr double kMaxAbsA_fit = 0.9;
static constexpr double kAlphaFitHi = 0.31;  // alpha fitting max > 0.30

// axes
static constexpr int kEtaAxis = 0;
static constexpr int kPtAvgAxis = 1;
static constexpr int kAlphaAxis = 2;
static constexpr int kAAxis = 3;

// methods 
static constexpr int     kNMethods      = 3;
static const char* const kMethodNames[] = { "gauss", "trunc90", "trunc95" };

// ---- result structs ----

struct GaussResult {
    double mean = 0, meanErr = 0, sigma = 0, chi2ndf = -1;
    bool valid = false;
};

struct TruncResult {
    double mean = 0, meanErr = 0, nEff = 0;
    bool valid = false;
};

// ---- Gaussian fit in [-kGausFitHW, +kGausFitHW] ----

static GaussResult FitGauss( TH1D* h ){
    GaussResult r;
    if( !CanFit( h, kMinEntries ) ) return r;

    TF1* g = new TF1( Form( "_gf_%s", h->GetName() ), "gaus", -kGausFitHW, kGausFitHW );
    g->SetParameter( 0, h->GetMaximum() );
    g->SetParameter( 1, h->GetMean() );
    g->SetParameter( 2, std::max( h->GetRMS(), 1e-3 ) );

    TFitResultPtr res = h->Fit( g, "NQSR" );
    if( res.Get() && res->IsValid() ){
        r.mean    = res->Parameter( 1 );
        r.meanErr = res->ParError( 1 );
        r.sigma   = res->Parameter( 2 );
        r.chi2ndf = ( res->Ndf() > 0 ) ? res->Chi2() / res->Ndf() : -1.0;
        r.valid   = true;
    }
    delete g;
    return r;
}

// ---- truncated mean ----

// Find the bin range that contains the central `fraction` of h's area.
// Returns {1, 0} (lo > hi = invalid) if h has fewer than kMinEntries.
static std::pair<int,int> FindTruncBins( TH1D* h, double fraction ){
    if( !CanFit( h, kMinEntries ) ) return {1, 0};
    double total = h->Integral();
    if( total <= 0 ) return {1, 0};
    double tailN = 0.5 * ( 1.0 - fraction ) * total;
    int nBins = h->GetNbinsX();
    double running = 0;
    int binLo = 1;
    for( int b = 1; b <= nBins; b++ ){
        running += h->GetBinContent( b );
        if( running > tailN ){ binLo = b; break; }
    }
    running = 0;
    int binHi = nBins;
    for( int b = nBins; b >= 1; b-- ){
        running += h->GetBinContent( b );
        if( running > tailN ){ binHi = b; break; }
    }
    return {binLo, binHi};
}

// Compute the mean of h restricted to bins [binLo, binHi].
// Converts nEff to entry-equivalent count so MC histograms filled with XS weights
// (where Integral() << GetEntries()) are handled correctly.
static TruncResult TruncMeanInRange( TH1D* h, int binLo, int binHi ){
    TruncResult r;
    if( binLo > binHi ) return r;
    double nEff = h->Integral( binLo, binHi );
    double hTotal = h->Integral();
    double nEffEntries = ( hTotal > 1e-10 ) ? h->GetEntries() * nEff / hTotal : nEff;
    if( nEffEntries < 10 ) return r;
    h->GetXaxis()->SetRange( binLo, binHi );
    double mean = h->GetMean();
    double rms  = h->GetRMS();
    h->GetXaxis()->SetRange( 0, 0 );
    r.mean    = mean;
    r.meanErr = rms / TMath::Sqrt( nEffEntries );
    r.nEff    = nEffEntries;
    r.valid   = true;
    return r;
}

// ---- <A> → R and propagated error ----

static double ToR( double A )               { return ( 1.0 + A ) / ( 1.0 - A ); }
static double ToRErr( double A, double dA ){ return dA * 2.0 / ( ( 1.0 - A ) * ( 1.0 - A ) ); }

static void ResetRange( THnSparse* h, int axis ){ h->GetAxis( axis )->SetRange( 0, 0 ); }

// ============================================================
// ExtractAndFit — runs the full extraction loop for one set of sparses.
//
// etaEdges:   kAbsEtaEdges (18 bins) or kEtaEdges (36 bins)
// nameSuffix: "" for |eta|, "_fulleta" for full eta
// pb:         progress bar updated once per pT slice
//
// Output names follow the global key order:
//   cone, object kind, eta mode/bin, ptavg, alpha, method/detail.
//
// Outputs per (method, ptSlice):
//   {cone}_intercept_{etaMode}_{ptSlice}_{method}
//   {cone}_intercept_{etaMode}_{ptSlice}_{method}_norm
//
// Outputs in dRvals per (method, alphaSlice, ptSlice):
//   {cone}_R_data_{etaMode}_{ptSlice}_{alphaSlice}_{method}
//   {cone}_R_mc_{etaMode}_{ptSlice}_{alphaSlice}_{method}
//
// Outputs in dGraphs per (method, ptSlice, etaBin):
//   TGraphErrors of R_MC/R_data vs alpha (all 9 bins), with fit function [0,0.31] embedded
// ============================================================

static void ExtractAndFit(
    THnSparse* hData, THnSparse* hMC,
    const TString& cone,
    const BinningConfig& bins,
    TDirectory* dQA_data, TDirectory* dQA_mc,
    TDirectory* dGraphs,
    TDirectory* dRvals,
    TDirectory* dOut,
    const std::vector<Double_t>& etaEdges,
    const TString& nameSuffix,
    ProgressBar& pb ){
    const int nPt    = ( int )bins.ptavgSlices.size();
    const int nAlpha = ( int )bins.alphaSlices.size();
    const int nEta   = hData->GetAxis( kEtaAxis )->GetNbins();
    const bool fullEta = !nameSuffix.IsNull();
    const TString etaMode = L2Name::EtaModeKey( fullEta );

    struct RPoint { double alpha, val, err; };
    // rpts[method][ipt][ieta] — all nAlpha bins; "QR" fit selects only those within [0, kAlphaFitHi]
    std::vector<std::vector<std::vector<std::vector<RPoint>>>>
        rpts( kNMethods,
            std::vector<std::vector<std::vector<RPoint>>>( nPt,
                std::vector<std::vector<RPoint>>( nEta ) ) );

    // R_data and R_mc TH1Ds per (method, ialpha, ipt) — eta on x-axis
    using RH = std::vector<std::vector<std::vector<TH1D*>>>;
    RH hRd( kNMethods, std::vector<std::vector<TH1D*>>( nAlpha, std::vector<TH1D*>( nPt, nullptr ) ) );
    RH hRm( kNMethods, std::vector<std::vector<TH1D*>>( nAlpha, std::vector<TH1D*>( nPt, nullptr ) ) );

    for( int m = 0; m < kNMethods; m++ ){
        for( int ia = 0; ia < nAlpha; ia++ ){
            for( int ip = 0; ip < nPt; ip++ ){
                TString bn = L2Name::ObjectName( cone, "R_data",
                    {etaMode, L2Name::PtKey( bins.ptavgSlices[ip] ), L2Name::AlphaKey( bins.alphaSlices[ia] )},
                    {kMethodNames[m]} );
                hRd[m][ia][ip] = new TH1D( bn, "", ( int )etaEdges.size() - 1, etaEdges.data() );
                bn = L2Name::ObjectName( cone, "R_mc",
                    {etaMode, L2Name::PtKey( bins.ptavgSlices[ip] ), L2Name::AlphaKey( bins.alphaSlices[ia] )},
                    {kMethodNames[m]} );
                hRm[m][ia][ip] = new TH1D( bn, "", ( int )etaEdges.size() - 1, etaEdges.data() );
            }
        }
    }

    // ---- main extraction loop ----

    for( int ipt = 0; ipt < nPt; ipt++ ){
        const auto& ptSlice = bins.ptavgSlices[ipt];

        for( int ialpha = 0; ialpha < nAlpha; ialpha++ ){
            const auto& aSlice = bins.alphaSlices[ialpha];

            hData->GetAxis( kPtAvgAxis )->SetRangeUser( ptSlice.lo, ptSlice.hi );
            hMC  ->GetAxis( kPtAvgAxis )->SetRangeUser( ptSlice.lo, ptSlice.hi );
            hData->GetAxis( kAlphaAxis )->SetRangeUser( aSlice.lo,  aSlice.hi );
            hMC  ->GetAxis( kAlphaAxis )->SetRangeUser( aSlice.lo,  aSlice.hi );

            // One 2D (eta, A) projection per (ipt, ialpha) instead of nEta 1D projections.
            // Projection(yDim, xDim) → TH2D with x=eta, y=A.
            TH2D* h2Data;
            TH2D* h2MC;
            {
                TDirectory::TContext nodir( nullptr );
                h2Data = hData->Projection( kAAxis, kEtaAxis );
                h2MC   = hMC  ->Projection( kAAxis, kEtaAxis );
            }

            for( int ieta = 0; ieta < nEta; ieta++ ){
                TString etaKey   = L2Name::EtaKey( ieta );
                TString ptKey    = L2Name::PtKey( ptSlice );
                TString alphaKey = L2Name::AlphaKey( aSlice );

                TH1D* hAData;
                TH1D* hAMC;
                {
                    TDirectory::TContext nodir( nullptr );
                    hAData = h2Data->ProjectionY(
                        L2Name::ObjectName( cone, "A_data", {etaMode, etaKey, ptKey, alphaKey} ),
                        ieta + 1, ieta + 1 );
                    hAMC = h2MC->ProjectionY(
                        L2Name::ObjectName( cone, "A_mc",   {etaMode, etaKey, ptKey, alphaKey} ),
                        ieta + 1, ieta + 1 );
                }

                dQA_data->cd(); hAData->Write();
                dQA_mc  ->cd(); hAMC  ->Write();

                GaussResult gd   = FitGauss( hAData );
                GaussResult gm   = FitGauss( hAMC );
                auto [dlo90, dhi90] = FindTruncBins( hAData, 0.90 );
                auto [mlo90, mhi90] = FindTruncBins( hAMC,   0.90 );
                auto [dlo95, dhi95] = FindTruncBins( hAData, 0.95 );
                auto [mlo95, mhi95] = FindTruncBins( hAMC,   0.95 );
                TruncResult td90 = TruncMeanInRange( hAData, dlo90, dhi90 );
                TruncResult tm90 = TruncMeanInRange( hAMC,   mlo90, mhi90 );
                TruncResult td95 = TruncMeanInRange( hAData, dlo95, dhi95 );
                TruncResult tm95 = TruncMeanInRange( hAMC,   mlo95, mhi95 );

                double alphaX = aSlice.hi;

                // acrunningulate all 9 alpha bins; "QR" fit later selects only those within [0, kAlphaFitHi]
                auto acrunning = [&]( int method, double Ad, double eAd,
                                             double Am, double eAm, bool ok ){
                    if( !ok ) return;
                    if( std::abs( Ad ) > kMaxAbsA_fit || std::abs( Am ) > kMaxAbsA_fit ) return;
                    double Rd = ToR( Ad ), Rm = ToR( Am );
                    double eRd = ToRErr( Ad, eAd ), eRm = ToRErr( Am, eAm );
                    if( std::abs( Rd ) < 1e-6 ) return;
                    double ratio  = Rm / Rd;
                    double eRatio = ratio * TMath::Sqrt(
                        ( eRd/Rd )*( eRd/Rd ) + ( eRm/Rm )*( eRm/Rm ) );
                    rpts[method][ipt][ieta].push_back( {alphaX, ratio, eRatio} );
                    hRd[method][ialpha][ipt]->SetBinContent( ieta + 1, Rd );
                    hRd[method][ialpha][ipt]->SetBinError  ( ieta + 1, eRd );
                    hRm[method][ialpha][ipt]->SetBinContent( ieta + 1, Rm );
                    hRm[method][ialpha][ipt]->SetBinError  ( ieta + 1, eRm );
                };

                acrunning( 0, gd.mean,   gd.meanErr,   gm.mean,   gm.meanErr,   gd.valid  && gm.valid );
                acrunning( 1, td90.mean, td90.meanErr, tm90.mean, tm90.meanErr, td90.valid && tm90.valid );
                acrunning( 2, td95.mean, td95.meanErr, tm95.mean, tm95.meanErr, td95.valid && tm95.valid );

                delete hAData;
                delete hAMC;
            }

            delete h2Data;
            delete h2MC;

            ResetRange( hData, kPtAvgAxis );
            ResetRange( hMC,   kPtAvgAxis );
            ResetRange( hData, kAlphaAxis );
            ResetRange( hMC,   kAlphaAxis );
        }

        pb.Update();
    }

    // ---- write R_data and R_mc histograms ----

    dRvals->cd();
    for( int m = 0; m < kNMethods; m++ ){
        for( int ia = 0; ia < nAlpha; ia++ ){
            for( int ip = 0; ip < nPt; ip++ ){
                hRd[m][ia][ip]->Write();
                hRm[m][ia][ip]->Write();
                delete hRd[m][ia][ip];
                delete hRm[m][ia][ip];
            }
        }
    }

    // ---- build TGraphErrors and fit R_ratio vs alpha ----

    for( int method = 0; method < kNMethods; method++ ){
        for( int ipt = 0; ipt < nPt; ipt++ ){
            const auto& ptSlice = bins.ptavgSlices[ipt];

            TString ptKey = L2Name::PtKey( ptSlice );
            TString corrName = L2Name::ObjectName( cone, "intercept",
                {etaMode, ptKey}, {kMethodNames[method]} );

            TH1D* hCorr = new TH1D( corrName, "",
                ( int )etaEdges.size() - 1, etaEdges.data() );
            hCorr->GetXaxis()->SetTitle( nameSuffix.IsNull() ? "|#eta|" : "#eta" );
            hCorr->GetYaxis()->SetTitle( "R_{MC}/R_{data} at #alpha=0" );

            TH1D* hCorrNorm = new TH1D( corrName + "_norm", "",
                ( int )etaEdges.size() - 1, etaEdges.data() );
            hCorrNorm->GetXaxis()->SetTitle( hCorr->GetXaxis()->GetTitle() );
            hCorrNorm->GetYaxis()->SetTitle( "k_{FSR} #cdot R_{MC}/R_{data}|_{#alpha=0.30}" );

            for( int ieta = 0; ieta < nEta; ieta++ ){
                const auto& pts = rpts[method][ipt][ieta];
                int n = ( int )pts.size();
                std::vector<double> x( n ), y( n ), ex( n, 0.0 ), ey( n );
                for( int k = 0; k < n; k++ ){
                    x[k]  = pts[k].alpha;
                    y[k]  = pts[k].val;
                    ey[k] = pts[k].err;
                }

                TString gname = L2Name::ObjectName( cone, "R",
                    {etaMode, L2Name::EtaKey( ieta ), ptKey}, {kMethodNames[method]} );

                TGraphErrors* gr = new TGraphErrors( n,
                    x.data(), y.data(), ex.data(), ey.data() );
                gr->SetName( gname );
                gr->SetTitle( ";#alpha threshold;R_{MC}/R_{data}" );
                gr->SetMarkerStyle( 20 );
                gr->SetMarkerColor( ptSlice.color );
                gr->SetLineColor( ptSlice.color );

                if( !CanFit( gr, {0.0, kAlphaFitHi}, 2 ) ){ delete gr; continue; }

                // "R" option: only points within [0, kAlphaFitHi] enter the chi2 — those above
                // 0.30 are displayed in the graph but excluded from the fit
                TF1* fitFn = new TF1( gname + "_fit", "[0]+[1]*x", 0.0, kAlphaFitHi );
                fitFn->SetParameter( 0, 1.0 );
                fitFn->SetParameter( 1, 0.0 );
                fitFn->SetLineColor( ptSlice.color );
                gr->Fit( fitFn, "QR" );

                hCorr->SetBinContent( ieta + 1, fitFn->GetParameter( 0 ) );
                hCorr->SetBinError  ( ieta + 1, fitFn->GetParError( 0 ) );

                // ---- normalized variant: divide each point by the value at alpha=0.30,
                //      fit, then multiply the intercept back. Errors differ from direct
                //      method because the normalization changes the fit input distribution. ----
                double val030 = 0, err030 = 0;
                for( int k = n - 1; k >= 0; k-- ){
                    if( pts[k].alpha <= kAlphaFitHi + 1e-4 && pts[k].val > 1e-6 ){
                        val030 = pts[k].val;
                        err030 = pts[k].err;
                        break;
                    }
                }
                if( val030 > 1e-6 ){
                    // count points within fit range
                    int nfit = 0;
                    for( int k = 0; k < n && pts[k].alpha <= kAlphaFitHi + 1e-4; k++ ) nfit++;

                    std::vector<double> xn( nfit ), yn( nfit ), exn( nfit, 0.0 ), eyn( nfit );
                    bool bad = false;
                    for( int k = 0; k < nfit; k++ ){
                        if( std::abs( pts[k].val ) < 1e-6 ){ bad = true; break; }
                        xn[k]  = pts[k].alpha;
                        yn[k]  = pts[k].val / val030;
                        eyn[k] = yn[k] * TMath::Sqrt(
                            TMath::Power( pts[k].err / pts[k].val, 2.0 ) +
                            TMath::Power( err030 / val030, 2.0 ) );
                    }
                    if( !bad && nfit >= 2 ){
                        TString gnorm = gname + "_norm";
                        TGraphErrors* grn = new TGraphErrors( nfit,
                            xn.data(), yn.data(), exn.data(), eyn.data() );
                        grn->SetName( gnorm );
                        grn->SetTitle( ";#alpha threshold;R_{MC}/R_{data} (norm.)" );
                        grn->SetMarkerStyle( 20 );
                        grn->SetMarkerColor( ptSlice.color );
                        grn->SetLineColor( ptSlice.color );
                        TF1* fn = new TF1( gnorm + "_fit", "[0]+[1]*x", 0.0, kAlphaFitHi );
                        fn->SetParameter( 0, 1.0 );
                        fn->SetParameter( 1, 0.0 );
                        fn->SetLineColor( ptSlice.color );
                        grn->Fit( fn, "QR" );

                        double c0n = fn->GetParameter( 0 );
                        double ec0n = fn->GetParError( 0 );
                        double c0 = c0n * val030;
                        double ec0 = ( std::abs( c0n ) > 1e-9 )
                            ? c0 * TMath::Sqrt(
                                TMath::Power( ec0n / c0n, 2.0 ) +
                                TMath::Power( err030 / val030, 2.0 ) )
                            : ec0n * val030;

                        hCorrNorm->SetBinContent( ieta + 1, c0 );
                        hCorrNorm->SetBinError( ieta + 1, ec0 );

                        dGraphs->cd();
                        grn->Write();
                        delete fn;
                        delete grn;
                    }
                }

                dGraphs->cd();
                gr->Write();   // fit function clone is embedded in graph by ROOT's Fit()
                delete fitFn;  // delete the original (clone in graph list is separately owned)
                delete gr;
            }

            dOut->cd();
            hCorr->Write();
            hCorrNorm->Write();
            delete hCorr;
            delete hCorrNorm;
        }
    }
}

// ============================================================

void runResiduals( TString dataFile, TString mcFile, TString outputFile ){

    const AnalysisConfig& cfg = Config();

    TFile* fData = TFile::Open( dataFile, "read" );
    TFile* fMC   = TFile::Open( mcFile, "read" );
    if( !fData || fData->IsZombie() ){ std::cerr << "Cannot open " << dataFile << "\n"; return; }
    if( !fMC || fMC->IsZombie() ){ std::cerr << "Cannot open " << mcFile   << "\n"; return; }

    TFile* fOut = new TFile( outputFile, "recreate" );

    BinningConfig bins;
    const int nPt = ( int )bins.ptavgSlices.size();

    // two ExtractAndFit calls per cone (|eta| and full eta), each steps through nPt pT slices
    const int totalSteps = ( int )cfg.coneLabels.size() * 2 * nPt;
    ProgressBar pb( "Extracting:", totalSteps );

    for( const TString& cone : cfg.coneLabels ){

        TDirectory* coneDataDir = ( TDirectory* )fData->Get( cone );
        TDirectory* coneMCDir = ( TDirectory* )fMC->Get( cone );
        THnSparse* hRawData = coneDataDir ? ( THnSparse* )coneDataDir->Get( cone + "_asym" ) : nullptr;
        THnSparse* hRawMC = coneMCDir ? ( THnSparse* )coneMCDir->Get( cone + "_asym" ) : nullptr;
        if( !hRawData ){ std::cerr << "Missing " << cone << "_asym in data\n"; continue; }
        if( !hRawMC ){ std::cerr << "Missing " << cone << "_asym in MC\n"; continue; }

        THnSparse* hData = FoldEtaAxis( hRawData, kEtaAxis, cone + "_asym_data_abseta" );
        THnSparse* hMC = FoldEtaAxis( hRawMC, kEtaAxis, cone + "_asym_mc_abseta" );

        TDirectory* coneDir = fOut->mkdir( cone.Data() );
        coneDir->cd();
        hData->Write();
        hMC->Write();

        TDirectory* dQA_data = coneDir->mkdir( "QA_data" );
        TDirectory* dQA_mc = coneDir->mkdir( "QA_mc" );
        TDirectory* dQA_data_full = coneDir->mkdir( "QA_data_fulleta" );
        TDirectory* dQA_mc_full = coneDir->mkdir( "QA_mc_fulleta" );
        TDirectory* dGraphs = coneDir->mkdir( "graphs" );
        TDirectory* dRvals = coneDir->mkdir( "Rvals" );
        TDirectory* dRvals_full = coneDir->mkdir( "Rvals_fulleta" );

        ExtractAndFit( hData, hMC, cone, bins,
                      dQA_data, dQA_mc, dGraphs, dRvals, coneDir,
                      kAbsEtaEdges, "", pb );

        ExtractAndFit( hRawData, hRawMC, cone, bins,
                      dQA_data_full, dQA_mc_full, dGraphs, dRvals_full, coneDir,
                      kEtaEdges, "_fulleta", pb );

        delete hData;
        delete hMC;
    }

    pb.Finish();
    fOut->Close();
    fData->Close();
    fMC->Close();
}
