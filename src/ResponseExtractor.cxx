#include "ResponseExtractor.h"

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFitResult.h"
#include "THnSparse.h"
#include "TDirectory.h"
#include "TString.h"
#include "TMath.h"

#include "Binning.h"
#include "Naming.h"
#include "Utilities.h"
#include "ProgressBar.h"
#include "AnalysisConfig.h"

#include <vector>
#include <iostream>
#include <cmath>

// axis order within {cone}_incl_resp/_tag_resp/_probe_resp -- matches
// DijetHistograms.h's ConeHistograms exactly (kept as a local redefinition,
// not a DijetHistograms.h include, the same decoupling ResidualsExtractor.cxx
// already uses for the main asymmetry sparse's axes).
static constexpr int kRespEtaRecoAxis = 0;
static constexpr int kRespPtRecoAxis = 1;
static constexpr int kRespEtaGenAxis = 2;
static constexpr int kRespPtGenAxis = 3;
static constexpr int kRespRAxis = 4;

static constexpr int kNCollections = 3;
static const char* const kCollectionKeys[kNCollections] = { "incl", "tag", "probe" };
static const char* const kCollectionSuffixes[kNCollections] = { "_incl_resp", "_tag_resp", "_probe_resp" };

// ---- Gaussian fit around 1.0 -> JES (mean) and JER (sigma/mean) ----
//
// Mirrors ResidualsExtractor.cxx's FitGauss in spirit (same NQSR fit, same
// CanFit entry-count guard) but centered at the response peak (1.0), not 0,
// and additionally returns the sigma/mean ratio directly.

struct ResponseFitResult {
    double jes = 0, jesErr = 0, jer = 0, jerErr = 0;
    bool valid = false;
};

static ResponseFitResult FitResponse( TH1D* h, double halfWidth, int minEntries ){
    ResponseFitResult r;
    if( !CanFit( h, minEntries ) ) return r;

    TF1* g = new TF1( Form( "_rf_%s", h->GetName() ), "gaus", 1.0 - halfWidth, 1.0 + halfWidth );
    g->SetParameter( 0, h->GetMaximum() );
    g->SetParameter( 1, h->GetMean() );
    g->SetParameter( 2, std::max( h->GetRMS(), 1e-3 ) );

    TFitResultPtr res = h->Fit( g, "NQSR" );
    if( res.Get() && res->IsValid() ){
        const double mean = res->Parameter( 1 );
        const double meanErr = res->ParError( 1 );
        const double sigma = res->Parameter( 2 );
        const double sigmaErr = res->ParError( 2 );
        r.jes = mean;
        r.jesErr = meanErr;
        if( std::abs( mean ) > 1e-6 ){
            r.jer = sigma / mean;
            r.jerErr = r.jer * TMath::Sqrt(
                TMath::Power( sigmaErr / sigma, 2.0 ) + TMath::Power( meanErr / mean, 2.0 ) );
            r.valid = true;
        }
    }
    delete g;
    return r;
}

// ---- vs eta_gen: fullEta selects kEtaEdges (36 bins, raw sparse) vs
// kAbsEtaEdges (18 bins, caller passes a sparse already folded on
// kRespEtaGenAxis via FoldEtaAxis). eta_reco/pt_reco/pt_gen are left fully
// unrestricted (marginal, not joint, extraction -- see ResponseExtractor.h).
// Writes the raw per-bin response projection (no embedded fit -- refit at
// plot time, mirroring Step 2's QA_data/QA_mc + PlotAsymDist's
// FitGaussianGuide) to dQA, and the JES/JER TH1Ds to dOut. ----

static void ExtractVsEtaGen(
    THnSparse* h, const TString& cone, const TString& collection, bool fullEta,
    TDirectory* dQA, TDirectory* dOut, double halfWidth, int minEntries ){

    const TString etaMode = L2Name::EtaModeKey( fullEta );
    const std::vector<Double_t>& edges = fullEta ? kEtaEdges : kAbsEtaEdges;
    const int nEta = ( int )edges.size() - 1;

    TString jesName = L2Name::ObjectName( cone, "JES", { etaMode, "vs_etagen" }, { collection } );
    TString jerName = L2Name::ObjectName( cone, "JER", { etaMode, "vs_etagen" }, { collection } );
    TH1D* hJES = new TH1D( jesName, "", nEta, edges.data() );
    TH1D* hJER = new TH1D( jerName, "", nEta, edges.data() );
    hJES->GetXaxis()->SetTitle( fullEta ? "#eta_{gen}" : "|#eta_{gen}|" );
    hJES->GetYaxis()->SetTitle( "JES = #LT p_{T}^{reco}/p_{T}^{gen} #GT" );
    hJER->GetXaxis()->SetTitle( hJES->GetXaxis()->GetTitle() );
    hJER->GetYaxis()->SetTitle( "JER = #sigma / #LT p_{T}^{reco}/p_{T}^{gen} #GT" );

    for( int ieta = 0; ieta < nEta; ieta++ ){
        TString distName = L2Name::ObjectName( cone, "response",
            { etaMode, L2Name::EtaKey( ieta, fullEta ) }, { collection } );

        TH1D* hProj = ProjectTHnSparse1D( h, kRespRAxis,
            { { kRespEtaGenAxis, edges[ieta], edges[ieta + 1] } } );
        hProj->SetName( distName );

        dQA->cd();
        hProj->Write();

        ResponseFitResult fr = FitResponse( hProj, halfWidth, minEntries );
        if( fr.valid ){
            hJES->SetBinContent( ieta + 1, fr.jes );
            hJES->SetBinError( ieta + 1, fr.jesErr );
            hJER->SetBinContent( ieta + 1, fr.jer );
            hJER->SetBinError( ieta + 1, fr.jerErr );
        }
        delete hProj;
    }

    dOut->cd();
    hJES->Write();
    hJER->Write();
    delete hJES;
    delete hJER;
}

// ---- vs pt_gen: marginal over both eta axes and pt_reco entirely, binned
// directly off the response sparse's own uniform pT axis granularity (no
// separate coarse slicing scheme, unlike Step 2's ptavgSlices). ----

static void ExtractVsPtGen(
    THnSparse* h, const TString& cone, const TString& collection,
    const AxisBins& ptBins, TDirectory* dQA, TDirectory* dOut,
    double halfWidth, int minEntries ){

    const int nPt = ptBins.nBins;
    const double lo = ptBins.lo;
    const double hi = ptBins.hi;
    const double width = ( hi - lo ) / nPt;

    TString jesName = L2Name::ObjectName( cone, "JES", { "vs_ptgen" }, { collection } );
    TString jerName = L2Name::ObjectName( cone, "JER", { "vs_ptgen" }, { collection } );
    TH1D* hJES = new TH1D( jesName, "", nPt, lo, hi );
    TH1D* hJER = new TH1D( jerName, "", nPt, lo, hi );
    hJES->GetXaxis()->SetTitle( "p_{T}^{gen} [GeV]" );
    hJES->GetYaxis()->SetTitle( "JES = #LT p_{T}^{reco}/p_{T}^{gen} #GT" );
    hJER->GetXaxis()->SetTitle( hJES->GetXaxis()->GetTitle() );
    hJER->GetYaxis()->SetTitle( "JER = #sigma / #LT p_{T}^{reco}/p_{T}^{gen} #GT" );

    for( int ip = 0; ip < nPt; ip++ ){
        const double ptLo = lo + ip * width;
        const double ptHi = ptLo + width;

        TString distName = L2Name::ObjectName( cone, "response",
            { L2Name::PtGenKey( ptLo, ptHi ) }, { collection } );

        TH1D* hProj = ProjectTHnSparse1D( h, kRespRAxis,
            { { kRespPtGenAxis, ptLo, ptHi } } );
        hProj->SetName( distName );

        dQA->cd();
        hProj->Write();

        ResponseFitResult fr = FitResponse( hProj, halfWidth, minEntries );
        if( fr.valid ){
            hJES->SetBinContent( ip + 1, fr.jes );
            hJES->SetBinError( ip + 1, fr.jesErr );
            hJER->SetBinContent( ip + 1, fr.jer );
            hJER->SetBinError( ip + 1, fr.jerErr );
        }
        delete hProj;
    }

    dOut->cd();
    hJES->Write();
    hJER->Write();
    delete hJES;
    delete hJER;
}

// ============================================================

void runResponse( TString inputFile, TString outputFile ){

    const AnalysisConfig& cfg = Config();
    PrintConfigSummary( cfg );

    const double halfWidth = cfg.responseGausFitHalfWidth;
    const int minEntries = ( cfg.minEntriesPerBin > 0 ) ? cfg.minEntriesPerBin : 100;

    TFile* fIn = TFile::Open( inputFile, "read" );
    if( !fIn || fIn->IsZombie() ){ std::cerr << "Cannot open " << inputFile << "\n"; return; }

    // sentinel pre-scan: incl/tag/probe response sparses are always created
    // together by ConeHistograms::Init when isMC (even if some end up
    // unfilled), so checking incl_resp alone is enough to know whether a
    // cone has response objects at all -- mirrors Kinematics.h's use of
    // "_incl" as its own Step-1-file sentinel.
    int totalSteps = 0;
    for( const TString& cone : cfg.coneLabels ){
        if( fIn->Get( cone + "/" + cone + "_incl_resp" ) ) totalSteps += kNCollections * 3;
    }
    if( totalSteps == 0 ){
        std::cerr << "No MC response histograms found in " << inputFile
                  << " -- is this a Step 1 file run with MODE=mc?\n";
        fIn->Close();
        return;
    }

    { Ssiz_t sl = outputFile.Last( '/' ); if( sl != kNPOS ){ gSystem->mkdir( TString( outputFile( 0, sl ) ), kTRUE ); } }
    TFile* fOut = new TFile( outputFile, "recreate" );

    BinningConfig bins;
    ProgressBar pb( "Extracting response:", totalSteps );

    for( const TString& cone : cfg.coneLabels ){

        TDirectory* coneDir = ( TDirectory* )fIn->Get( cone );
        THnSparse* hSentinel = coneDir ? ( THnSparse* )coneDir->Get( cone + "_incl_resp" ) : nullptr;
        if( !hSentinel ) continue;

        TDirectory* coneDirOut = fOut->mkdir( cone.Data() );
        TDirectory* dQA_abseta = coneDirOut->mkdir( "QA_response_abseta" );
        TDirectory* dQA_fulleta = coneDirOut->mkdir( "QA_response_fulleta" );
        TDirectory* dQA_ptgen = coneDirOut->mkdir( "QA_response_ptgen" );

        for( int ic = 0; ic < kNCollections; ic++ ){
            const TString collection = kCollectionKeys[ic];
            THnSparse* hRaw = ( THnSparse* )coneDir->Get( cone + kCollectionSuffixes[ic] );
            if( !hRaw ){
                std::cerr << "Missing " << cone << kCollectionSuffixes[ic] << " in " << inputFile << "\n";
                continue;
            }

            THnSparse* hFolded = FoldEtaAxis( hRaw, kRespEtaGenAxis,
                cone + kCollectionSuffixes[ic] + "_abseta" );

            ExtractVsEtaGen( hFolded, cone, collection, false, dQA_abseta, coneDirOut, halfWidth, minEntries );
            pb.Update();
            ExtractVsEtaGen( hRaw, cone, collection, true, dQA_fulleta, coneDirOut, halfWidth, minEntries );
            pb.Update();
            ExtractVsPtGen( hRaw, cone, collection, bins.pt, dQA_ptgen, coneDirOut, halfWidth, minEntries );
            pb.Update();

            delete hFolded;
        }
    }

    pb.Finish();
    fOut->Close();
    fIn->Close();
}
