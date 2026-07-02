#include <iostream>
#include <cmath>
#include <cstdio>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"

#include "Binning.h"
#include "DijetHistograms.h"
#include "ResponseExtractor.h"

static int nPass = 0;
static int nFail = 0;
static constexpr double kEps = 1e-6;

void Check( bool cond, const char* msg ){
    if( cond ){ std::cout << "  PASS  " << msg << std::endl; nPass++; }
    else { std::cout << "  FAIL  " << msg << std::endl; nFail++; }
}

// Response values chosen symmetric around 1.0 so the expected fitted mean
// (JES) is exactly 1.0 regardless of fit-window/binning details; the spread
// gives FitResponse something real to measure for JER.
static const double kResponseValues[] = { 0.90, 0.93, 0.96, 0.99, 1.00, 1.01, 1.04, 1.07, 1.10 };
static constexpr int kNResponseValues = 9;
static constexpr int kRepeatsPerValue = 20;  // 180 total entries, well above min_entries_per_bin=100

int main(){
    gROOT->SetBatch( true );

    std::cout << "=== TestResponseExtractor ===" << std::endl;

    const char* inputPath = "/tmp/tre_input.root";
    const char* outputPath = "/tmp/tre_output.root";
    std::remove( inputPath );
    std::remove( outputPath );

    // ---- build a synthetic Step 1 MC output file (ak4PF, incl collection only) ----
    {
        BinningConfig bins;
        ConeHistograms h;
        h.Init( "ak4PF", bins, true );

        // eta=refEta=0.13 -> |eta_gen| bin 1 [0,0.261], full-eta bin covering [0,0.261]
        // refPt=105 -> pt_gen bin [100,110)
        const float eta = 0.13f;
        const float refEta = 0.13f;
        const float refPt = 105.0f;
        for( int iv = 0; iv < kNResponseValues; iv++ ){
            const float corrPt = ( float )( kResponseValues[iv] * refPt );
            for( int r = 0; r < kRepeatsPerValue; r++ ){
                h.FillInclJetResp( corrPt, eta, refPt, refEta, 1.0f );
            }
        }

        // a deliberately under-filled pt_gen bin (only 5 entries, well below
        // min_entries_per_bin=100) -- should NOT produce a valid fit.
        for( int r = 0; r < 5; r++ ){
            h.FillInclJetResp( 300.0f, eta, 300.0f, refEta, 1.0f );
        }

        TFile* fo = new TFile( inputPath, "recreate" );
        TDirectory* dir = fo->mkdir( "ak4PF" );
        h.Write( dir );
        fo->Close();
        delete fo;
    }

    // ---- run the real extraction ----
    runResponse( inputPath, outputPath );

    // ---- verify output ----
    TFile* fIn = TFile::Open( outputPath, "read" );
    Check( fIn && !fIn->IsZombie(), "output file opens" );
    if( !fIn || fIn->IsZombie() ){
        std::cout << "\n=== " << nPass << " passed, " << nFail << " failed ===" << std::endl;
        return 1;
    }

    std::cout << "\n[1] vs pt_gen" << std::endl;
    {
        TH1D* hJES = ( TH1D* )fIn->Get( "ak4PF/ak4PF_JES_vs_ptgen_incl" );
        TH1D* hJER = ( TH1D* )fIn->Get( "ak4PF/ak4PF_JER_vs_ptgen_incl" );
        Check( hJES != nullptr, "JES_vs_ptgen_incl exists" );
        Check( hJER != nullptr, "JER_vs_ptgen_incl exists" );
        if( hJES && hJER ){
            const int bin = hJES->GetXaxis()->FindBin( 105.0 );  // pt_gen=105 -> [100,110) bin
            const double jes = hJES->GetBinContent( bin );
            const double jer = hJER->GetBinContent( bin );
            Check( std::fabs( jes - 1.0 ) < 0.02, "JES(pt_gen~105) close to 1.0 (symmetric synthetic response)" );
            Check( jer > 0.01 && jer < 0.20, "JER(pt_gen~105) in a sane positive range" );

            const int emptyBin = hJES->GetXaxis()->FindBin( 305.0 );  // under-filled bin
            Check( std::fabs( hJES->GetBinContent( emptyBin ) ) < kEps,
                "under-filled pt_gen bin (5 entries < min_entries_per_bin) stays at 0, no fit" );
        }
    }

    std::cout << "\n[2] vs eta_gen (abseta and fulleta)" << std::endl;
    {
        TH1D* hJESAbs = ( TH1D* )fIn->Get( "ak4PF/ak4PF_JES_abseta_vs_etagen_incl" );
        TH1D* hJESFull = ( TH1D* )fIn->Get( "ak4PF/ak4PF_JES_fulleta_vs_etagen_incl" );
        Check( hJESAbs != nullptr, "JES_abseta_vs_etagen_incl exists" );
        Check( hJESFull != nullptr, "JES_fulleta_vs_etagen_incl exists" );
        if( hJESAbs ){
            const double jes = hJESAbs->GetBinContent( hJESAbs->GetXaxis()->FindBin( 0.13 ) );
            Check( std::fabs( jes - 1.0 ) < 0.02, "JES(|eta_gen|~0.13) close to 1.0" );
        }
        if( hJESFull ){
            const double jes = hJESFull->GetBinContent( hJESFull->GetXaxis()->FindBin( 0.13 ) );
            Check( std::fabs( jes - 1.0 ) < 0.02, "JES(eta_gen~0.13) close to 1.0" );
        }
    }

    std::cout << "\n[3] tag/probe collections were never filled -> objects exist but bins stay zero" << std::endl;
    {
        // ConeHistograms::Init always constructs all three response sparses
        // when isMC -- unfilled ones are empty, not absent, so extraction
        // still writes JES/JER TH1Ds for them (matching the unity/zero
        // fallback convention used elsewhere, e.g. TextFileWriter.cxx).
        TH1D* hJESTag = ( TH1D* )fIn->Get( "ak4PF/ak4PF_JES_vs_ptgen_tag" );
        Check( hJESTag != nullptr, "JES_vs_ptgen_tag object still written (unfilled, not absent)" );
        if( hJESTag ){
            bool allZero = true;
            for( int b = 1; b <= hJESTag->GetNbinsX(); b++ )
                if( std::fabs( hJESTag->GetBinContent( b ) ) > kEps ) allZero = false;
            Check( allZero, "tag collection (never filled) has no valid fits anywhere -> all bins zero" );
        }
    }

    fIn->Close();
    std::remove( inputPath );
    std::remove( outputPath );

    std::cout << "\n=== " << nPass << " passed, " << nFail << " failed ===" << std::endl;
    return nFail > 0 ? 1 : 0;
}
