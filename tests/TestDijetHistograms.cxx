#include <iostream>
#include <cmath>

#include "TROOT.h"
#include "THnSparse.h"

#include "Binning.h"
#include "DijetHistograms.h"
#include "Dijet.h"

static int nPass = 0;
static int nFail = 0;
static constexpr double kEps = 1e-6;

void Check( bool cond, const char* msg ){
    if( cond ){ std::cout << "  PASS  " << msg << std::endl; nPass++; }
    else { std::cout << "  FAIL  " << msg << std::endl; nFail++; }
}

int main(){
    gROOT->SetBatch( true );

    std::cout << "=== TestDijetHistograms ===" << std::endl;

    BinningConfig bins;

    // ── data mode: no response histograms ───────────────────────────────────
    std::cout << "\n[1] Data mode: no response histograms" << std::endl;
    {
        ConeHistograms h;
        h.Init( "test_data", bins, false );
        Check( h.hInclJetResp == nullptr, "hInclJetResp null in data mode" );
        Check( h.hTagJetResp == nullptr, "hTagJetResp null in data mode" );
        Check( h.hProbeJetResp == nullptr, "hProbeJetResp null in data mode" );
    }

    // ── MC mode: response histograms are created ────────────────────────────
    std::cout << "\n[2] MC mode: response histograms created" << std::endl;
    ConeHistograms h;
    h.Init( "test_mc", bins, true );
    Check( h.hInclJetResp != nullptr, "hInclJetResp created in MC mode" );
    Check( h.hTagJetResp != nullptr, "hTagJetResp created in MC mode" );
    Check( h.hProbeJetResp != nullptr, "hProbeJetResp created in MC mode" );

    // ── FillInclJetResp: matched jet lands at pT/refPt ──────────────────────
    // Values are chosen away from axis bin edges (phi bins are 0.1-wide, pT
    // bins are 10-wide, response bins are 0.01-wide) and exactly representable
    // in float, so there's no float->double rounding ambiguity at a boundary
    // between this fill and the test's independent FindBin lookup below.
    std::cout << "\n[3] FillInclJetResp: matched jet" << std::endl;
    {
        // corrPt=105, eta=0.13, phi=0.25, refPt=80 -> response = 1.3125
        h.FillInclJetResp( 105.0f, 0.13f, 0.25f, 80.0f, 1.0f );
        Double_t x[4] = { 0.13, 0.25, 105.0, 1.3125 };
        Int_t coords[4];
        for( Int_t i = 0; i < 4; i++ ) coords[i] = h.hInclJetResp->GetAxis( i )->FindBin( x[i] );
        Double_t content = h.hInclJetResp->GetBinContent( coords );
        Check( std::fabs( content - 1.0 ) < kEps, "matched jet fills response = corrPt/refPt = 1.3125" );
    }

    // ── FillInclJetResp: unmatched jet (refPt < 0) is dropped ───────────────
    std::cout << "\n[4] FillInclJetResp: unmatched jet skipped" << std::endl;
    {
        Long64_t before = h.hInclJetResp->GetNbins();
        h.FillInclJetResp( 50.0f, 1.0f, 1.0f, -1.0f, 1.0f );
        Long64_t after = h.hInclJetResp->GetNbins();
        Check( before == after, "unmatched jet (refPt<0) does not add a filled bin" );
    }

    // ── FillResp: tag matched, probe unmatched — independent per-leg ────────
    std::cout << "\n[5] FillResp: tag matched, probe unmatched" << std::endl;
    {
        float pt[]    = { 117.0f, 80.0f };
        float eta[]   = { 0.3f, 2.5f };
        float phi[]   = { 0.05f, 3.0f };
        float refPt[] = { 96.0f, -1.0f };   // tag matched, probe unmatched

        DijetResult d;
        d.tagIdx = 0;
        d.probeIdx = 1;

        Long64_t beforeTag = h.hTagJetResp->GetNbins();
        Long64_t beforeProbe = h.hProbeJetResp->GetNbins();

        h.FillResp( d, pt, eta, phi, refPt, 1.0f );

        Long64_t afterTag = h.hTagJetResp->GetNbins();
        Long64_t afterProbe = h.hProbeJetResp->GetNbins();

        Check( afterTag == beforeTag + 1, "tag jet (matched) fills hTagJetResp" );
        Check( afterProbe == beforeProbe, "probe jet (unmatched) leaves hProbeJetResp untouched" );

        Double_t x[4] = { eta[0], phi[0], pt[0], pt[0] / refPt[0] };
        Int_t coords[4];
        for( Int_t i = 0; i < 4; i++ ) coords[i] = h.hTagJetResp->GetAxis( i )->FindBin( x[i] );
        Double_t content = h.hTagJetResp->GetBinContent( coords );
        Check( std::fabs( content - 1.0 ) < kEps, "tag response value = pT/refPt = 1.21875" );
    }

    std::cout << "\n=== " << nPass << " passed, " << nFail << " failed ===" << std::endl;
    return nFail > 0 ? 1 : 0;
}
