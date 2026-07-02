#ifndef DIJETHISTOGRAMS_H
#define DIJETHISTOGRAMS_H

#include "TH3.h"
#include "THnSparse.h"
#include "TString.h"

#include "Binning.h"
#include "Utilities.h"
#include "Dijet.h"

// One instance per cone size. Each instance owns:
//   hAsym    — 4D THnSparse (eta_probe, pT_avg, alpha, A), one fill per valid dijet event
//   hInclJet — TH3D (eta, phi, pT) for all corrected jets passing cfg.minPt
//   hTagJet  — TH3D (eta, phi, pT) for the tag jet of each valid dijet
//   hProbeJet — TH3D (eta, phi, pT) for the probe jet of each valid dijet
//
// The TH3Ds use variable-width CMS JEC eta bins on X, enabling eta-phi maps
// at arbitrary pT thresholds via SetRangeUser on Z + Project3D("yx").
//
// MC-only (Init's isMC flag), JES/JER inputs — parallel to, not a replacement
// for, the TH3Ds above (those keep the same behavior/stats for both data and
// MC):
//   hInclJetResp  — 4D THnSparse (eta, phi, pT, pT^reco/pT^gen) for corrected
//                    jets passing cfg.minPt, ref-matched only
//   hTagJetResp   — same, for the tag jet of each valid dijet, ref-matched only
//   hProbeJetResp — same, for the probe jet of each valid dijet, ref-matched only
// "ref-matched" means the ntuple's refpt[j] (index-aligned to reco jet j) is
// >= 0; a negative refpt marks no matching GenJet, and that jet is dropped
// from these histograms entirely (never filled with a bogus ratio).

struct ConeHistograms {

    static constexpr int kEtaAxis = 0;
    static constexpr int kPtAvgAxis = 1;
    static constexpr int kAlphaAxis = 2;
    static constexpr int kAAxis = 3;

    // axis order within hInclJetResp/hTagJetResp/hProbeJetResp
    static constexpr int kRespEtaAxis = 0;
    static constexpr int kRespPhiAxis = 1;
    static constexpr int kRespPtAxis = 2;
    static constexpr int kRespRAxis = 3;

    THnSparse* hAsym = nullptr;
    TH3D* hInclJet = nullptr;
    TH3D* hTagJet = nullptr;
    TH3D* hProbeJet = nullptr;

    THnSparse* hInclJetResp = nullptr;
    THnSparse* hTagJetResp = nullptr;
    THnSparse* hProbeJetResp = nullptr;

    void Init( const TString& prefix, const BinningConfig& bins, bool isMC = false ){
        hAsym = MakeTHnSparse<THnSparseD>( prefix + "_asym", "",
            { bins.eta, bins.ptavg, bins.alpha, bins.asymmetry } );
        SetEtaBins( hAsym, kEtaAxis );
        hAsym->Sumw2();

        hInclJet = MakeTH3DEtaPhiPt( prefix + "_incl", kEtaEdges, bins.phi, bins.pt );
        hTagJet = MakeTH3DEtaPhiPt( prefix + "_tag", kEtaEdges, bins.phi, bins.pt );
        hProbeJet = MakeTH3DEtaPhiPt( prefix + "_probe", kEtaEdges, bins.phi, bins.pt );

        if( isMC ){
            hInclJetResp = MakeRespSparse( prefix + "_incl_resp", bins );
            hTagJetResp = MakeRespSparse( prefix + "_tag_resp", bins );
            hProbeJetResp = MakeRespSparse( prefix + "_probe_resp", bins );
        }
    }

    void Fill( const DijetResult& d, const float* pt, const float* eta,
              const float* phi, float weight ){
        double x[4] = { eta[d.probeIdx], d.ptavg, d.alpha, d.A };
        hAsym->Fill( x, weight );

        hTagJet ->Fill( eta[d.tagIdx], phi[d.tagIdx], pt[d.tagIdx], weight );
        hProbeJet->Fill( eta[d.probeIdx], phi[d.probeIdx], pt[d.probeIdx], weight );
    }

    void FillInclJet( float corrPt, float jetEta, float jetPhi, float weight ){
        hInclJet->Fill( jetEta, jetPhi, corrPt, weight );
    }

    // MC only — refPt is the ntuple's refpt[j] for this jet; refPt < 0 means
    // no matched GenJet, so the jet is dropped rather than filled with a
    // meaningless ratio.
    void FillInclJetResp( float corrPt, float jetEta, float jetPhi, float refPt, float weight ){
        if( refPt < 0 ) return;
        double x[4] = { jetEta, jetPhi, corrPt, corrPt / refPt };
        hInclJetResp->Fill( x, weight );
    }

    // MC only — refPt is index-aligned to pt/eta/phi (the ntuple's refpt[]
    // array). Tag and probe are matched independently: one being unmatched
    // does not drop the other.
    void FillResp( const DijetResult& d, const float* pt, const float* eta,
                  const float* phi, const float* refPt, float weight ){
        if( refPt[d.tagIdx] >= 0 ){
            double xt[4] = { eta[d.tagIdx], phi[d.tagIdx], pt[d.tagIdx], pt[d.tagIdx] / refPt[d.tagIdx] };
            hTagJetResp->Fill( xt, weight );
        }
        if( refPt[d.probeIdx] >= 0 ){
            double xp[4] = { eta[d.probeIdx], phi[d.probeIdx], pt[d.probeIdx], pt[d.probeIdx] / refPt[d.probeIdx] };
            hProbeJetResp->Fill( xp, weight );
        }
    }

    void Write( TDirectory* dir = nullptr ){
        TDirectory* saved = gDirectory;
        if( dir ) dir->cd();
        WriteAll( hAsym );
        WriteAll( hInclJet );
        WriteAll( hTagJet );
        WriteAll( hProbeJet );
        WriteAll( hInclJetResp );
        WriteAll( hTagJetResp );
        WriteAll( hProbeJetResp );
        if( dir && saved ) saved->cd();
    }

private:
    // eta, phi, pT, response — variable-width CMS JEC eta bins on axis 0,
    // matching the TH3Ds' X axis (kEtaEdges).
    static THnSparse* MakeRespSparse( const TString& name, const BinningConfig& bins ){
        THnSparse* h = MakeTHnSparse<THnSparseD>( name, "",
            { bins.eta, bins.phi, bins.pt, bins.response } );
        SetEtaBins( h, kRespEtaAxis );
        h->Sumw2();
        return h;
    }
};

#endif
