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
//   hInclJet — TH3D (eta, phi, pT) for all corrected jets passing kMinPt
//   hTagJet  — TH3D (eta, phi, pT) for the tag jet of each valid dijet
//   hProbeJet — TH3D (eta, phi, pT) for the probe jet of each valid dijet
//
// The TH3Ds use variable-width CMS JEC eta bins on X, enabling eta-phi maps
// at arbitrary pT thresholds via SetRangeUser on Z + Project3D("yx").

struct ConeHistograms {

    static constexpr int kEtaAxis   = 0;
    static constexpr int kPtAvgAxis = 1;
    static constexpr int kAlphaAxis = 2;
    static constexpr int kAAxis     = 3;

    THnSparse* hAsym     = nullptr;
    TH3D*      hInclJet  = nullptr;
    TH3D*      hTagJet   = nullptr;
    TH3D*      hProbeJet = nullptr;

    void Init(const TString& prefix, const BinningConfig& bins) {
        hAsym = MakeTHnSparse<THnSparseD>(prefix + "_asym", "",
            {bins.eta, bins.ptavg, bins.alpha, bins.asymmetry});
        SetEtaBins(hAsym, kEtaAxis);
        hAsym->Sumw2();

        hInclJet  = MakeTH3DEtaPhiPt(prefix + "_incl",  kEtaEdges, bins.phi, bins.pt);
        hTagJet   = MakeTH3DEtaPhiPt(prefix + "_tag",   kEtaEdges, bins.phi, bins.pt);
        hProbeJet = MakeTH3DEtaPhiPt(prefix + "_probe", kEtaEdges, bins.phi, bins.pt);
    }

    void Fill(const DijetResult& d, const float* pt, const float* eta,
              const float* phi, float weight) {
        double x[4] = {eta[d.probeIdx], d.ptavg, d.alpha, d.A};
        hAsym->Fill(x, weight);

        hTagJet  ->Fill(eta[d.tagIdx],   phi[d.tagIdx],   pt[d.tagIdx],   weight);
        hProbeJet->Fill(eta[d.probeIdx], phi[d.probeIdx], pt[d.probeIdx], weight);
    }

    void FillInclJet(float corrPt, float jetEta, float jetPhi, float weight) {
        hInclJet->Fill(jetEta, jetPhi, corrPt, weight);
    }

    void Write() {
        WriteAll(hAsym);
        WriteAll(hInclJet);
        WriteAll(hTagJet);
        WriteAll(hProbeJet);
    }
};

#endif
