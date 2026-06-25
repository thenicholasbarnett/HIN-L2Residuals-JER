#ifndef DIJETHISTOGRAMS_H
#define DIJETHISTOGRAMS_H

#include "TH1D.h"
#include "THnSparse.h"
#include "TString.h"

#include "Binning.h"
#include "Utilities.h"
#include "Dijet.h"

// One instance per cone size. Each instance owns a single 4D THnSparse filled
// once per valid dijet event at the event's actual (etaProbe, pT_avg, alpha, A).
// Alpha slices are extracted cumulatively at analysis time via SetRangeUser —
// see BinningConfig::alphaSlices.

struct ConeHistograms {

    // axis indices into hAsym — used by extraction code for SetRangeUser and Projection
    static constexpr int kEtaAxis   = 0;
    static constexpr int kPtAvgAxis = 1;
    static constexpr int kAlphaAxis = 2;
    static constexpr int kAAxis     = 3;

    // main physics histogram: 4D (eta_probe, pT_avg, alpha, A)
    THnSparse* hAsym = nullptr;

    // control histograms — inclusive jets (all corrected jets in event), tag, probe
    TH1D* hPtIncl    = nullptr;
    TH1D* hEtaIncl   = nullptr;
    TH1D* hPhiIncl   = nullptr;
    TH1D* hPtTag     = nullptr;
    TH1D* hEtaTag    = nullptr;
    TH1D* hPhiTag    = nullptr;
    TH1D* hPtProbe   = nullptr;
    TH1D* hEtaProbe  = nullptr;
    TH1D* hPhiProbe  = nullptr;

    // prefix is the cone size label, e.g. "ak4PF" — prepended to all histogram names
    void Init(const TString& prefix, const BinningConfig& bins) {
        hAsym = MakeTHnSparse<THnSparseD>(prefix + "_asym", "",
            {bins.eta, bins.ptavg, bins.alpha, bins.asymmetry});
        SetEtaBins(hAsym, kEtaAxis);
        hAsym->Sumw2();

        hPtIncl   = MakeTH1<TH1D>(prefix + "_pt_incl",   bins.pt);
        hEtaIncl  = MakeTH1<TH1D>(prefix + "_eta_incl",  bins.eta);
        hPhiIncl  = MakeTH1<TH1D>(prefix + "_phi_incl",  bins.phi);
        hPtTag    = MakeTH1<TH1D>(prefix + "_pt_tag",    bins.pt);
        hEtaTag   = MakeTH1<TH1D>(prefix + "_eta_tag",   bins.eta);
        hPhiTag   = MakeTH1<TH1D>(prefix + "_phi_tag",   bins.phi);
        hPtProbe  = MakeTH1<TH1D>(prefix + "_pt_probe",  bins.pt);
        hEtaProbe = MakeTH1<TH1D>(prefix + "_eta_probe", bins.eta);
        hPhiProbe = MakeTH1<TH1D>(prefix + "_phi_probe", bins.phi);
    }

    // called once per valid dijet event
    void Fill(const DijetResult& d, const float* pt, const float* eta, const float* phi, float weight) {
        double x[4] = {eta[d.probeIdx], d.ptavg, d.alpha, d.A};
        hAsym->Fill(x, weight);

        hPtTag   ->Fill(pt [d.tagIdx],   weight);
        hEtaTag  ->Fill(eta[d.tagIdx],   weight);
        hPhiTag  ->Fill(phi[d.tagIdx],   weight);
        hPtProbe ->Fill(pt [d.probeIdx], weight);
        hEtaProbe->Fill(eta[d.probeIdx], weight);
        hPhiProbe->Fill(phi[d.probeIdx], weight);
    }

    // called per corrected jet in the inclusive jet loop
    void FillInclJet(float corrPt, float jetEta, float jetPhi, float weight) {
        hPtIncl ->Fill(corrPt,  weight);
        hEtaIncl->Fill(jetEta,  weight);
        hPhiIncl->Fill(jetPhi,  weight);
    }

    void Write() {
        WriteAll(hAsym);
        WriteAll(hPtIncl);   WriteAll(hEtaIncl);   WriteAll(hPhiIncl);
        WriteAll(hPtTag);    WriteAll(hEtaTag);     WriteAll(hPhiTag);
        WriteAll(hPtProbe);  WriteAll(hEtaProbe);   WriteAll(hPhiProbe);
    }
};

#endif
