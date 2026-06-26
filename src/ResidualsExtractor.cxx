#include "ResidualsExtractor.h"

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "THnSparse.h"
#include "TDirectory.h"
#include "TString.h"
#include "TMath.h"

#include "Binning.h"
#include "Utilities.h"
#include "2024ppRef.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

// ---- constants ----

static constexpr int    kMinEntries   = 100;
static constexpr int    kNAlphaFit    = 6;     // alpha thresholds 0.05–0.30 only
static constexpr double kGausFitHW    = 0.5;
static constexpr double kMaxAbsA_fit  = 0.9;
static constexpr double kAlphaFitHi   = 0.31;

// ---- axis indices — same before and after FoldEtaAxis ----
static constexpr int kEtaAxis   = 0;
static constexpr int kPtAvgAxis = 1;
static constexpr int kAlphaAxis = 2;
static constexpr int kAAxis     = 3;

// ---- methods ----
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

static GaussResult FitGauss(TH1D* h) {
    GaussResult r;
    if (!h || h->GetEntries() < kMinEntries) return r;

    TF1* g = new TF1(Form("_gf_%s", h->GetName()), "gaus", -kGausFitHW, kGausFitHW);
    g->SetParameter(0, h->GetMaximum());
    g->SetParameter(1, h->GetMean());
    g->SetParameter(2, std::max(h->GetRMS(), 1e-3));

    TFitResultPtr res = h->Fit(g, "NQSR");
    if (res.Get() && res->IsValid()) {
        r.mean    = res->Parameter(1);
        r.meanErr = res->ParError(1);
        r.sigma   = res->Parameter(2);
        r.chi2ndf = (res->Ndf() > 0) ? res->Chi2() / res->Ndf() : -1.0;
        r.valid   = true;
    }
    delete g;
    return r;
}

// ---- truncated mean: trim equal tails so central `fraction` of area survives ----

static TruncResult TruncMean(TH1D* h, double fraction) {
    TruncResult r;
    if (!h || h->GetEntries() < kMinEntries) return r;

    double total = h->Integral();
    if (total <= 0) return r;
    double tailN = 0.5 * (1.0 - fraction) * total;

    int nBins = h->GetNbinsX();

    double cum = 0;
    int binLo = 1;
    for (int b = 1; b <= nBins; b++) {
        cum += h->GetBinContent(b);
        if (cum > tailN) { binLo = b; break; }
    }

    cum = 0;
    int binHi = nBins;
    for (int b = nBins; b >= 1; b--) {
        cum += h->GetBinContent(b);
        if (cum > tailN) { binHi = b; break; }
    }

    if (binLo > binHi) return r;

    double nEff = h->Integral(binLo, binHi);
    if (nEff < 10) return r;

    h->GetXaxis()->SetRange(binLo, binHi);
    double mean = h->GetMean();
    double rms  = h->GetRMS();
    h->GetXaxis()->SetRange(0, 0);

    r.mean    = mean;
    r.meanErr = rms / TMath::Sqrt(nEff);
    r.nEff    = nEff;
    r.valid   = true;
    return r;
}

// ---- <A> → R and propagated error ----

static double ToR(double A)               { return (1.0 + A) / (1.0 - A); }
static double ToRErr(double A, double dA) { return dA * 2.0 / ((1.0 - A) * (1.0 - A)); }

static void ResetRange(THnSparse* h, int axis) { h->GetAxis(axis)->SetRange(0, 0); }

// ============================================================
// ExtractAndFit — runs the full extraction loop for one set of sparses.
//
// etaEdges:   kAbsEtaEdges (18 bins) or kEtaEdges (36 bins)
// nameSuffix: "" for |eta|, "_fulleta" for full eta
//             appended to graph names and intercept TH1D names
// ============================================================

static void ExtractAndFit(
    THnSparse* hData, THnSparse* hMC,
    const TString& cone,
    const BinningConfig& bins,
    TDirectory* dQA_data, TDirectory* dQA_mc,
    TDirectory* dGraphs,
    TFile* fOut,
    const std::vector<Double_t>& etaEdges,
    const TString& nameSuffix)
{
    const int nPt    = (int)bins.ptavgSlices.size();
    const int nAlpha = (int)bins.alphaSlices.size();
    const int nEta   = hData->GetAxis(kEtaAxis)->GetNbins();

    struct RPoint { double alpha, val, err; };
    std::vector<std::vector<std::vector<std::vector<RPoint>>>>
        rpts(kNMethods,
            std::vector<std::vector<std::vector<RPoint>>>(nPt,
                std::vector<std::vector<RPoint>>(nEta)));

    for (int ipt = 0; ipt < nPt; ipt++) {
        const auto& ptSlice = bins.ptavgSlices[ipt];
        std::cout << "    " << ptSlice.title << "\n";

        for (int ialpha = 0; ialpha < nAlpha; ialpha++) {
            const auto& aSlice = bins.alphaSlices[ialpha];

            hData->GetAxis(kPtAvgAxis)->SetRangeUser(ptSlice.lo, ptSlice.hi);
            hMC  ->GetAxis(kPtAvgAxis)->SetRangeUser(ptSlice.lo, ptSlice.hi);
            hData->GetAxis(kAlphaAxis)->SetRangeUser(aSlice.lo,  aSlice.hi);
            hMC  ->GetAxis(kAlphaAxis)->SetRangeUser(aSlice.lo,  aSlice.hi);

            for (int ieta = 0; ieta < nEta; ieta++) {
                int etaBin = ieta + 1;

                hData->GetAxis(kEtaAxis)->SetRange(etaBin, etaBin);
                hMC  ->GetAxis(kEtaAxis)->SetRange(etaBin, etaBin);

                TH1D* hAData = (TH1D*)hData->Projection(kAAxis);
                TH1D* hAMC   = (TH1D*)hMC  ->Projection(kAAxis);

                TString suffix = Form("%s%s_eta%02d",
                    ptSlice.shortName.Data(), aSlice.shortName.Data(), ieta);
                hAData->SetName(cone + nameSuffix + "_A_data_" + suffix);
                hAMC  ->SetName(cone + nameSuffix + "_A_mc_"   + suffix);

                dQA_data->cd(); hAData->Write();
                dQA_mc  ->cd(); hAMC  ->Write();

                GaussResult gd   = FitGauss(hAData);
                GaussResult gm   = FitGauss(hAMC);
                TruncResult td90 = TruncMean(hAData, 0.90);
                TruncResult tm90 = TruncMean(hAMC,   0.90);
                TruncResult td95 = TruncMean(hAData, 0.95);
                TruncResult tm95 = TruncMean(hAMC,   0.95);

                if (ialpha < kNAlphaFit) {
                    double alphaX = aSlice.hi;

                    auto accum = [&](int method, double Ad, double eAd,
                                                 double Am, double eAm, bool ok) {
                        if (!ok) return;
                        if (std::abs(Ad) > kMaxAbsA_fit || std::abs(Am) > kMaxAbsA_fit) return;
                        double Rd = ToR(Ad), Rm = ToR(Am);
                        double eRd = ToRErr(Ad, eAd), eRm = ToRErr(Am, eAm);
                        if (std::abs(Rm) < 1e-6) return;
                        double ratio  = Rd / Rm;
                        double eRatio = ratio * TMath::Sqrt(
                            (eRd/Rd)*(eRd/Rd) + (eRm/Rm)*(eRm/Rm));
                        rpts[method][ipt][ieta].push_back({alphaX, ratio, eRatio});
                    };

                    accum(0, gd.mean,   gd.meanErr,   gm.mean,   gm.meanErr,   gd.valid  && gm.valid);
                    accum(1, td90.mean, td90.meanErr, tm90.mean, tm90.meanErr, td90.valid && tm90.valid);
                    accum(2, td95.mean, td95.meanErr, tm95.mean, tm95.meanErr, td95.valid && tm95.valid);
                }

                ResetRange(hData, kEtaAxis);
                ResetRange(hMC,   kEtaAxis);
                delete hAData;
                delete hAMC;
            }

            ResetRange(hData, kPtAvgAxis);
            ResetRange(hMC,   kPtAvgAxis);
            ResetRange(hData, kAlphaAxis);
            ResetRange(hMC,   kAlphaAxis);
        }
    }

    // ---- build TGraphErrors and fit R_ratio vs alpha ----

    for (int method = 0; method < kNMethods; method++) {
        for (int ipt = 0; ipt < nPt; ipt++) {
            const auto& ptSlice = bins.ptavgSlices[ipt];

            TString corrName = Form("%s_intercept_%s%s%s",
                cone.Data(), kMethodNames[method],
                ptSlice.shortName.Data(), nameSuffix.Data());
            TH1D* hCorr = new TH1D(corrName, "",
                (int)etaEdges.size() - 1, etaEdges.data());
            hCorr->GetXaxis()->SetTitle(nameSuffix.IsNull() ? "|#eta|" : "#eta");
            hCorr->GetYaxis()->SetTitle("R_{data}/R_{MC} at #alpha=0");

            for (int ieta = 0; ieta < nEta; ieta++) {
                const auto& pts = rpts[method][ipt][ieta];
                if ((int)pts.size() < 2) continue;

                int n = (int)pts.size();
                std::vector<double> x(n), y(n), ex(n, 0.0), ey(n);
                for (int k = 0; k < n; k++) {
                    x[k] = pts[k].alpha; y[k] = pts[k].val; ey[k] = pts[k].err;
                }

                TString gname = Form("%s_R_%s%s_eta%02d%s",
                    cone.Data(), kMethodNames[method],
                    ptSlice.shortName.Data(), ieta, nameSuffix.Data());

                TGraphErrors* gr = new TGraphErrors(n,
                    x.data(), y.data(), ex.data(), ey.data());
                gr->SetName(gname);
                gr->SetTitle(";#alpha threshold;R_{data}/R_{MC}");
                gr->SetMarkerStyle(20);
                gr->SetMarkerColor(ptSlice.color);
                gr->SetLineColor(ptSlice.color);

                TF1* fitFn = new TF1(gname + "_fit", "[0]+[1]*x", 0.0, kAlphaFitHi);
                fitFn->SetParameter(0, 1.0);
                fitFn->SetParameter(1, 0.0);
                fitFn->SetLineColor(ptSlice.color);
                gr->Fit(fitFn, "Q");

                hCorr->SetBinContent(ieta + 1, fitFn->GetParameter(0));
                hCorr->SetBinError  (ieta + 1, fitFn->GetParError(0));

                dGraphs->cd();
                gr->Write();
                gr->GetListOfFunctions()->Remove(fitFn);
                delete fitFn;
                delete gr;
            }

            fOut->cd();
            hCorr->Write();
            delete hCorr;
        }
    }
}

// ============================================================

void runResiduals(TString dataFile, TString mcFile, TString outputFile) {

    std::cout << "Data:   " << dataFile   << "\n"
              << "MC:     " << mcFile     << "\n"
              << "Output: " << outputFile << "\n";

    TFile* fData = TFile::Open(dataFile, "read");
    TFile* fMC   = TFile::Open(mcFile,   "read");
    if (!fData || fData->IsZombie()) { std::cerr << "Cannot open " << dataFile << "\n"; return; }
    if (!fMC   || fMC->IsZombie())   { std::cerr << "Cannot open " << mcFile   << "\n"; return; }

    TFile* fOut = new TFile(outputFile, "recreate");

    BinningConfig bins;

    for (const TString& cone : kConeLabels) {

        THnSparse* hRawData = (THnSparse*)fData->Get(cone + "_asym");
        THnSparse* hRawMC   = (THnSparse*)fMC  ->Get(cone + "_asym");
        if (!hRawData) { std::cerr << "Missing " << cone << "_asym in data\n"; continue; }
        if (!hRawMC)   { std::cerr << "Missing " << cone << "_asym in MC\n";   continue; }

        std::cout << "\n=== " << cone << " ===\n";

        // fold full-eta → |eta| once; both folded sparses are written to fOut
        THnSparse* hData = FoldEtaAxis(hRawData, kEtaAxis, cone + "_asym_data_abseta");
        THnSparse* hMC   = FoldEtaAxis(hRawMC,   kEtaAxis, cone + "_asym_mc_abseta");
        fOut->cd();
        hData->Write();
        hMC  ->Write();

        TDirectory* dQA_data      = fOut->mkdir(cone + "_QA_data");
        TDirectory* dQA_mc        = fOut->mkdir(cone + "_QA_mc");
        TDirectory* dQA_data_full = fOut->mkdir(cone + "_QA_data_fulleta");
        TDirectory* dQA_mc_full   = fOut->mkdir(cone + "_QA_mc_fulleta");
        TDirectory* dGraphs       = fOut->mkdir(cone + "_graphs");

        std::cout << "  |eta| extraction\n";
        ExtractAndFit(hData, hMC, cone, bins,
                      dQA_data, dQA_mc, dGraphs, fOut,
                      kAbsEtaEdges, "");

        std::cout << "  full eta extraction\n";
        ExtractAndFit(hRawData, hRawMC, cone, bins,
                      dQA_data_full, dQA_mc_full, dGraphs, fOut,
                      kEtaEdges, "_fulleta");

        delete hData;
        delete hMC;
    }

    fOut->Close();
    fData->Close();
    fMC  ->Close();

    std::cout << "\nDone. Output: " << outputFile << "\n";
}
