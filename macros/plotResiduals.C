// Compiled:    ./bin/plotResiduals <residuals.root> [out_dir]
// Interpreted: root -l -b -q 'macros/plotResiduals.C("residuals.root")'
//              (run from repo root so relative paths resolve)

#ifdef __CLING__
R__ADD_INCLUDE_PATH(include)
R__ADD_INCLUDE_PATH(cfg)
#if defined(__APPLE__)
R__LOAD_LIBRARY(lib/libl2residuals.dylib)
#else
R__LOAD_LIBRARY(lib/libl2residuals.so)
#endif
#endif

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"

#include "Binning.h"
#include "Colors.h"
#include "Naming.h"
#include "Utilities.h"
#include "ProgressBar.h"
#include "2024ppRef.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

// ---- method palette ----

static const char* const kMethodKeys[]   = { "gauss",       "trunc90",          "trunc95"          };
static const char* const kMethodLabels[] = { "Gauss fit",   "Trunc. mean 90%",  "Trunc. mean 95%"  };
static const Color_t     kMethodColors[] = { HiroshigeNightBlue(), HiroshigeOrange(), HiroshigeLightRed() };
static const int         kMethodStyles[] = { 20, 21, 22 };   // circle / square / triangle-up
static constexpr int     kNMethods = 3;

// colors for the eta-symmetry comparison
static const Color_t kColFull = HiroshigeLightRed();
static const Color_t kColRefl = HiroshigeNightBlue();

// ============================================================
// Canvas helpers
// ============================================================

static constexpr Float_t kAspectRatio = 4.0f / 3.0f;
static void RealAspectRatio(TCanvas* c) {
    if (!gROOT->IsBatch()) c->SetRealAspectRatio(kAspectRatio);
}

struct TwoPad {
    TCanvas* c     = nullptr;
    TPad*    main  = nullptr;
    TPad*    ratio = nullptr;
};

// Main panel: y [0.30, 1.0].  Ratio panel: y [0.00, 0.32].
// Pads are added to the canvas's primitive list; delete canvas to cascade-delete both.
static TwoPad MakeTwoPad(const TString& name) {
    TwoPad cv;
    cv.c = new TCanvas(name, "", 800, 600);
    RealAspectRatio(cv.c);
    cv.c->cd();

    cv.main  = new TPad(name + "_m", "", 0.0, 0.30, 1.0, 1.00);
    cv.main->SetBottomMargin(0.025);
    cv.main->SetTopMargin   (0.10);
    cv.main->SetLeftMargin  (0.13);
    cv.main->SetRightMargin (0.04);
    cv.main->Draw();

    cv.ratio = new TPad(name + "_r", "", 0.0, 0.00, 1.0, 0.32);
    cv.ratio->SetTopMargin   (0.05);
    cv.ratio->SetBottomMargin(0.33);
    cv.ratio->SetLeftMargin  (0.13);
    cv.ratio->SetRightMargin (0.04);
    cv.ratio->Draw();

    return cv;
}

// ============================================================
// Histogram helpers
// ============================================================

// Clone h from file, disassociate from directory so the canvas can own it cleanly.
static TH1D* GetH(TFile* f, const TString& name) {
    TH1D* src = (TH1D*)f->Get(name);
    if (!src) return nullptr;
    TH1D* h = (TH1D*)src->Clone(name + "_c");
    h->SetDirectory(0);
    return h;
}

static TH1D* GetHAny(TFile* f, const std::vector<TString>& names) {
    for (const auto& name : names) {
        TH1D* h = GetH(f, name);
        if (h) return h;
    }
    return nullptr;
}

static TH1D* GetHAny(TDirectory* d, const std::vector<TString>& names) {
    if (!d) return nullptr;
    for (const auto& name : names) {
        TH1D* src = (TH1D*)d->Get(name);
        if (!src) continue;
        TH1D* h = (TH1D*)src->Clone(name + "_c");
        h->SetDirectory(0);
        return h;
    }
    return nullptr;
}

static TGraphErrors* GetGraphAny(TDirectory* d, const std::vector<TString>& names) {
    if (!d) return nullptr;
    for (const auto& name : names) {
        TGraphErrors* g = (TGraphErrors*)d->Get(name);
        if (g) return g;
    }
    return nullptr;
}

static TString OldEtaSuffix(bool fullEta) {
    return fullEta ? "_fulleta" : "";
}

// Mirror an |eta| TH1D (kAbsEtaEdges, 18 bins) onto a full-eta TH1D (kEtaEdges, 36 bins).
// Full-eta bin i (1-indexed): i <= 18 → |eta| bin (18-i+1),  i > 18 → |eta| bin (i-18).
static TH1D* Reflect(TH1D* hAbs, const TString& name) {
    const int n = hAbs->GetNbinsX();   // 18
    TH1D* h = new TH1D(name, "", (int)kEtaEdges.size() - 1, kEtaEdges.data());
    h->SetDirectory(0);
    for (int i = 1; i <= 2 * n; i++) {
        const int a = (i <= n) ? (n - i + 1) : (i - n);
        h->SetBinContent(i, hAbs->GetBinContent(a));
        h->SetBinError  (i, hAbs->GetBinError(a));
    }
    return h;
}

// h1 / h2 bin-by-bin; leaves bins zero where either is zero.
static TH1D* RatioH(TH1D* h1, TH1D* h2, const TString& name) {
    TH1D* r = (TH1D*)h1->Clone(name);
    r->SetDirectory(0);
    r->Reset();
    const int n = r->GetNbinsX();
    for (int i = 1; i <= n; i++) {
        const double v1 = h1->GetBinContent(i), e1 = h1->GetBinError(i);
        const double v2 = h2->GetBinContent(i), e2 = h2->GetBinError(i);
        if (std::abs(v2) < 1e-9 || std::abs(v1) < 1e-9) continue;
        const double rv = v1 / v2;
        r->SetBinContent(i, rv);
        r->SetBinError  (i, rv * std::hypot(e1 / v1, e2 / v2));
    }
    return r;
}

static void StyleH(TH1D* h, Color_t col, int mstyle, float lw = 2.0f) {
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    h->SetMarkerStyle(mstyle);
    h->SetMarkerSize(0.85f);
    h->SetLineWidth((Width_t)lw);
}

// Horizontal dashed reference line at y=ref drawn on pad p.
static void RefLine(TPad* p, double x0, double x1, double ref = 1.0) {
    p->cd();
    TLine* l = new TLine(x0, ref, x1, ref);
    l->SetLineStyle(2);
    l->SetLineColor(kGray + 2);
    l->SetLineWidth(1);
    l->Draw();
}

// Tune axis sizes for the ratio pad (30% of canvas height).
// Text sizes are larger to appear the same physical size as the main pad.
static void TuneRatio(TH1D* h, const TString& xTitle, const TString& yTitle,
                      double ylo, double yhi) {
    h->SetTitle("");
    h->GetXaxis()->SetTitle(xTitle);
    h->GetXaxis()->SetTitleSize(0.145);
    h->GetXaxis()->SetTitleOffset(0.85);
    h->GetXaxis()->SetLabelSize(0.120);

    h->GetYaxis()->SetTitle(yTitle);
    h->GetYaxis()->SetTitleSize(0.120);
    h->GetYaxis()->SetTitleOffset(0.44);
    h->GetYaxis()->SetLabelSize(0.100);
    h->GetYaxis()->SetNdivisions(504);
    h->GetYaxis()->SetRangeUser(ylo, yhi);
}

// Auto y-range across a set of histograms, ignoring empty bins, with padding.
static std::pair<double, double> YRange(const std::vector<TH1D*>& hv,
                                        double pad = 0.15) {
    double lo = 1e9, hi = -1e9;
    for (auto* h : hv) {
        if (!h) continue;
        for (int i = 1; i <= h->GetNbinsX(); i++) {
            const double v = h->GetBinContent(i), e = h->GetBinError(i);
            if (v == 0 && e == 0) continue;
            lo = std::min(lo, v - e);
            hi = std::max(hi, v + e);
        }
    }
    if (lo > hi) return {0.88, 1.12};
    const double span = hi - lo;
    return {lo - pad * span, hi + pad * span};
}

static TString PlotDir(const TString& outDir, const TString& cone,
                       const TString& plotType,
                       const std::vector<TString>& orderedKeys = {}) {
    std::vector<TString> parts = {outDir, cone, plotType};
    parts.insert(parts.end(), orderedKeys.begin(), orderedKeys.end());
    TString dir = L2Name::Join(parts, "/");
    if (gSystem->mkdir(dir, true) < 0 && gSystem->AccessPathName(dir))
        std::cerr << "Warning: cannot create plot subdirectory: " << dir << "\n";
    return dir;
}

static void SavePlot(TCanvas* c, const TString& outDir, const TString& cone,
                     const TString& plotType,
                     const std::vector<TString>& orderedKeys,
                     const TString& fileBase) {
    TString dir = PlotDir(outDir, cone, plotType, orderedKeys);
    c->SaveAs(Form("%s/%s.png", dir.Data(), fileBase.Data()));
}

static constexpr int kMinEntriesPlot = 10;

// Return the x-values at the low and high truncation boundaries for `fraction` of the area.
static std::pair<double, double> TruncBounds(TH1D* h, double fraction) {
    if (!h || h->Integral() <= 0) return {0, 0};
    double total = h->Integral();
    double tailN = 0.5 * (1.0 - fraction) * total;
    int nBins = h->GetNbinsX();
    double running = 0; int binLo = 1;
    for (int b = 1; b <= nBins; b++) { running += h->GetBinContent(b); if (running > tailN) { binLo = b; break; } }
    running = 0; int binHi = nBins;
    for (int b = nBins; b >= 1; b--) { running += h->GetBinContent(b); if (running > tailN) { binHi = b; break; } }
    return { h->GetBinLowEdge(binLo), h->GetBinLowEdge(binHi) + h->GetBinWidth(binHi) };
}

static TLine* VLine(double x, Color_t col, int style = 3) {
    // vertical line placeholder — actual y range set by caller after drawing histogram
    TLine* l = new TLine(x, 0, x, 1);
    l->SetLineColor(col);
    l->SetLineStyle(style);
    l->SetLineWidth(2);
    return l;
}

static void NormalizeDensity(TH1D* h) {
    if (!h) return;
    const double norm = h->Integral("width");
    if (norm > 0) h->Scale(1.0 / norm);
}

static double MinPositiveBin(const std::vector<TH1D*>& hv) {
    double minPos = 1e9;
    for (auto* h : hv) {
        if (!h) continue;
        for (int i = 1; i <= h->GetNbinsX(); i++) {
            const double v = h->GetBinContent(i);
            if (v > 0) minPos = std::min(minPos, v);
        }
    }
    return (minPos < 1e9) ? minPos : 1e-6;
}

static void StyleGuideLine(TLine* l, Color_t col) {
    l->SetLineColorAlpha(col, 0.65);
    l->SetLineStyle(3);
    l->SetLineWidth(1);
}

static void StyleFit(TF1* f, Color_t col) {
    f->SetLineColorAlpha(col, 0.70);
    f->SetLineStyle(3);
    f->SetLineWidth(2);
}

static void DrawTruncLines(TH1D* h, double fraction, Color_t col,
                           double yMin, double yMax,
                           TLine*& loLine, TLine*& hiLine) {
    auto [xlo, xhi] = TruncBounds(h, fraction);
    loLine = new TLine(xlo, yMin, xlo, yMax);
    hiLine = new TLine(xhi, yMin, xhi, yMax);
    StyleGuideLine(loLine, col);
    StyleGuideLine(hiLine, col);
    loLine->Draw();
    hiLine->Draw();
}

static TF1* FitGaussianGuide(TH1D* h, const TString& name, Color_t col) {
    if (!h || h->GetEntries() < kMinEntriesPlot) return nullptr;
    TF1* fit = new TF1(name, "gaus", -0.5, 0.5);
    fit->SetParameter(0, h->GetMaximum());
    fit->SetParameter(1, h->GetMean());
    fit->SetParameter(2, std::max(h->GetRMS(), 1e-3));
    StyleFit(fit, col);
    h->Fit(fit, "NQSR");
    return fit;
}

static void DrawAsymBase(TH1D* hData, TH1D* hMC,
                         const TString& xTitle,
                         double yMin, double yMax) {
    StyleH(hData, kBlack, 20, 1.5f);
    StyleH(hMC, kRed + 1, 24, 1.5f);

    hData->SetTitle("");
    hData->GetXaxis()->SetTitle(xTitle);
    hData->GetXaxis()->CenterTitle();
    hData->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{dA}");
    hData->GetYaxis()->CenterTitle();
    hData->SetMinimum(yMin);
    hData->SetMaximum(yMax);

    hData->Draw("E1");
    hMC->Draw("E1 same");
}

static void DrawAsymHeader() {
    TLatex* cms = new TLatex(0.13, 0.905, "#bf{CMS} #it{Internal}");
    cms->SetNDC();
    cms->SetTextFont(42);
    cms->SetTextAlign(11);
    cms->SetTextSize(0.035);
    cms->Draw();

    TLatex* lumi = new TLatex(0.95, 0.905, "2024 pp (5.36 TeV)");
    lumi->SetNDC();
    lumi->SetTextFont(42);
    lumi->SetTextAlign(31);
    lumi->SetTextSize(0.035);
    lumi->Draw();
}

// ============================================================
// Plot type 3: Asymmetry distributions
//
// For each (cone, pT slice, alpha slice, eta bin): three canvases with
//   data (black) and MC (red) overlaid, log-y scale:
//     trunc90/trunc95: data and MC truncation bounds
//     gauss: data and MC Gaussian fit guides
// Skips bins with fewer than kMinEntriesPlot entries.
// ============================================================

static void PlotAsymDist(TFile* fIn, const TString& outDir,
                         const TString& cone, const BinningConfig& bins,
                         ProgressBar& pb) {
    TDirectory* dData = (TDirectory*)fIn->Get(cone + "/QA_data");
    if (!dData) dData = (TDirectory*)fIn->Get(cone + "_QA_data");
    TDirectory* dMC   = (TDirectory*)fIn->Get(cone + "/QA_mc");
    if (!dMC)   dMC   = (TDirectory*)fIn->Get(cone + "_QA_mc");

    const int nPt    = (int)bins.ptavgSlices.size();
    const int nAlpha = (int)bins.alphaSlices.size();
    const int nEta   = (int)kAbsEtaEdges.size() - 1;

    for (int ip = 0; ip < nPt; ip++) {
        const auto& ptSl = bins.ptavgSlices[ip];
        for (int ia = 0; ia < nAlpha; ia++) {
            const auto& aSl = bins.alphaSlices[ia];
            for (int ie = 0; ie < nEta; ie++) {
                TString etaKey = L2Name::EtaKey(ie);
                TString ptKey = L2Name::PtKey(ptSl);
                TString alphaKey = L2Name::AlphaKey(aSl);
                TString oldSfx = Form("%s%s_eta%02d",
                    ptSl.shortName.Data(), aSl.shortName.Data(), ie);
                TString dname = L2Name::ObjectName(cone, "A_data",
                    {L2Name::EtaModeKey(false), etaKey, ptKey, alphaKey});
                TString mname = L2Name::ObjectName(cone, "A_mc",
                    {L2Name::EtaModeKey(false), etaKey, ptKey, alphaKey});

                TH1D* hd = GetHAny(dData, {dname, cone + "_A_data_" + oldSfx});
                TH1D* hm = GetHAny(dMC,   {mname, cone + "_A_mc_"   + oldSfx});

                if (!hd || !hm || hd->GetEntries() < kMinEntriesPlot) {
                    if (hd) delete hd;
                    if (hm) delete hm;
                    pb.Update();
                    pb.Update();
                    pb.Update();
                    continue;
                }

                Long64_t nData = (Long64_t)hd->GetEntries();
                Long64_t nMC   = (Long64_t)hm->GetEntries();

                TH1D* hdc = (TH1D*)hd->Clone(dname + "_c"); hdc->SetDirectory(0);
                TH1D* hmc = (TH1D*)hm->Clone(mname + "_c"); hmc->SetDirectory(0);

                NormalizeDensity(hdc);
                NormalizeDensity(hmc);

                const double ymax = std::max(hdc->GetMaximum(), hmc->GetMaximum()) * 5.0;
                const double ymin = std::max(1e-6, MinPositiveBin({hdc, hmc}) * 0.5);

                double etalo = kAbsEtaEdges[ie];
                double etahi = kAbsEtaEdges[ie + 1];

                auto drawInfo = [&]() {
                    TLatex* tex = new TLatex();
                    tex->SetNDC();
                    tex->SetTextSize(0.031);
                    tex->SetTextFont(42);
                    tex->DrawLatex(0.57, 0.865, cone);
                    tex->DrawLatex(0.57, 0.810, Form("%.3f < |#eta^{probe}| < %.3f", etalo, etahi));
                    tex->DrawLatex(0.57, 0.755, ptSl.title.Data());
                    tex->DrawLatex(0.57, 0.700, aSl.title.Data());
                };

                auto makeLegend = [&]() {
                    TLegend* leg = new TLegend(0.16, 0.60, 0.48, 0.84);
                    leg->SetBorderSize(0);
                    leg->SetFillColorAlpha(kWhite, 0.75);
                    leg->SetFillStyle(1001);
                    leg->SetTextSize(0.031);
                    leg->AddEntry(hdc, Form("Data  (N = %lld)", nData), "lp");
                    leg->AddEntry(hmc, Form("MC    (N = %lld)", nMC), "lp");
                    return leg;
                };

                auto drawTruncPlot = [&](double fraction, const TString& tag, const TString& label) {
                    TString cvName = Form("adist_%s_%s_%s_%s_%s",
                        cone.Data(), etaKey.Data(), ptKey.Data(), alphaKey.Data(), tag.Data());
                    TCanvas* c = new TCanvas(cvName, "", 800, 600);
                    c->SetRealAspectRatio(kAspectRatio);
                    c->SetLogy();
                    c->SetLeftMargin(0.13);
                    c->SetRightMargin(0.05);

                    DrawAsymBase(hdc, hmc, "A", ymin, ymax);
                    DrawAsymHeader();
                    drawInfo();

                    TLine *ldLo = nullptr, *ldHi = nullptr, *lmLo = nullptr, *lmHi = nullptr;
                    DrawTruncLines(hdc, fraction, kBlack, ymin, ymax, ldLo, ldHi);
                    DrawTruncLines(hmc, fraction, kRed + 1, ymin, ymax, lmLo, lmHi);

                    TLegend* leg = makeLegend();
                    leg->AddEntry(ldLo, Form("Data %s", label.Data()), "l");
                    leg->AddEntry(lmLo, Form("MC %s", label.Data()), "l");
                    leg->Draw();

                    SavePlot(c, outDir, cone, "adist", {etaKey, ptKey, alphaKey}, cvName);
                    pb.Update();
                    delete c;
                };

                drawTruncPlot(0.90, "trunc90", "trunc. 90%");
                drawTruncPlot(0.95, "trunc95", "trunc. 95%");

                {
                    TString cvName = Form("adist_%s_%s_%s_%s_gauss",
                        cone.Data(), etaKey.Data(), ptKey.Data(), alphaKey.Data());
                    TCanvas* c = new TCanvas(cvName, "", 800, 600);
                    c->SetRealAspectRatio(kAspectRatio);
                    c->SetLogy();
                    c->SetLeftMargin(0.13);
                    c->SetRightMargin(0.05);

                    DrawAsymBase(hdc, hmc, "A", ymin, ymax);
                    DrawAsymHeader();
                    drawInfo();

                    TF1* fd = FitGaussianGuide(hdc, cvName + "_data_fit", kBlack);
                    TF1* fm = FitGaussianGuide(hmc, cvName + "_mc_fit", kRed + 1);
                    if (fd) fd->Draw("same");
                    if (fm) fm->Draw("same");

                    TLegend* leg = makeLegend();
                    if (fd) leg->AddEntry(fd, "Data Gaussian fit", "l");
                    if (fm) leg->AddEntry(fm, "MC Gaussian fit", "l");
                    leg->Draw();

                    SavePlot(c, outDir, cone, "adist", {etaKey, ptKey, alphaKey}, cvName);
                    pb.Update();
                    delete fd;
                    delete fm;
                    delete c;
                }

                delete hdc;
                delete hmc;
                delete hd;
                delete hm;
            }
        }
    }
}

// ============================================================
// Plot type 4: R_data and R_MC overlay with ratio panel
//
// For each (cone, method, alpha slice, pT slice): one two-panel canvas with
//   top panel  — R_data (blue) and R_MC (red) vs |eta|, reference at 1
//   bottom panel — R_data/R_MC vs |eta|, reference at 1
// ============================================================

static void PlotROverlay(TFile* fIn, const TString& outDir,
                         const TString& cone, const BinningConfig& bins,
                         ProgressBar& pb) {
    TDirectory* dRvals = (TDirectory*)fIn->Get(cone + "/Rvals");
    if (!dRvals) dRvals = (TDirectory*)fIn->Get(cone + "_Rvals");

    const int nPt    = (int)bins.ptavgSlices.size();
    const int nAlpha = (int)bins.alphaSlices.size();

    for (int m = 0; m < kNMethods; m++) {
        for (int ia = 0; ia < nAlpha; ia++) {
            const auto& aSl = bins.alphaSlices[ia];
            for (int ip = 0; ip < nPt; ip++) {
                const auto& ptSl = bins.ptavgSlices[ip];

                TString ptKey = L2Name::PtKey(ptSl);
                TString alphaKey = L2Name::AlphaKey(aSl);
                TString rdName = L2Name::ObjectName(cone, "R_data",
                    {L2Name::EtaModeKey(false), ptKey, alphaKey}, {kMethodKeys[m]});
                TString rmName = L2Name::ObjectName(cone, "R_mc",
                    {L2Name::EtaModeKey(false), ptKey, alphaKey}, {kMethodKeys[m]});
                TString oldRdName = Form("%s_R_data_%s%s%s",
                    cone.Data(), kMethodKeys[m], ptSl.shortName.Data(), aSl.shortName.Data());
                TString oldRmName = Form("%s_R_mc_%s%s%s",
                    cone.Data(), kMethodKeys[m], ptSl.shortName.Data(), aSl.shortName.Data());

                TH1D* hRd = GetHAny(dRvals, {rdName, oldRdName});
                TH1D* hRm = GetHAny(dRvals, {rmName, oldRmName});

                if (!hRd || !hRm) {
                    if (hRd) delete hRd;
                    if (hRm) delete hRm;
                    pb.Update();
                    continue;
                }

                TH1D* hRdc = (TH1D*)hRd->Clone(rdName + "_c"); hRdc->SetDirectory(0);
                TH1D* hRmc = (TH1D*)hRm->Clone(rmName + "_c"); hRmc->SetDirectory(0);
                TH1D* hRat = RatioH(hRdc, hRmc, rdName + "_rat");

                const TString cvName = Form("roverlay_%s_%s_%s_%s",
                    cone.Data(), kMethodKeys[m], ptKey.Data(), alphaKey.Data());
                TwoPad cv = MakeTwoPad(cvName);

                // ---- main pad ----
                cv.main->cd();
                cv.main->SetGridx(); cv.main->SetGridy();

                StyleH(hRdc, kBlue+1,  20, 1.5f);
                StyleH(hRmc, kRed+1,   21, 1.5f);

                auto [ylo, yhi] = YRange({hRdc, hRmc});
                hRmc->GetYaxis()->SetRangeUser(ylo, yhi);
                hRmc->GetYaxis()->SetTitle("R");
                hRmc->GetYaxis()->SetTitleSize(0.065);
                hRmc->GetYaxis()->SetTitleOffset(1.0);
                hRmc->GetYaxis()->SetLabelSize(0.055);
                hRmc->GetXaxis()->SetLabelSize(0.0);
                hRmc->GetXaxis()->SetTitle("");
                hRmc->SetTitle("");

                hRmc->Draw("E1");
                hRdc->Draw("E1 same");
                RefLine(cv.main, (double)kAbsEtaEdges.front(), (double)kAbsEtaEdges.back(), 1.0);

                TLegend* leg = new TLegend(0.16, 0.14, 0.50, 0.28);
                leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.048);
                leg->AddEntry(hRmc, "R_{MC}",   "lp");
                leg->AddEntry(hRdc, "R_{data}", "lp");
                leg->Draw();

                TLatex* tex = new TLatex();
                tex->SetNDC(); tex->SetTextSize(0.050); tex->SetTextFont(62);
                tex->DrawLatex(0.16, 0.91, Form("%s  |  %s  |  %s  |  %s",
                    cone.Data(), kMethodLabels[m], ptSl.title.Data(), aSl.title.Data()));

                // ---- ratio pad ----
                cv.ratio->cd();
                cv.ratio->SetGridx(); cv.ratio->SetGridy();

                StyleH(hRat, kBlack, 20, 1.5f);
                TuneRatio(hRat, "|#eta|", "Data/MC", 0.93, 1.07);
                hRat->Draw("E1");
                RefLine(cv.ratio, (double)kAbsEtaEdges.front(), (double)kAbsEtaEdges.back(), 1.0);

                cv.c->cd();
                SavePlot(cv.c, outDir, cone, "roverlay", {ptKey, alphaKey}, cvName);
                pb.Update();

                delete hRd; delete hRm;
                delete cv.c;   // cascade-deletes hRdc, hRmc, hRat, leg, tex
            }
        }
    }
}

// ============================================================
// Plot type 5: Alpha fit plots
//
// For each (cone, method, pT slice, eta bin): one canvas showing
//   all 9 alpha threshold points with fit line drawn only through [0, 0.31],
//   points at alpha > 0.30 are shown but outside the fit line.
// ============================================================

static void PlotAlphaFit(TFile* fIn, const TString& outDir,
                         const TString& cone, const BinningConfig& bins,
                         ProgressBar& pb) {
    TDirectory* dGraphs = (TDirectory*)fIn->Get(cone + "/graphs");
    if (!dGraphs) dGraphs = (TDirectory*)fIn->Get(cone + "_graphs");

    const int nPt  = (int)bins.ptavgSlices.size();
    const int nEta = (int)kAbsEtaEdges.size() - 1;

    for (int m = 0; m < kNMethods; m++) {
        for (int ip = 0; ip < nPt; ip++) {
            const auto& ptSl = bins.ptavgSlices[ip];
            for (int ie = 0; ie < nEta; ie++) {
                TString etaKey = L2Name::EtaKey(ie);
                TString ptKey = L2Name::PtKey(ptSl);
                TString gname = L2Name::ObjectName(cone, "R",
                    {L2Name::EtaModeKey(false), etaKey, ptKey}, {kMethodKeys[m]});
                TString oldGName = Form("%s_R_%s%s_eta%02d",
                    cone.Data(), kMethodKeys[m], ptSl.shortName.Data(), ie);

                TGraphErrors* gr = GetGraphAny(dGraphs, {gname, oldGName});

                if (!gr || gr->GetN() < 2) {
                    pb.Update();
                    continue;
                }
                TGraphErrors* gc = (TGraphErrors*)gr->Clone(gname + "_c");

                const TString cvName = Form("alphafit_%s_%s_%s_%s",
                    cone.Data(), kMethodKeys[m], etaKey.Data(), ptKey.Data());
                TCanvas* c = new TCanvas(cvName, "", 800, 600);
            RealAspectRatio(c);
                c->SetLeftMargin(0.13);
                c->SetGridx(); c->SetGridy();

                gc->SetMarkerStyle(20);
                gc->SetMarkerColor(ptSl.color);
                gc->SetLineColor(ptSl.color);
                gc->SetMarkerSize(0.9);

                gc->GetXaxis()->SetTitle("#alpha threshold");
                gc->GetYaxis()->SetTitle("R_{data}/R_{MC}");
                gc->GetXaxis()->CenterTitle();
                gc->GetYaxis()->CenterTitle();
                gc->GetXaxis()->SetLimits(0.0, 0.50);
                gc->SetTitle("");

                gc->Draw("AP");   // embedded fit function draws automatically

                // vertical reference at x=0.30 to mark the fit boundary
                double ylo = gc->GetHistogram()->GetMinimum();
                double yhi = gc->GetHistogram()->GetMaximum();
                TLine* vl = new TLine(0.30, ylo, 0.30, yhi);
                vl->SetLineStyle(2); vl->SetLineColor(kGray+2); vl->SetLineWidth(1);
                vl->Draw();

                // horizontal reference at y=1
                TLine* hl = new TLine(0.0, 1.0, 0.50, 1.0);
                hl->SetLineStyle(2); hl->SetLineColor(kGray+2); hl->SetLineWidth(1);
                hl->Draw();

                TLatex* tex = new TLatex();
                tex->SetNDC(); tex->SetTextSize(0.042); tex->SetTextFont(62);
                tex->DrawLatex(0.14, 0.92, Form("%s  |  %s  |  %s  |  |#eta| bin %d",
                    cone.Data(), kMethodLabels[m], ptSl.title.Data(), ie));

                SavePlot(c, outDir, cone, "alpha", {etaKey, ptKey}, cvName);
                pb.Update();

                delete gc;
                delete c;
            }
        }
    }
}
//
// For each (cone, method, ptavg slice): one canvas with
//   top panel  — full-eta corrections + |eta| reflected, reference at 1
//   bottom panel — (full eta) / (reflected |eta|)
// ============================================================

static void PlotEtaSym(TFile* fIn, const TString& outDir,
                       const TString& cone, const BinningConfig& bins,
                       ProgressBar& pb) {
    for (int m = 0; m < kNMethods; m++) {
        for (int ip = 0; ip < (int)bins.ptavgSlices.size(); ip++) {
            const auto& sl = bins.ptavgSlices[ip];
            TString ptKey = L2Name::PtKey(sl);

            const TString nameAbs = L2Name::ObjectName(cone, "intercept",
                {L2Name::EtaModeKey(false), ptKey}, {kMethodKeys[m]});
            const TString nameFull = L2Name::ObjectName(cone, "intercept",
                {L2Name::EtaModeKey(true), ptKey}, {kMethodKeys[m]});
            const TString oldNameAbs = Form("%s_intercept_%s%s",
                cone.Data(), kMethodKeys[m], sl.shortName.Data());
            const TString oldNameFull = oldNameAbs + "_fulleta";

            TH1D* hAbs  = GetHAny(fIn, {cone + "/" + nameAbs,  nameAbs,  oldNameAbs});
            TH1D* hFull = GetHAny(fIn, {cone + "/" + nameFull, nameFull, oldNameFull});

            if (!hAbs || !hFull) {
                delete hAbs; delete hFull;
                pb.Update();
                continue;
            }

            // Reflect |eta| onto full-eta axis then compute ratio
            TH1D* hRefl  = Reflect(hAbs, nameAbs + "_refl");
            TH1D* hRatio = RatioH(hFull, hRefl, nameAbs + "_rat");

            const TString cvName = Form("etasym_%s_%s_%s",
                cone.Data(), kMethodKeys[m], ptKey.Data());
            TwoPad cv = MakeTwoPad(cvName);

            // ---- main pad ----
            cv.main->cd();
            cv.main->SetGridx();
            cv.main->SetGridy();

            StyleH(hFull, kColFull, 20, 2.f);
            StyleH(hRefl, kColRefl, 21, 2.f);
            hRefl->SetLineStyle(2);

            auto [ylo, yhi] = YRange({hFull, hRefl});
            hFull->GetYaxis()->SetRangeUser(ylo, yhi);
            hFull->GetYaxis()->SetTitle("R_{data}/R_{MC} at #alpha#rightarrow0");
            hFull->GetYaxis()->SetTitleSize(0.055);
            hFull->GetYaxis()->SetTitleOffset(1.10);
            hFull->GetYaxis()->SetLabelSize(0.050);
            hFull->GetXaxis()->SetLabelSize(0.0);
            hFull->GetXaxis()->SetTitle("");
            hFull->SetTitle("");

            hFull->Draw("E1");
            hRefl->Draw("E1 same");
            RefLine(cv.main, kEtaEdges.front(), kEtaEdges.back(), 1.0);

            TLegend* leg = new TLegend(0.16, 0.14, 0.56, 0.28);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.048);
            leg->AddEntry(hFull, "Full #eta",        "lp");
            leg->AddEntry(hRefl, "|#eta| reflected", "lp");
            leg->Draw();

            TLatex* tex = new TLatex();
            tex->SetNDC();
            tex->SetTextSize(0.051);
            tex->SetTextFont(62);
            tex->DrawLatex(0.16, 0.91, Form("%s   |   %s   |   %s",
                cone.Data(), kMethodLabels[m], sl.title.Data()));

            // ---- ratio pad ----
            cv.ratio->cd();
            cv.ratio->SetGridx();
            cv.ratio->SetGridy();

            StyleH(hRatio, kBlack, 20, 1.5f);
            TuneRatio(hRatio, "#eta", "Full / Refl.", 0.975, 1.025);
            hRatio->Draw("E1");
            RefLine(cv.ratio, kEtaEdges.front(), kEtaEdges.back(), 1.0);

            cv.c->cd();
            SavePlot(cv.c, outDir, cone, "etasym", {ptKey}, cvName);
            pb.Update();

            // hAbs was not drawn — delete manually
            delete hAbs;
            // hFull, hRefl, hRatio, leg, tex drawn on pads — canvas cascade deletes them
            delete cv.c;
        }
    }
}

// ============================================================
// Plot type 2: method comparison (gauss vs trunc90 vs trunc95)
//
// For each (cone, ptavg slice, eta type): one canvas with
//   top panel  — all three methods overlaid
//   bottom panel — trunc90/gauss and trunc95/gauss
// ============================================================

static void PlotMethodComp(TFile* fIn, const TString& outDir,
                           const TString& cone, const BinningConfig& bins,
                           bool fullEta, ProgressBar& pb) {
    const TString suffix   = fullEta ? "_fulleta" : "";
    const TString etaMode  = L2Name::EtaModeKey(fullEta);
    const TString etaLabel = fullEta ? "Full #eta" : "|#eta|";
    const double  xMin     = fullEta ? kEtaEdges.front()    : (double)kAbsEtaEdges.front();
    const double  xMax     = fullEta ? kEtaEdges.back()     : (double)kAbsEtaEdges.back();
    const TString xTitle   = fullEta ? "#eta" : "|#eta|";

    for (int ip = 0; ip < (int)bins.ptavgSlices.size(); ip++) {
        const auto& sl = bins.ptavgSlices[ip];

        std::vector<TH1D*> hists(kNMethods, nullptr);
        for (int m = 0; m < kNMethods; m++) {
            TString name = L2Name::ObjectName(cone, "intercept",
                {etaMode, L2Name::PtKey(sl)}, {kMethodKeys[m]});
            TString oldName = Form("%s_intercept_%s%s%s",
                cone.Data(), kMethodKeys[m], sl.shortName.Data(), suffix.Data());
            hists[m] = GetHAny(fIn, {cone + "/" + name, name, oldName});
        }

        if (!hists[0]) {
            for (auto* h : hists) delete h;
            pb.Update();
            continue;
        }

        // ratios to gauss for non-null trunc methods
        std::vector<TH1D*> ratios;
        std::vector<int>   ratioIdx;
        for (int m = 1; m < kNMethods; m++) {
            if (!hists[m]) continue;
            TString rname = Form("%s_mcomp%s%s_r%d",
                cone.Data(), sl.shortName.Data(), suffix.Data(), m);
            ratios.push_back(RatioH(hists[m], hists[0], rname));
            ratioIdx.push_back(m);
        }

        TString ptKey = L2Name::PtKey(sl);
        const TString cvName = Form("methods_%s_%s_%s",
            cone.Data(), etaMode.Data(), ptKey.Data());
        TwoPad cv = MakeTwoPad(cvName);

        // ---- main pad ----
        cv.main->cd();
        cv.main->SetGridx();
        cv.main->SetGridy();

        auto [ylo, yhi] = YRange({hists[0], hists[1], hists[2]});
        bool first = true;

        TLegend* leg = new TLegend(0.16, 0.14, 0.58, 0.14 + 0.065 * kNMethods);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.046);

        for (int m = 0; m < kNMethods; m++) {
            if (!hists[m]) continue;
            StyleH(hists[m], kMethodColors[m], kMethodStyles[m], 2.f);
            hists[m]->GetYaxis()->SetRangeUser(ylo, yhi);
            hists[m]->GetYaxis()->SetTitle("R_{data}/R_{MC} at #alpha#rightarrow0");
            hists[m]->GetYaxis()->SetTitleSize(0.055);
            hists[m]->GetYaxis()->SetTitleOffset(1.10);
            hists[m]->GetYaxis()->SetLabelSize(0.050);
            hists[m]->GetXaxis()->SetLabelSize(0.0);
            hists[m]->GetXaxis()->SetTitle("");
            hists[m]->SetTitle("");
            hists[m]->Draw(first ? "E1" : "E1 same");
            first = false;
            leg->AddEntry(hists[m], kMethodLabels[m], "lp");
        }

        RefLine(cv.main, xMin, xMax, 1.0);
        leg->Draw();

        TLatex* tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.051);
        tex->SetTextFont(62);
        tex->DrawLatex(0.16, 0.91, Form("%s   |   %s   |   %s",
            cone.Data(), etaLabel.Data(), sl.title.Data()));

        // ---- ratio pad ----
        cv.ratio->cd();
        cv.ratio->SetGridx();
        cv.ratio->SetGridy();

        TLegend* rleg = new TLegend(0.16, 0.62, 0.58, 0.95);
        rleg->SetBorderSize(0);
        rleg->SetFillStyle(0);
        rleg->SetTextSize(0.115);

        bool firstR = true;
        for (int k = 0; k < (int)ratios.size(); k++) {
            const int m = ratioIdx[k];
            StyleH(ratios[k], kMethodColors[m], kMethodStyles[m], 1.5f);
            TuneRatio(ratios[k], xTitle, "/ Gauss", 0.993, 1.007);
            ratios[k]->Draw(firstR ? "E1" : "E1 same");
            firstR = false;
            rleg->AddEntry(ratios[k], kMethodLabels[m], "lp");
        }
        if (!firstR) {
            RefLine(cv.ratio, xMin, xMax, 1.0);
            rleg->Draw();
        }

        cv.c->cd();
        SavePlot(cv.c, outDir, cone, "methods", {etaMode, ptKey}, cvName);
        pb.Update();

        // ratios and hists were drawn on pads — canvas cascade deletes them
        // rleg, leg, tex also drawn — also cascade-deleted
        delete cv.c;
    }
}

// ============================================================
// Plot type 6: Final extrapolated values — all pT slices overlaid
//
// For each (cone, method): two canvases
//   finals_{cone}_{method}_abseta  — R_data/R_MC at alpha→0 vs |eta|, all pT bins overlaid
//   finals_{cone}_{method}_fulleta — same vs full eta
// PlotEtaSym does the |eta|-vs-fulleta symmetry check per pT slice; this overlays pT slices.
// ============================================================

static void PlotFinals(TFile* fIn, const TString& outDir,
                       const TString& cone, const BinningConfig& bins,
                       ProgressBar& pb) {
    for (int m = 0; m < kNMethods; m++) {
        for (int ieta = 0; ieta < 2; ieta++) {   // 0 = |eta|, 1 = full eta
            const bool   fullEta  = (ieta == 1);
            const TString suffix  = fullEta ? "_fulleta" : "";
            const TString xTitle  = fullEta ? "#eta" : "|#eta|";
            const double  xMin    = fullEta ? kEtaEdges.front()    : (double)kAbsEtaEdges.front();
            const double  xMax    = fullEta ? kEtaEdges.back()     : (double)kAbsEtaEdges.back();

            std::vector<TH1D*> hists;
            for (const auto& ptSl : bins.ptavgSlices) {
                TString name = L2Name::ObjectName(cone, "intercept",
                    {L2Name::EtaModeKey(fullEta), L2Name::PtKey(ptSl)}, {kMethodKeys[m]});
                TString oldName = Form("%s_intercept_%s%s%s",
                    cone.Data(), kMethodKeys[m], ptSl.shortName.Data(), suffix.Data());
                hists.push_back(GetHAny(fIn, {cone + "/" + name, name, oldName}));
            }

            bool anyValid = false;
            for (auto* h : hists) if (h) { anyValid = true; break; }
            if (!anyValid) {
                for (auto* h : hists) delete h;
                pb.Update();
                continue;
            }

            const TString etaMode = L2Name::EtaModeKey(fullEta);
            const TString cvName = Form("finals_%s_%s_%s",
                cone.Data(), kMethodKeys[m], etaMode.Data());
            TCanvas* c = new TCanvas(cvName, "", 800, 600);
            RealAspectRatio(c);
            c->SetLeftMargin(0.13);
            c->SetGridx();
            c->SetGridy();

            auto [ylo, yhi] = YRange(hists);

            TLegend* leg = new TLegend(0.60, 0.68, 0.93, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.038);

            bool first = true;
            for (int ip = 0; ip < (int)bins.ptavgSlices.size(); ip++) {
                if (!hists[ip]) continue;
                const auto& ptSl = bins.ptavgSlices[ip];
                StyleH(hists[ip], ptSl.color, kMethodStyles[m], 1.5f);
                hists[ip]->GetYaxis()->SetRangeUser(ylo, yhi);
                hists[ip]->GetYaxis()->SetTitle("R_{data}/R_{MC} at #alpha#rightarrow0");
                hists[ip]->GetYaxis()->SetTitleSize(0.052);
                hists[ip]->GetYaxis()->SetTitleOffset(1.15);
                hists[ip]->GetYaxis()->SetLabelSize(0.048);
                hists[ip]->GetXaxis()->SetTitle(xTitle);
                hists[ip]->GetXaxis()->SetTitleSize(0.052);
                hists[ip]->GetXaxis()->SetLabelSize(0.048);
                hists[ip]->GetXaxis()->CenterTitle();
                hists[ip]->GetYaxis()->CenterTitle();
                hists[ip]->SetTitle("");
                hists[ip]->Draw(first ? "E1" : "E1 same");
                first = false;
                leg->AddEntry(hists[ip], ptSl.title, "lp");
            }

            TLine* rl = new TLine(xMin, 1.0, xMax, 1.0);
            rl->SetLineStyle(2);
            rl->SetLineColor(kGray + 2);
            rl->SetLineWidth(1);
            rl->Draw();

            leg->Draw();

            TLatex* tex = new TLatex();
            tex->SetNDC();
            tex->SetTextSize(0.048);
            tex->SetTextFont(62);
            tex->DrawLatex(0.14, 0.92, Form("%s  |  %s  |  %s",
                cone.Data(), kMethodLabels[m], xTitle.Data()));

            SavePlot(c, outDir, cone, "finals", {etaMode}, cvName);
            pb.Update();

            delete c;   // cascade-deletes hists, leg, tex, rl
        }
    }
}

// ============================================================
// Plot type 7: Step-1 jet kinematics
//
// Reads the TH3D(eta, phi, pT) control histograms written by runAsymmetry:
//   {cone}_incl, {cone}_tag, {cone}_probe
//
// For each cone and collection, writes pT, eta, phi projections plus eta-phi
// maps above a few pT thresholds. This mode expects a Step-1 runAsymmetry file,
// not a Step-2 residuals file.
// ============================================================

static constexpr int kNKinematicsCollections = 3;
static const char* const kKinematicsCollections[] = { "incl", "tag", "probe" };

static constexpr int kNKinematicsPtMins = 3;
static const double kKinematicsPtMins[] = { 40.0, 100.0, 200.0 };

static TH1D* ProjectTH3D1D(TH3D* h3, const char* axis, const TString& name) {
    TDirectory::TContext nodir(nullptr);
    TH1D* h = (TH1D*)h3->Project3D(axis);
    h->SetName(name);
    h->SetDirectory(0);
    return h;
}

static TH2D* ProjectEtaPhi(TH3D* h3, double ptMin, const TString& name) {
    h3->GetZaxis()->SetRangeUser(ptMin, h3->GetZaxis()->GetXmax());
    TDirectory::TContext nodir(nullptr);
    TH2D* h = (TH2D*)h3->Project3D("yx");
    h->SetName(name);
    h->SetDirectory(0);
    h3->GetZaxis()->SetRange(0, 0);
    return h;
}

static void DrawKinematics1D(TH1D* h, const TString& xTitle, bool logy) {
    h->SetTitle("");
    h->GetXaxis()->SetTitle(xTitle);
    h->GetYaxis()->SetTitle("Jets");
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->SetLineColor(HiroshigeNightBlue());
    h->SetMarkerColor(HiroshigeNightBlue());
    h->SetMarkerStyle(20);
    h->SetLineWidth(2);
    h->Draw("hist");
    if (logy && h->GetMaximum() > 0) {
        h->SetMinimum(std::max(0.5, h->GetMinimum(0.0) * 0.5));
        gPad->SetLogy();
    }
}

static TString PtMinKey(double ptMin) {
    return Form("ptmin_%g", ptMin);
}

static void PlotKinematics(TFile* fIn, const TString& outDir,
                           const TString& cone, ProgressBar& pb) {
    const int plotsPerCollection = 3 + kNKinematicsPtMins;
    {
        TH3D* hTest = (TH3D*)fIn->Get(cone + "/" + cone + "_incl");
        if (!hTest) hTest = (TH3D*)fIn->Get(cone + "_incl");
        if (!hTest) {
            for (int i = 0; i < kNKinematicsCollections * plotsPerCollection; i++) pb.Update();
            return;
        }
    }
    for (int ic = 0; ic < kNKinematicsCollections; ic++) {
        const TString coll = kKinematicsCollections[ic];
        TH3D* h3 = (TH3D*)fIn->Get(cone + "/" + cone + "_" + coll);
        if (!h3) h3 = (TH3D*)fIn->Get(cone + "_" + coll);

        if (!h3) {
            for (int i = 0; i < plotsPerCollection; i++) pb.Update();
            continue;
        }

        h3->GetXaxis()->SetRange(0, 0);
        h3->GetYaxis()->SetRange(0, 0);
        h3->GetZaxis()->SetRange(0, 0);

        {
            TString cvName = Form("kinematics_%s_%s_pt", cone.Data(), coll.Data());
            TH1D* h = ProjectTH3D1D(h3, "z", cvName + "_h");
            TCanvas* c = new TCanvas(cvName, "", 800, 600);
            RealAspectRatio(c);
            c->SetLeftMargin(0.13);
            DrawKinematics1D(h, "p_{T} [GeV/c]", true);
            SavePlot(c, outDir, cone, "kinematics", {coll}, cvName);
            delete c;
            delete h;
            pb.Update();
        }

        {
            TString cvName = Form("kinematics_%s_%s_eta", cone.Data(), coll.Data());
            TH1D* h = ProjectTH3D1D(h3, "x", cvName + "_h");
            TCanvas* c = new TCanvas(cvName, "", 800, 600);
            RealAspectRatio(c);
            c->SetLeftMargin(0.13);
            DrawKinematics1D(h, "#eta", false);
            SavePlot(c, outDir, cone, "kinematics", {coll}, cvName);
            delete c;
            delete h;
            pb.Update();
        }

        {
            TString cvName = Form("kinematics_%s_%s_phi", cone.Data(), coll.Data());
            TH1D* h = ProjectTH3D1D(h3, "y", cvName + "_h");
            TCanvas* c = new TCanvas(cvName, "", 800, 600);
            RealAspectRatio(c);
            c->SetLeftMargin(0.13);
            DrawKinematics1D(h, "#phi", false);
            SavePlot(c, outDir, cone, "kinematics", {coll}, cvName);
            delete c;
            delete h;
            pb.Update();
        }

        for (int ip = 0; ip < kNKinematicsPtMins; ip++) {
            const double ptMin = kKinematicsPtMins[ip];
            const TString ptKey = PtMinKey(ptMin);
            TString cvName = Form("kinematics_%s_%s_eta_phi_%s",
                cone.Data(), coll.Data(), ptKey.Data());

            TH2D* h = ProjectEtaPhi(h3, ptMin, cvName + "_h");
            TCanvas* c = new TCanvas(cvName, "", 800, 600);
            RealAspectRatio(c);
            c->SetLeftMargin(0.12);
            c->SetRightMargin(0.16);
            h->SetTitle("");
            h->GetXaxis()->SetTitle("#eta");
            h->GetYaxis()->SetTitle("#phi");
            h->GetZaxis()->SetTitle("Jets");
            h->GetXaxis()->CenterTitle();
            h->GetYaxis()->CenterTitle();
            h->Draw("colz");

            TLatex* tex = new TLatex();
            tex->SetNDC();
            tex->SetTextSize(0.040);
            tex->SetTextFont(62);
            tex->DrawLatex(0.13, 0.93, Form("%s  |  %s jets  |  p_{T} > %.0f GeV/c",
                cone.Data(), coll.Data(), ptMin));

            SavePlot(c, outDir, cone, "kinematics", {coll, ptKey}, cvName);
            delete c;
            delete h;
            pb.Update();
        }
    }
}

// ============================================================
// Plot type: Event-level QA
//
// Reads event-level control histograms written by runAsymmetry:
//   hvz_all   — vz before any cuts
//   hvz       — vz after vz + vertex filter cuts
//   hfilt     — pprimaryVertexFilter (DATA only)
//   h_hlt_j80 — HLT_AK4PFJet80 bits (hard-probe DATA only)
//
// Expects a Step-1 runAsymmetry output file; gracefully skips
// histograms that are absent (e.g. hfilt/h_hlt_j80 on MC).
// ============================================================

static void PlotEvent(TFile* fIn, const TString& outDir, ProgressBar& pb) {

    // ---- vz: all events vs after cuts ----
    {
        TH1D* hall = (TH1D*)fIn->Get("hvz_all");
        TH1D* hvz  = (TH1D*)fIn->Get("hvz");
        if (hall && hvz) {
            TH1D* hallc = (TH1D*)hall->Clone("hvz_all_c"); hallc->SetDirectory(0);
            TH1D* hvzc  = (TH1D*)hvz ->Clone("hvz_c");     hvzc ->SetDirectory(0);

            hallc->SetLineColor(kGray + 2);
            hallc->SetFillColor(kGray + 1);
            hallc->SetFillStyle(1001);
            hallc->SetLineWidth(1);
            StyleH(hvzc, HiroshigeNightBlue(), 20, 2.0f);

            double ymax = std::max(hallc->GetMaximum(), hvzc->GetMaximum()) * 1.2;
            hallc->SetMaximum(ymax);
            hallc->SetTitle("");
            hallc->GetXaxis()->SetTitle("v_{z} [cm]");
            hallc->GetYaxis()->SetTitle("Events");
            hallc->GetXaxis()->CenterTitle();
            hallc->GetYaxis()->CenterTitle();

            TCanvas* c = new TCanvas("event_vz", "", 800, 600);
            RealAspectRatio(c);
            c->SetLeftMargin(0.13); c->SetGridx(); c->SetGridy();
            hallc->Draw("hist");
            hvzc->Draw("E1 same");

            TLegend* leg = new TLegend(0.62, 0.76, 0.93, 0.90);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
            leg->AddEntry(hallc, Form("All events (N = %lld)", (Long64_t)hall->GetEntries()), "f");
            leg->AddEntry(hvzc,  Form("After cuts  (N = %lld)", (Long64_t)hvz->GetEntries()),  "lp");
            leg->Draw();

            SavePlot(c, outDir, "", "event", {}, "event_vz");
            delete c;
        }
        pb.Update();
    }

    // ---- primary vertex filter ----
    {
        TH1I* hfilt = (TH1I*)fIn->Get("hfilt");
        if (hfilt) {
            TCanvas* c = new TCanvas("event_ppvF", "", 800, 600);
            RealAspectRatio(c);
            c->SetLeftMargin(0.15);
            hfilt->SetTitle("");
            hfilt->GetXaxis()->SetTitle("pprimaryVertexFilter");
            hfilt->GetYaxis()->SetTitle("Events");
            hfilt->GetXaxis()->CenterTitle();
            hfilt->GetYaxis()->CenterTitle();
            hfilt->SetLineColor(HiroshigeNightBlue());
            hfilt->SetFillColor(HiroshigeLightBlue());
            hfilt->SetFillStyle(1001);
            hfilt->Draw("hist");
            SavePlot(c, outDir, "", "event", {}, "event_ppvF");
            delete c;
        }
        pb.Update();
    }

    // ---- HLT trigger ----
    {
        TH1I* htrig = (TH1I*)fIn->Get("h_hlt_j80");
        if (htrig) {
            TCanvas* c = new TCanvas("event_hlt", "", 800, 600);
            RealAspectRatio(c);
            c->SetLeftMargin(0.15);
            htrig->SetTitle("");
            htrig->GetXaxis()->SetTitle("HLT_AK4PFJet80");
            htrig->GetYaxis()->SetTitle("Events");
            htrig->GetXaxis()->CenterTitle();
            htrig->GetYaxis()->CenterTitle();
            htrig->SetLineColor(HiroshigeNightBlue());
            htrig->SetFillColor(HiroshigeLightBlue());
            htrig->SetFillStyle(1001);
            htrig->Draw("hist");
            SavePlot(c, outDir, "", "event", {}, "event_hlt_j80");
            delete c;
        }
        pb.Update();
    }
}

// ============================================================
// Entry point
//
// flags (space-separated keywords, default "all"):
//   "etasym"     — full-eta vs |eta| reflected symmetry check (PlotEtaSym)
//   "methods"    — method comparison: gauss vs trunc90 vs trunc95 (PlotMethodComp)
//   "finals"     — final R_data/R_MC at alpha→0, all pT slices overlaid (PlotFinals)
//   "adist"      — asymmetry distributions per bin with log-y and truncation lines
//   "roverlay"   — R_data and R_MC overlay with ratio panel per alpha/pT
//   "alpha"      — alpha fit plots: all 9 points, fit line through 0.05–0.30
//   "kinematics" — Step-1 inclusive/tag/probe jet kinematics from runAsymmetry output
//   "event"      — Step-1 event-level QA: vz, primary vertex filter, HLT trigger
//   "all"        — run all plots including kinematics and event (default)
// ============================================================

void plotResiduals(TString residualsFile, TString outDir = "", TString flags = "all") {
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetGridColor(kGray + 1);
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(1);
    gErrorIgnoreLevel = kWarning;

    TFile* fIn = TFile::Open(residualsFile, "read");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "Cannot open " << residualsFile << "\n";
        return;
    }

    if (outDir.IsNull()) outDir = MakePlotDir("plots_residuals");
    if (gSystem->mkdir(outDir, true) < 0 && gSystem->AccessPathName(outDir)) {
        std::cerr << "Cannot create output directory: " << outDir << "\n";
        return;
    }

    BinningConfig bins;
    const int nCones    = (int)kConeLabels.size();
    const int nPtSlices = (int)bins.ptavgSlices.size();
    const int nAlpha    = (int)bins.alphaSlices.size();
    const int nEta      = (int)kAbsEtaEdges.size() - 1;

    const bool doAll     = flags.IsNull() || flags == "all";
    const bool doEtaSym  = doAll || flags.Contains("etasym");
    const bool doMethods = doAll || flags.Contains("methods");
    const bool doFinals  = doAll || flags.Contains("finals");
    const bool doAdist   = doAll || flags.Contains("adist");
    const bool doRover   = doAll || flags.Contains("roverlay");
    const bool doAlpha   = doAll || flags.Contains("alpha");
    const bool doKine    = doAll || flags.Contains("kinematics");
    const bool doEvent   = doAll || flags.Contains("event");

    int totalPlots = 0;
    if (doEtaSym)  totalPlots += nCones * kNMethods * nPtSlices;
    if (doMethods) totalPlots += nCones * 2 * nPtSlices;
    if (doFinals)  totalPlots += nCones * kNMethods * 2;
    if (doAdist)   totalPlots += 3 * nCones * nAlpha * nPtSlices * nEta;
    if (doRover)   totalPlots += nCones * kNMethods * nAlpha * nPtSlices;
    if (doAlpha)   totalPlots += nCones * kNMethods * nPtSlices * nEta;
    if (doKine)    totalPlots += nCones * kNKinematicsCollections * (3 + kNKinematicsPtMins);
    if (doEvent)   totalPlots += 3;

    ProgressBar pb("Saving plots:", totalPlots);

    if (doEvent) PlotEvent(fIn, outDir, pb);

    for (const TString& cone : kConeLabels) {
        if (doEtaSym)  PlotEtaSym   (fIn, outDir, cone, bins,        pb);
        if (doMethods) PlotMethodComp(fIn, outDir, cone, bins, false, pb);
        if (doMethods) PlotMethodComp(fIn, outDir, cone, bins, true,  pb);
        if (doFinals)  PlotFinals   (fIn, outDir, cone, bins,        pb);
        if (doAdist)   PlotAsymDist (fIn, outDir, cone, bins,        pb);
        if (doRover)   PlotROverlay (fIn, outDir, cone, bins,        pb);
        if (doAlpha)   PlotAlphaFit (fIn, outDir, cone, bins,        pb);
        if (doKine)    PlotKinematics(fIn, outDir, cone,             pb);
    }

    pb.Finish();
    fIn->Close();
}

#ifndef __CLING__
#include <iostream>
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: plotResiduals <residuals.root> [out_dir] [flags]\n"
                  << "  flags: all etasym methods finals adist roverlay alpha kinematics event (space-separated)\n";
        return 1;
    }
    plotResiduals(argv[1], argc >= 3 ? argv[2] : "", argc >= 4 ? argv[3] : "all");
    return 0;
}
#endif
