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
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"

#include "Binning.h"
#include "Colors.h"
#include "Utilities.h"
#include "2024ppRef.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

// ---- method palette ----

static const char* const kMethodKeys[]   = { "gauss",       "trunc90",          "trunc95"          };
static const char* const kMethodLabels[] = { "Gauss fit",   "Trunc. mean 90%",  "Trunc. mean 95%"  };
static const Color_t     kMethodColors[] = { HiroshigeNightBlue, HiroshigeOrange, HiroshigeLightRed };
static const int         kMethodStyles[] = { 20, 21, 22 };   // circle / square / triangle-up
static constexpr int     kNMethods = 3;

// colors for the eta-symmetry comparison
static const Color_t kColFull = HiroshigeLightRed;
static const Color_t kColRefl = HiroshigeNightBlue;

// ============================================================
// Canvas helpers
// ============================================================

struct TwoPad {
    TCanvas* c     = nullptr;
    TPad*    main  = nullptr;
    TPad*    ratio = nullptr;
};

// Main panel: y [0.30, 1.0].  Ratio panel: y [0.00, 0.32].
// Pads are added to the canvas's primitive list; delete canvas to cascade-delete both.
static TwoPad MakeTwoPad(const TString& name) {
    TwoPad cv;
    cv.c = new TCanvas(name, "", 900, 700);
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

// Strip leading underscore from ptSlice.shortName for clean filenames.
static TString SafeKey(const TString& s) {
    return s.BeginsWith("_") ? s(1, s.Length() - 1) : s;
}

// ============================================================
// Plot type 1: |eta| reflected vs full eta
//
// For each (cone, method, ptavg slice): one canvas with
//   top panel  — full-eta corrections + |eta| reflected, reference at 1
//   bottom panel — (full eta) / (reflected |eta|)
// ============================================================

static void PlotEtaSym(TFile* fIn, const TString& outDir,
                       const TString& cone, const BinningConfig& bins) {
    for (int m = 0; m < kNMethods; m++) {
        for (int ip = 0; ip < (int)bins.ptavgSlices.size(); ip++) {
            const auto& sl = bins.ptavgSlices[ip];

            const TString nameAbs  = Form("%s_intercept_%s%s",
                cone.Data(), kMethodKeys[m], sl.shortName.Data());
            const TString nameFull = nameAbs + "_fulleta";

            TH1D* hAbs  = GetH(fIn, nameAbs);
            TH1D* hFull = GetH(fIn, nameFull);

            if (!hAbs || !hFull) {
                std::cout << "  skip eta-sym: " << nameAbs << "\n";
                delete hAbs; delete hFull;
                continue;
            }

            // Reflect |eta| onto full-eta axis then compute ratio
            TH1D* hRefl  = Reflect(hAbs, nameAbs + "_refl");
            TH1D* hRatio = RatioH(hFull, hRefl, nameAbs + "_rat");

            const TString cvName = Form("etasym_%s_%s_%s",
                cone.Data(), kMethodKeys[m], SafeKey(sl.shortName).Data());
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
            cv.c->SaveAs(Form("%s/%s.png", outDir.Data(), cvName.Data()));

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
                           bool fullEta) {
    const TString suffix   = fullEta ? "_fulleta" : "";
    const TString etaLabel = fullEta ? "Full #eta" : "|#eta|";
    const double  xMin     = fullEta ? kEtaEdges.front()    : (double)kAbsEtaEdges.front();
    const double  xMax     = fullEta ? kEtaEdges.back()     : (double)kAbsEtaEdges.back();
    const TString xTitle   = fullEta ? "#eta" : "|#eta|";

    for (int ip = 0; ip < (int)bins.ptavgSlices.size(); ip++) {
        const auto& sl = bins.ptavgSlices[ip];

        std::vector<TH1D*> hists(kNMethods, nullptr);
        for (int m = 0; m < kNMethods; m++) {
            TString name = Form("%s_intercept_%s%s%s",
                cone.Data(), kMethodKeys[m], sl.shortName.Data(), suffix.Data());
            hists[m] = GetH(fIn, name);
        }

        if (!hists[0]) {
            std::cout << "  skip method-comp: gauss missing for " << cone
                      << " " << sl.shortName << suffix << "\n";
            for (auto* h : hists) delete h;
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

        const TString cvName = Form("methods_%s_%s%s",
            cone.Data(), SafeKey(sl.shortName).Data(), suffix.IsNull() ? "_abseta" : suffix.Data());
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
        cv.c->SaveAs(Form("%s/%s.png", outDir.Data(), cvName.Data()));

        // ratios and hists were drawn on pads — canvas cascade deletes them
        // rleg, leg, tex also drawn — also cascade-deleted
        delete cv.c;
    }
}

// ============================================================
// Entry point
// ============================================================

void plotResiduals(TString residualsFile, TString outDir = "") {
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetGridColor(kGray + 1);
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(1);

    TFile* fIn = TFile::Open(residualsFile, "read");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "Cannot open " << residualsFile << "\n";
        return;
    }

    if (outDir.IsNull()) outDir = MakePlotDir("plots_residuals");
    gSystem->mkdir(outDir, true);
    std::cout << "Writing plots to: " << outDir << "\n";

    BinningConfig bins;

    for (const TString& cone : kConeLabels) {
        std::cout << "\n=== " << cone << " ===\n";

        std::cout << "  eta symmetry...\n";
        PlotEtaSym(fIn, outDir, cone, bins);

        std::cout << "  method comparison (|eta|)...\n";
        PlotMethodComp(fIn, outDir, cone, bins, false);

        std::cout << "  method comparison (full eta)...\n";
        PlotMethodComp(fIn, outDir, cone, bins, true);
    }

    fIn->Close();
    std::cout << "\nDone. Plots in: " << outDir << "\n";
}

#ifndef __CLING__
#include <iostream>
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: plotResiduals <residuals.root> [out_dir]\n";
        return 1;
    }
    plotResiduals(argv[1], argc >= 3 ? argv[2] : "");
    return 0;
}
#endif
