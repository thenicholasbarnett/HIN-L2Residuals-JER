#ifndef L2RESIDUALS_PLOTTING_STYLE_H
#define L2RESIDUALS_PLOTTING_STYLE_H

#include "TH1D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"

#include "Colors.h"

#include <cmath>

// ============================================================
// Method palette — shared across every plot family that compares
// Gauss / double-Gauss / trunc90 / trunc95.
// ============================================================

static const char *const kMethodKeys[] = {"gauss", "doubleGauss", "trunc90",
                                          "trunc95"};
static const char *const kMethodLabels[] = {
    "Gauss fit", "Double-Gauss fit", "Trunc. mean 90%", "Trunc. mean 95%"};
static const Color_t kMethodColors[] = {HiroshigeNightBlue(), KlimtGreen(),
                                        HiroshigeOrange(), HiroshigeLightRed()};
static const int kMethodStyles[] = {
    20, 34, 21, 22}; // circle / cross / square / triangle-up
static constexpr int kNMethods = 4;

// colors for the eta-symmetry comparison
static const Color_t kColFull = HiroshigeLightRed();
static const Color_t kColRefl = HiroshigeNightBlue();

// ============================================================
// JEC (default, mean-derived) vs JER SF (stddev-derived) selector for
// Step 2 plots reading intercept/R/R_data/R_mc-kind objects -- see
// runCalibration's CALIBRATION=JEC|JER token. Object names differ by a
// "_jer" suffix (CalibrationExtractor.cxx writes both kinds side by side);
// CalibKind builds the right one from the base JEC-kind string.
// ============================================================

inline TString CalibKind(const TString &base, bool useJer) {
  return useJer ? base + "_jer" : base;
}

inline TString CalibTag(bool useJer) { return useJer ? "JER SF" : "JEC"; }

inline TString CalibYTitle(bool useJer) {
  return useJer ? "JER SF at #alpha#rightarrow0"
                : "R_{MC}/R_{data} at #alpha#rightarrow0";
}

// ============================================================
// Global plotting style
// ============================================================

inline void SetupPlotStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetGridColor(kGray + 1);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  gErrorIgnoreLevel = kWarning;
}

// ============================================================
// Histogram / line style helpers
// ============================================================

inline void StyleH(TH1D *h, Color_t col, int mstyle, float lw = 2.0f) {
  h->SetLineColor(col);
  h->SetMarkerColor(col);
  h->SetMarkerStyle(mstyle);
  h->SetMarkerSize(0.85f);
  h->SetLineWidth((Width_t)lw);
}

// Horizontal dashed reference line at y=ref drawn on pad p.
inline void RefLine(TPad *p, double x0, double x1, double ref = 1.0) {
  p->cd();
  TLine *l = new TLine(x0, ref, x1, ref);
  l->SetLineStyle(2);
  l->SetLineColor(kGray + 2);
  l->SetLineWidth(1);
  l->Draw();
}

// Vertical line placeholder — actual y range set by caller after drawing histogram.
inline TLine *VLine(double x, Color_t col, int style = 3) {
  TLine *l = new TLine(x, 0, x, 1);
  l->SetLineColor(col);
  l->SetLineStyle(style);
  l->SetLineWidth(2);
  return l;
}

// Tune axis sizes for the ratio pad (30% of canvas height).
// Text sizes are larger to appear the same physical size as the main pad.
inline void TuneRatio(TH1D *h, const TString &xTitle, const TString &yTitle,
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

// ============================================================
// CMS/internal labels and common text labels
// ============================================================

inline void DrawCMSInternalHeader(double xLeft = 0.13, double xRight = 0.95,
                                  double y = 0.905) {
  TLatex *cms = new TLatex(xLeft, y, "#bf{CMS} #it{Internal}");
  cms->SetNDC();
  cms->SetTextFont(42);
  cms->SetTextAlign(11);
  cms->SetTextSize(0.035);
  cms->Draw();

  TLatex *lumi = new TLatex(xRight, y, "2024 pp (5.36 TeV)");
  lumi->SetNDC();
  lumi->SetTextFont(42);
  lumi->SetTextAlign(31);
  lumi->SetTextSize(0.035);
  lumi->Draw();
}

inline void DrawAsymHeader() { DrawCMSInternalHeader(); }

inline TString FormatEntriesText(Long64_t entries) {
  if (entries < 100)
    return Form("%lld Entries", entries);
  const int exponent = (int)std::floor(std::log10((double)entries));
  const double mantissa = entries / std::pow(10.0, exponent);
  return Form("%.3g #times 10^{%d} Entries", mantissa, exponent);
}

inline void DrawEntriesLabel(Long64_t entries, double x = 0.88,
                             double y = 0.84) {
  TLatex *tex = new TLatex(x, y, FormatEntriesText(entries));
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextAlign(31);
  tex->SetTextSize(0.035);
  tex->Draw();
}

#endif
