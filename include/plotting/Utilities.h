#ifndef L2RESIDUALS_PLOTTING_UTILITIES_H
#define L2RESIDUALS_PLOTTING_UTILITIES_H

#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TDatime.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

#include "Binning.h"
#include "Naming.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

// canvas / pad helpers

inline void RealAspectRatio(TCanvas *c) {
  if (!gROOT->IsBatch())
    c->SetRealAspectRatio();
}

struct TwoPad {
  TCanvas *c = nullptr;
  TPad *main = nullptr;
  TPad *ratio = nullptr;
};

// Main panel: y [0.30, 1.0].  Ratio panel: y [0.00, 0.32].
// Pads are added to the canvas's primitive list; delete canvas to cascade-delete both.
inline TwoPad MakeTwoPad(const TString &name) {
  TwoPad cv;
  cv.c = new TCanvas(name, "", 800, 800);
  RealAspectRatio(cv.c);
  cv.c->cd();

  cv.main = new TPad(name + "_m", "", 0.0, 0.30, 1.0, 1.00);
  cv.main->SetBottomMargin(0.025);
  cv.main->SetTopMargin(0.10);
  cv.main->SetLeftMargin(0.13);
  cv.main->SetRightMargin(0.04);
  cv.main->Draw();

  cv.ratio = new TPad(name + "_r", "", 0.0, 0.00, 1.0, 0.32);
  cv.ratio->SetTopMargin(0.05);
  cv.ratio->SetBottomMargin(0.33);
  cv.ratio->SetLeftMargin(0.13);
  cv.ratio->SetRightMargin(0.04);
  cv.ratio->Draw();

  return cv;
}

// ROOT object retrieval / cloning

// Clone h from file, disassociate from directory so the canvas can own it cleanly.
inline TH1D *GetH(TFile *f, const TString &name) {
  TH1D *src = (TH1D *)f->Get(name);
  if (!src)
    return nullptr;
  TH1D *h = (TH1D *)src->Clone(name + "_c");
  h->SetDirectory(0);
  return h;
}

inline TH1D *GetHAny(TFile *f, const std::vector<TString> &names) {
  for (const auto &name : names) {
    TH1D *h = GetH(f, name);
    if (h)
      return h;
  }
  return nullptr;
}

inline bool HasH(TFile *f, const TString &name) { return f && f->Get(name); }

inline bool HasH(TDirectory *d, const TString &name) {
  return d && d->Get(name);
}

inline bool HasHAny(TFile *f, const std::vector<TString> &names) {
  for (const auto &name : names)
    if (HasH(f, name))
      return true;
  return false;
}

inline bool HasGraphWithN(TDirectory *d, const std::vector<TString> &names,
                          int minN) {
  if (!d)
    return false;
  for (const auto &name : names) {
    TGraphErrors *g = (TGraphErrors *)d->Get(name);
    if (g && g->GetN() >= minN)
      return true;
  }
  return false;
}

inline TH1D *GetHAny(TDirectory *d, const std::vector<TString> &names) {
  if (!d)
    return nullptr;
  for (const auto &name : names) {
    TH1D *src = (TH1D *)d->Get(name);
    if (!src)
      continue;
    TH1D *h = (TH1D *)src->Clone(name + "_c");
    h->SetDirectory(0);
    return h;
  }
  return nullptr;
}

inline TGraphErrors *GetGraphAny(TDirectory *d,
                                 const std::vector<TString> &names) {
  if (!d)
    return nullptr;
  for (const auto &name : names) {
    TGraphErrors *g = (TGraphErrors *)d->Get(name);
    if (g)
      return g;
  }
  return nullptr;
}

// A TGraphErrors built by CalibrationExtractor.cxx/TextFileWriter.cxx carries
// at most one embedded fit function; this is that function, if any.
inline TF1 *GetGraphFit(TGraphErrors *g) {
  if (!g || !g->GetListOfFunctions() || g->GetListOfFunctions()->GetSize() == 0)
    return nullptr;
  return (TF1 *)g->GetListOfFunctions()->At(0);
}

// Clone a TH2D from file, disassociate from directory. 2D analog of GetHAny(TFile*, ...).
inline TH2D *GetH2Any(TFile *f, const std::vector<TString> &names) {
  for (const auto &name : names) {
    TH2D *src = (TH2D *)f->Get(name);
    if (!src)
      continue;
    TH2D *h = (TH2D *)src->Clone(name + "_c");
    h->SetDirectory(0);
    return h;
  }
  return nullptr;
}

// histogram math helpers

// Mirror an |eta| TH1D (kAbsEtaEdges, 18 bins) onto a full-eta TH1D (kEtaEdges, 36 bins).
// Full-eta bin i (1-indexed): i <= 18 → |eta| bin (18-i+1),  i > 18 → |eta| bin (i-18).
inline TH1D *Reflect(TH1D *hAbs, const TString &name) {
  const int n = hAbs->GetNbinsX(); // 18
  TH1D *h = new TH1D(name, "", (int)kEtaEdges.size() - 1, kEtaEdges.data());
  h->SetDirectory(0);
  for (int i = 1; i <= 2 * n; i++) {
    const int a = (i <= n) ? (n - i + 1) : (i - n);
    h->SetBinContent(i, hAbs->GetBinContent(a));
    h->SetBinError(i, hAbs->GetBinError(a));
  }
  return h;
}

// h1 / h2 bin-by-bin; leaves bins zero where either is zero.
inline TH1D *RatioH(TH1D *h1, TH1D *h2, const TString &name) {
  TH1D *r = (TH1D *)h1->Clone(name);
  r->SetDirectory(0);
  r->Reset();
  const int n = r->GetNbinsX();
  for (int i = 1; i <= n; i++) {
    const double v1 = h1->GetBinContent(i), e1 = h1->GetBinError(i);
    const double v2 = h2->GetBinContent(i), e2 = h2->GetBinError(i);
    if (std::abs(v2) < 1e-9 || std::abs(v1) < 1e-9)
      continue;
    const double rv = v1 / v2;
    r->SetBinContent(i, rv);
    r->SetBinError(i, rv * std::hypot(e1 / v1, e2 / v2));
  }
  return r;
}

// Auto y-range across a set of histograms, ignoring empty bins, with padding.
inline std::pair<double, double> YRange(const std::vector<TH1D *> &hv,
                                        double pad = 0.15) {
  double lo = 1e9, hi = -1e9;
  for (auto *h : hv) {
    if (!h)
      continue;
    for (int i = 1; i <= h->GetNbinsX(); i++) {
      const double v = h->GetBinContent(i), e = h->GetBinError(i);
      if (v == 0 && e == 0)
        continue;
      lo = std::min(lo, v - e);
      hi = std::max(hi, v + e);
    }
  }
  if (lo > hi)
    return {0.88, 1.12};
  const double span = hi - lo;
  return {lo - pad * span, hi + pad * span};
}

// TH2D analog of YRange -- structural cells with no data (e.g. a pT slice
// with no jet reach at a given eta) sit at exactly 0 and otherwise drag
// COLZ's default [0, max] autorange down to 0, washing out a real spread
// clustered far from zero (e.g. JER SF values near 1.0).
inline std::pair<double, double> ZRange(TH2D *h, double pad = 0.15) {
  double lo = 1e9, hi = -1e9;
  for (int ix = 1; ix <= h->GetNbinsX(); ix++) {
    for (int iy = 1; iy <= h->GetNbinsY(); iy++) {
      const double v = h->GetBinContent(ix, iy), e = h->GetBinError(ix, iy);
      if (v == 0 && e == 0)
        continue;
      lo = std::min(lo, v);
      hi = std::max(hi, v);
    }
  }
  if (lo > hi)
    return {0.0, 1.0};
  const double span = hi - lo;
  if (span < 1e-9)
    return {lo - 0.05 * std::max(std::fabs(lo), 1.0),
            hi + 0.05 * std::max(std::fabs(hi), 1.0)};
  return {lo - pad * span, hi + pad * span};
}

// Widest occupied bin range across a set of histograms sharing the same
// axis binning, extended by one empty bin on each side that has content --
// so the axis ends just past the real data instead of running out to the
// full binning's extent when the data itself stops well short of it (e.g.
// a central pT slice with limited jet reach in eta). Returns {0,0} if none
// of the histograms have any content.
//
// yLo/yHi optionally restrict "occupied" to bins whose value actually falls
// inside that window -- for a fixed, tight y-range (closures) a point can
// sit well outside the visible band while still counting as occupied,
// stretching the x-axis out past where anything is actually drawn. Defaults
// span the full real line, so callers with an auto-scaled y-range see no
// change in behavior.
inline std::pair<double, double> OccupiedRangeWithMargin(
    const std::vector<TH1D *> &hv,
    double yLo = -std::numeric_limits<double>::infinity(),
    double yHi = std::numeric_limits<double>::infinity()) {
  int loBin = -1, hiBin = -1;
  TH1D *ref = nullptr;
  for (auto *h : hv) {
    if (!h)
      continue;
    if (!ref)
      ref = h;
    for (int i = 1; i <= h->GetNbinsX(); i++) {
      const double v = h->GetBinContent(i), e = h->GetBinError(i);
      if ((v == 0 && e == 0) || v < yLo || v > yHi)
        continue;
      if (loBin < 0 || i < loBin)
        loBin = i;
      if (i > hiBin)
        hiBin = i;
    }
  }
  if (!ref || loBin < 0)
    return {0.0, 0.0};
  const int nBins = ref->GetNbinsX();
  const int loMargin = std::max(1, loBin - 1);
  const int hiMargin = std::min(nBins, hiBin + 1);
  return {ref->GetXaxis()->GetBinLowEdge(loMargin),
          ref->GetXaxis()->GetBinUpEdge(hiMargin)};
}

// plot path helpers

inline TString PlotDir(const TString &outDir, const TString &cone,
                       const TString &plotType,
                       const std::vector<TString> &orderedKeys = {}) {
  std::vector<TString> parts = {outDir, cone, plotType};
  parts.insert(parts.end(), orderedKeys.begin(), orderedKeys.end());
  TString dir = L2Name::Join(parts, "/");
  if (gSystem->mkdir(dir, true) < 0 && gSystem->AccessPathName(dir))
    std::cerr << "Warning: cannot create plot subdirectory: " << dir << "\n";
  return dir;
}

inline void SavePlot(TCanvas *c, const TString &outDir, const TString &cone,
                     const TString &plotType,
                     const std::vector<TString> &orderedKeys,
                     const TString &fileBase) {
  TString dir = PlotDir(outDir, cone, plotType, orderedKeys);
  c->SaveAs(Form("%s/%s.png", dir.Data(), fileBase.Data()));
}

// Small borderless/fill-less legend of label-only lines. Replaces the old
// convention of one oversized bold TLatex line concatenating cone/method/
// slice/etc with "|" -- that line ran off the canvas edge for anything with
// more than two or three fields, since a single TLatex has no wrapping or
// bounds checking the way a TLegend's box does.
inline TLegend *DrawInfoLegend(double x1, double y1, double x2, double y2,
                               const std::vector<TString> &lines,
                               double textSize = 0.032) {
  TLegend *leg = new TLegend(x1, y1, x2, y2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(textSize);
  for (const auto &line : lines)
    leg->AddEntry((TObject *)nullptr, line.Data(), "");
  leg->Draw();
  return leg;
}

inline TString MakePlotDir(const TString &prefix = "plots") {
  TDatime now;
  TString timestamp =
      Form("%04d-%02d-%02d_%02d-%02d-%02d", now.GetYear(), now.GetMonth(),
           now.GetDay(), now.GetHour(), now.GetMinute(), now.GetSecond());
  TString dir = prefix + "_" + timestamp;
  gSystem->mkdir(dir, true);
  return dir;
}

// plotting-presentation helpers transplanted from top-level Utilities.h

struct PlotConfig {
  TString runNumber = "";
  TString globalTag = "";
  TString jetAlgo = "";

  float xmin = 15.0;
  float xmax = 200.0;
  float ymax = 1.1;
  float ymin = 0.0;

  int canvasSize = 2400;
  float imageScaling = 1.0;

  float legfrac = 0.2;
};

inline TCanvas *MakeSinglePadCanvas(const TString &name, const PlotConfig &cfg,
                                    bool grid = false) {
  TCanvas *c = new TCanvas(name, "", cfg.canvasSize, cfg.canvasSize);
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.12);
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.05);
  if (!gROOT->IsBatch())
    c->SetRealAspectRatio();
  if (grid) {
    c->SetGridx();
    c->SetGridy();
  }
  return c;
}

inline TCanvas *MakeColzCanvas(const TString &name, const PlotConfig &cfg) {
  TCanvas *c = new TCanvas(name, "", cfg.canvasSize, cfg.canvasSize);
  c->SetTopMargin(0.08);
  c->SetBottomMargin(0.12);
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.15);
  if (!gROOT->IsBatch())
    c->SetRealAspectRatio();
  return c;
}

struct SplitCanvas {
  TCanvas *c = nullptr;
  TPad *plotpad = nullptr;
  TPad *legpad = nullptr;
};

inline SplitCanvas MakeSplitPadCanvas(const TString &name,
                                      const PlotConfig &cfg) {
  SplitCanvas sc;
  sc.c = new TCanvas(name, "", cfg.canvasSize, cfg.canvasSize);
  RealAspectRatio(sc.c);
  sc.c->cd();

  sc.legpad = new TPad(name + "_legpad", "", 0.0, 1.0 - cfg.legfrac, 1.0, 1.0);
  sc.legpad->SetTopMargin(0.05);
  sc.legpad->SetBottomMargin(0.0);
  sc.legpad->SetLeftMargin(0.0);
  sc.legpad->SetRightMargin(0.0);
  sc.legpad->SetFillStyle(0);
  sc.legpad->Draw();

  sc.plotpad =
      new TPad(name + "_plotpad", "", 0.0, 0.0, 1.0, 1.0 - cfg.legfrac);
  sc.plotpad->SetTopMargin(0.02);
  sc.plotpad->SetBottomMargin(0.12);
  sc.plotpad->SetLeftMargin(0.12);
  sc.plotpad->SetRightMargin(0.05);
  sc.plotpad->SetGridx();
  sc.plotpad->SetGridy();
  sc.plotpad->Draw();

  return sc;
}

inline TLegend *MakeLegend(float xmin = 0.575, float ymin = 0.75,
                           float xmax = 0.8, float ymax = 0.95) {
  TLegend *l = new TLegend(xmin, ymin, xmax, ymax);
  l->SetTextSize(0);
  l->SetTextAlign(22);
  l->SetMargin(0);
  l->SetFillColor(kWhite);
  l->SetFillStyle(1001);
  l->SetLineColor(kBlack);
  l->SetLineWidth(3);
  return l;
}

inline void AddInfoEntries(TLegend *l, const PlotConfig &cfg) {
  if (!cfg.runNumber.IsNull())
    l->AddEntry((TObject *)nullptr, "Run " + cfg.runNumber, "");
  if (!cfg.globalTag.IsNull())
    l->AddEntry((TObject *)nullptr, cfg.globalTag, "");
  if (!cfg.jetAlgo.IsNull())
    l->AddEntry((TObject *)nullptr, cfg.jetAlgo, "");
}

inline void SetPlotStyle(const PlotConfig &cfg) {
  gStyle->SetImageScaling(cfg.imageScaling);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.04, "XY");
  gStyle->SetLabelSize(0.04, "XY");
  gStyle->SetGridColor(kGray + 2);
  gStyle->SetGridWidth(1);
  gStyle->SetGridStyle(3);
  gStyle->SetPalette(kBird);
}

inline void StyleTH1(TH1 *h, Color_t color) {
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.8);
  h->SetLineWidth(4);
}

inline void DrawLabel(const TString &text, float x, float y,
                      float textSize = 0.035, Int_t align = 11) {
  if (text.IsNull())
    return;
  TLatex *tex = new TLatex(x, y, text);
  tex->SetNDC();
  tex->SetTextAlign(align);
  tex->SetTextSize(textSize);
  tex->Draw();
}

inline void DrawCMSLabel(const TString &cms_type = "", float x = 0.13,
                         float y = 0.965, float textSize = 0.035) {
  DrawLabel(Form("CMS #bf{#it{%s}}", cms_type.Data()), x, y, textSize);
}

inline void DrawJetAlgoLabel(const TString &jetAlgo = "akCs4PF", float x = 0.8,
                             float y = 0.04, float textSize = 0.035) {
  DrawLabel(Form("#bf{%s}", jetAlgo.Data()), x, y, textSize);
}

inline void DrawRunNumber(const TString &RunNumber = "", float x = 0.14,
                          float y = 0.9, float textSize = 0.035) {
  DrawLabel(Form("#bf{%s}", RunNumber.Data()), x, y, textSize);
}

#endif
