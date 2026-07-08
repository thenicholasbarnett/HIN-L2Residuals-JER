#ifndef L2RESIDUALS_PLOTTING_RESPONSE_H
#define L2RESIDUALS_PLOTTING_RESPONSE_H

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"

#include "plotting/Style.h"
#include "plotting/Utilities.h"
#include "plotting/AsymmetryDistributions.h" // StyleFit
#include "Binning.h"
#include "Naming.h"
#include "ProgressBar.h"
#include "jetmet/Variables.h" // JME-standard axis titles for the pT_gen binning variable

#include <algorithm>
#include <vector>

// runResponse (JES/JER) QA + summary plots. Three variants per
// (cone, collection = incl/tag/probe): corr/reco/raw -- see ResponseExtractor.h.
//   response_dist:     raw distribution + Gaussian guide fit, refit here
//                       (extraction only writes the raw histogram).
//   response_summary:  JES/JER vs pT_gen, incl/tag/probe overlaid.
//   response_variants: corr/reco/raw overlaid, one per collection.
//   response_etarange: JES/JER vs pT_gen, coarse |eta| detector regions
//                       overlaid, incl jets (corr variant).
// No-op if the input has no response objects.

static const char *const kResponseCollections[] = {"incl", "tag", "probe"};
static constexpr int kNResponseCollections = 3;

static const char *const kResponseVariants[] = {"corr", "reco", "raw"};
static const char *const kResponseVariantLabels[] = {"p_{T}^{corr}/p_{T}^{gen}",
                                                     "p_{T}^{reco}/p_{T}^{gen}",
                                                     "p_{T}^{raw}/p_{T}^{gen}"};
static constexpr int kNResponseVariants = 3;

// pT_gen axis title uses JetMET's refpt convention; response-ratio labels
// stay local (no JME equivalent)
static const TString kPtGenAxisTitle =
    VARIABLES::getVariableAxisTitleString(VARIABLES::refpt, true);

// generic ratio label for plots overlaying all 3 variants at once (the
// legend, not the axis, names which pT is which there)
static const TString kGenericRatioLabel = "p_{T}/p_{T}^{gen}";

// JES = mean of the response ratio, JER = its fractional resolution
// (sigma/mean) -- distinct from runCalibration's data/MC JER scale factor
// (Style.h::CalibYTitle), this is a pure-MC quantity.
inline TString ResponseYTitle(const TString &quantity,
                              const TString &ratioLabel) {
  if (quantity == "JES")
    return "#LT " + ratioLabel + " #GT";
  return "#sigma(" + ratioLabel + ") / #LT " + ratioLabel + " #GT";
}

// below cfg.minJetPt there's no real jet to speak of; start the pT_gen axis
// there instead of at 0 (also keeps the log-x axis away from a zero edge)
inline double ResponsePtGenMin(const TAxis *xaxis) {
  double cutoff = Config().minJetPt;
  return (cutoff > 0.0) ? std::max(cutoff, xaxis->GetXmin())
                        : xaxis->GetXmin();
}

// Two-pass fit: an initial Gaussian within [1-halfWidth, 1+halfWidth] locates
// the real peak, then a second Gaussian refit narrows to +-2 sigma of that
// first pass's own mean -- keeps a single unconstrained pass from chasing
// the long non-Gaussian tail toward the response axis's [0,2] edge.
inline TF1 *FitResponseGuide(TH1D *h, const TString &name, Color_t col,
                             double halfWidth, int minEntries) {
  if (!h || h->GetEntries() < minEntries)
    return nullptr;

  TF1 pass1(name + "_pass1", "gaus", 1.0 - halfWidth, 1.0 + halfWidth);
  pass1.SetParameter(0, h->GetMaximum());
  pass1.SetParameter(1, h->GetMean());
  pass1.SetParameter(2, std::max(h->GetRMS(), 1e-3));
  h->Fit(&pass1, "NQR");

  const double c = pass1.GetParameter(0);
  const double m = pass1.GetParameter(1);
  const double s = std::max(pass1.GetParameter(2), 1e-3);

  TF1 *fit = new TF1(name, "gaus", m - 2.0 * s, m + 2.0 * s);
  fit->SetParameters(c, m, s);
  StyleFit(fit, col);
  h->Fit(fit, "NQSR");
  return fit;
}

inline void DrawResponseDist(TFile *fIn, const TString &outDir,
                             const TString &cone, const TString &dirName,
                             const TString &objName, const TString &label,
                             const std::vector<TString> &plotKeys,
                             const TString &xTitle, double halfWidth,
                             int minEntries, ProgressBar &pb) {
  TDirectory *d = (TDirectory *)fIn->Get(cone + "/" + dirName);
  TH1D *h = GetHAny(d, {objName});
  if (!h || h->GetEntries() < minEntries) {
    if (h)
      delete h;
    return;
  }

  const Long64_t nEntries = (Long64_t)h->GetEntries();
  NormalizeDensity(h); // 1/N dN/d(ratio) -- Scale() leaves GetEntries() alone

  // single series (markers + one fit line): always kBlue/kRed
  StyleH(h, kBlue, 20, 1.5f);
  h->SetTitle("");
  h->GetXaxis()->SetTitle(xTitle);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d(" + xTitle + ")}");
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetTitleOffset(1.25);
  h->GetYaxis()->SetTitleOffset(1.7);

  TString cvName = "response_" + objName;
  TCanvas *c = new TCanvas(cvName, "", 800, 800);
  RealAspectRatio(c);
  c->SetLeftMargin(0.17);
  c->SetLogy();

  h->Draw("E1");

  TF1 *fit =
      FitResponseGuide(h, cvName + "_fit", kRed, halfWidth, minEntries);
  if (fit)
    fit->Draw("same");

  // zoom around the fit peak -- the response axis's hard [0,2] edge leaves a
  // long, sparsely-filled tail that otherwise dominates the drawn range
  double xlo = h->GetXaxis()->GetXmin();
  double xhi = h->GetXaxis()->GetXmax();
  if (fit) {
    const double m = fit->GetParameter(1);
    const double s = std::max(fit->GetParameter(2), 1e-3);
    xlo = std::max(xlo, m - 6.0 * s);
    xhi = std::min(xhi, m + 6.0 * s);
  }
  h->GetXaxis()->SetRangeUser(xlo, xhi);
  h->SetMinimum(0.1); // normalized density -- the far tail below this is noise

  DrawCMSInternalHeader(0.17, 0.90);

  // one legend, hugging the left frame border, upper-left (clear of the
  // Gaussian's peak, which sits centered)
  TLegend *leg = new TLegend(0.17, 0.66, 0.52, 0.855);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.030);
  leg->SetTextFont(42);
  leg->AddEntry((TObject *)nullptr, label.Data(), "");
  leg->AddEntry((TObject *)nullptr, cone.Data(), "");
  leg->AddEntry((TObject *)nullptr, FormatEntriesText(nEntries).Data(), "");
  if (fit) {
    const double chi2Ndf =
        (fit->GetNDF() > 0) ? fit->GetChisquare() / fit->GetNDF() : 0.0;
    leg->AddEntry(fit, Form("#chi^{2}/ndf = %.1f", chi2Ndf), "l");
  }
  leg->Draw();

  SavePlot(c, outDir, cone, "response_dist", plotKeys, cvName);
  pb.Update();

  delete fit;
  delete h;
  delete c;
}

inline void DrawResponseSummary(TFile *fIn, const TString &outDir,
                                const TString &cone,
                                const TString &quantity, // "JES" or "JER"
                                const std::vector<TString> &orderedKeys,
                                const TString &xTitle,
                                const TString &ratioLabel, const TString &tag,
                                ProgressBar &pb) {
  TH1D *h[kNResponseCollections] = {};
  bool any = false;
  for (int ic = 0; ic < kNResponseCollections; ic++) {
    TString name = L2Name::ObjectName(cone, quantity, orderedKeys,
                                      {kResponseCollections[ic]});
    h[ic] = GetHAny(fIn, {cone + "/" + name});
    if (h[ic])
      any = true;
  }
  if (!any)
    return;

  static const Color_t cols[kNResponseCollections] = {
      HiroshigeNightBlue(), HiroshigeBlue(), KlimtRed()};
  static const int markers[kNResponseCollections] = {20, 21, 22};

  TH1D *frame = nullptr;
  for (int ic = 0; ic < kNResponseCollections; ic++)
    if (h[ic]) {
      frame = h[ic];
      break;
    }

  auto [ylo, yhi] = YRange({h[0], h[1], h[2]});
  if (quantity == "JES") {
    ylo = 0.97;
    yhi = 1.02;
  }
  const double xMin = ResponsePtGenMin(frame->GetXaxis());
  const double xMax = frame->GetXaxis()->GetXmax();

  TCanvas *c =
      new TCanvas("response_summary_" + cone + "_" + tag, "", 800, 800);
  RealAspectRatio(c);
  c->SetLogx();
  c->SetLeftMargin(0.14);
  c->SetGridx();
  c->SetGridy();

  frame->SetTitle("");
  frame->GetXaxis()->SetTitle(xTitle);
  frame->GetYaxis()->SetTitle(ResponseYTitle(quantity, ratioLabel));
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitleOffset(1.25);
  frame->GetXaxis()->SetRangeUser(xMin, xMax);
  frame->GetYaxis()->SetRangeUser(ylo, yhi);

  bool first = true;
  TLegend *leg = new TLegend(0.60, 0.68, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.034);
  // cone (the clustering algorithm) rides in the same legend as the
  // incl/tag/probe marker entries instead of its own separate box
  leg->AddEntry((TObject *)nullptr, cone.Data(), "");
  for (int ic = 0; ic < kNResponseCollections; ic++) {
    if (!h[ic])
      continue;
    StyleH(h[ic], cols[ic], markers[ic], 1.5f);
    h[ic]->Draw(first ? "E1" : "E1 same");
    first = false;
    leg->AddEntry(h[ic], kResponseCollections[ic], "lp");
  }
  leg->Draw();

  if (quantity == "JES") {
    TLine *rl = new TLine(xMin, 1.0, xMax, 1.0);
    rl->SetLineStyle(2);
    rl->SetLineColor(kGray + 2);
    rl->SetLineWidth(1);
    rl->Draw();

    TLine *rl99 = new TLine(xMin, 0.99, xMax, 0.99);
    rl99->SetLineStyle(3);
    rl99->SetLineColor(kBlack);
    rl99->SetLineWidth(2);
    rl99->Draw();

    TLine *rl101 = new TLine(xMin, 1.01, xMax, 1.01);
    rl101->SetLineStyle(3);
    rl101->SetLineColor(kBlack);
    rl101->SetLineWidth(2);
    rl101->Draw();
  }

  DrawCMSInternalHeader(0.14, 0.90);

  SavePlot(c, outDir, cone, "response_summary", {tag}, c->GetName());
  pb.Update();

  delete c; // cascade-deletes h[*], leg, lab (and rl if drawn)
}

// overlays corr/reco/raw for one collection: JES/JER vs pT_gen
inline void DrawVariantComparison(TFile *fIn, const TString &outDir,
                                  const TString &cone,
                                  const TString &collection,
                                  const TString &quantity, // "JES" or "JER"
                                  ProgressBar &pb) {
  TH1D *h[kNResponseVariants] = {};
  bool any = false;
  for (int iv = 0; iv < kNResponseVariants; iv++) {
    TString name = L2Name::ObjectName(
        cone, quantity, {kResponseVariants[iv], "vs_ptgen"}, {collection});
    h[iv] = GetHAny(fIn, {cone + "/" + name});
    if (h[iv])
      any = true;
  }
  if (!any)
    return;

  static const Color_t cols[kNResponseVariants] = {HiroshigeNightBlue(),
                                                   HiroshigeBlue(), KlimtRed()};
  static const int markers[kNResponseVariants] = {20, 21, 22};

  TH1D *frame = nullptr;
  for (int iv = 0; iv < kNResponseVariants; iv++)
    if (h[iv]) {
      frame = h[iv];
      break;
    }

  auto [ylo, yhi] = YRange({h[0], h[1], h[2]});
  if (quantity == "JES") {
    ylo = 0.97;
    yhi = 1.02;
  }
  const double xMin = ResponsePtGenMin(frame->GetXaxis());
  const double xMax = frame->GetXaxis()->GetXmax();

  TCanvas *c = new TCanvas("response_variants_" + cone + "_" + collection +
                               "_" + quantity,
                           "", 800, 800);
  RealAspectRatio(c);
  c->SetLogx();
  c->SetLeftMargin(0.14);
  c->SetGridx();
  c->SetGridy();

  frame->SetTitle("");
  frame->GetXaxis()->SetTitle(kPtGenAxisTitle);
  frame->GetYaxis()->SetTitle(ResponseYTitle(quantity, kGenericRatioLabel));
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitleOffset(1.25);
  frame->GetXaxis()->SetRangeUser(xMin, xMax);
  frame->GetYaxis()->SetRangeUser(ylo, yhi);

  bool first = true;
  TLegend *leg = new TLegend(0.68, 0.68, 0.895, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.034);
  // cone (the clustering algorithm) rides in the same legend as the
  // corr/reco/raw marker entries instead of its own separate box
  leg->AddEntry((TObject *)nullptr, cone.Data(), "");
  for (int iv = 0; iv < kNResponseVariants; iv++) {
    if (!h[iv])
      continue;
    StyleH(h[iv], cols[iv], markers[iv], 1.5f);
    h[iv]->Draw(first ? "E1" : "E1 same");
    first = false;
    leg->AddEntry(h[iv], kResponseVariants[iv], "lp");
  }
  leg->Draw();

  if (quantity == "JES") {
    TLine *rl = new TLine(xMin, 1.0, xMax, 1.0);
    rl->SetLineStyle(2);
    rl->SetLineColor(kGray + 2);
    rl->SetLineWidth(1);
    rl->Draw();

    TLine *rl99 = new TLine(xMin, 0.99, xMax, 0.99);
    rl99->SetLineStyle(3);
    rl99->SetLineColor(kBlack);
    rl99->SetLineWidth(2);
    rl99->Draw();

    TLine *rl101 = new TLine(xMin, 1.01, xMax, 1.01);
    rl101->SetLineStyle(3);
    rl101->SetLineColor(kBlack);
    rl101->SetLineWidth(2);
    rl101->Draw();
  }

  DrawCMSInternalHeader(0.14, 0.90);
  DrawInfoLegend(0.16, 0.80, 0.36, 0.885, {collection});

  SavePlot(c, outDir, cone, "response_variants", {collection, quantity},
           c->GetName());
  pb.Update();

  delete c; // cascade-deletes h[*], leg, lab (and rl if drawn)
}

// overlays JES/JER vs pT_gen across the coarse |eta| detector regions
// (kEtaRangeEdges), incl jets only -- the collection with real eta reach
// (tag is barrel-confined, probe limited). Reads runResponse's
// JER_per_etarange corr-variant output.
inline void DrawEtaRangeComparison(TFile *fIn, const TString &outDir,
                                   const TString &cone,
                                   const TString &quantity, // "JES" or "JER"
                                   ProgressBar &pb) {
  const std::vector<RangeBin> slices = BuildEtaRangeSlices();
  const int nSl = (int)slices.size();
  std::vector<TH1D *> h(nSl, nullptr);
  bool any = false;
  for (int is = 0; is < nSl; is++) {
    TString name = L2Name::ObjectName(
        cone, quantity,
        {"corr", "vs_ptgen", L2Name::EtaKey(slices[is].lo, slices[is].hi)},
        {"incl"});
    h[is] = GetHAny(fIn, {cone + "/JER_per_etarange/" + name});
    if (h[is]) {
      any = true;
    }
  }
  if (!any) {
    return;
  }

  static const int markers[] = {20, 21, 22, 23, 33, 34};

  TH1D *frame = nullptr;
  for (int is = 0; is < nSl; is++) {
    if (h[is]) {
      frame = h[is];
      break;
    }
  }

  auto [ylo, yhi] = YRange(h);
  const double xMin = ResponsePtGenMin(frame->GetXaxis());
  const double xMax = frame->GetXaxis()->GetXmax();

  TCanvas *c =
      new TCanvas("response_etarange_" + cone + "_" + quantity, "", 800, 800);
  RealAspectRatio(c);
  c->SetLogx();
  c->SetLeftMargin(0.14);
  c->SetGridx();
  c->SetGridy();

  frame->SetTitle("");
  frame->GetXaxis()->SetTitle(kPtGenAxisTitle);
  frame->GetYaxis()->SetTitle(
      ResponseYTitle(quantity, kResponseVariantLabels[0]));
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitleOffset(1.25);
  frame->GetXaxis()->SetRangeUser(xMin, xMax);
  frame->GetYaxis()->SetRangeUser(ylo, yhi);

  bool first = true;
  TLegend *leg = new TLegend(0.52, 0.66, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.032);
  leg->SetTextFont(42);
  leg->AddEntry((TObject *)nullptr, cone + ", incl jets", "");
  for (int is = 0; is < nSl; is++) {
    if (!h[is]) {
      continue;
    }
    StyleH(h[is], slices[is].color, markers[is % 6], 1.5f);
    h[is]->Draw(first ? "E1" : "E1 same");
    first = false;
    leg->AddEntry(h[is], slices[is].title, "lp");
  }
  leg->Draw();

  if (quantity == "JES") {
    TLine *rl = new TLine(xMin, 1.0, xMax, 1.0);
    rl->SetLineStyle(2);
    rl->SetLineColor(kGray + 2);
    rl->SetLineWidth(1);
    rl->Draw();
  }

  DrawCMSInternalHeader(0.14, 0.90);

  SavePlot(c, outDir, cone, "response_etarange", {quantity}, c->GetName());
  pb.Update();

  delete c; // cascade-deletes h[*], leg (and rl if drawn)
}

inline void PlotResponse(TFile *fIn, const TString &outDir, const TString &cone,
                         double halfWidth, int minEntries, ProgressBar &pb) {
  // sentinel: does this cone have runResponse output at all?
  bool any = false;
  for (int ic = 0; ic < kNResponseCollections && !any; ic++) {
    TString name = L2Name::ObjectName(cone, "JES", {"corr", "vs_ptgen"},
                                      {kResponseCollections[ic]});
    if (HasHAny(fIn, {cone + "/" + name}))
      any = true;
  }
  if (!any)
    return;

  for (int ic = 0; ic < kNResponseCollections; ic++) {
    const TString collection = kResponseCollections[ic];

    // QA dists vs pt_gen; bin edges read from each variant's own JES output,
    // not re-derived
    for (int iv = 0; iv < kNResponseVariants; iv++) {
      const TString variant = kResponseVariants[iv];
      TString jesName =
          L2Name::ObjectName(cone, "JES", {variant, "vs_ptgen"}, {collection});
      TH1D *hJes = GetHAny(fIn, {cone + "/" + jesName});
      if (hJes) {
        for (int ip = 1; ip <= hJes->GetNbinsX(); ip++) {
          if (!pb.ShouldKeep())
            continue;
          const double lo = hJes->GetXaxis()->GetBinLowEdge(ip);
          const double hi = hJes->GetXaxis()->GetBinUpEdge(ip);
          TString objName =
              L2Name::ObjectName(cone, "response_" + variant,
                                 {L2Name::PtGenKey(lo, hi)}, {collection});
          TString label =
              Form("%.0f #leq p_{T}^{gen} < %.0f GeV", lo, hi);
          DrawResponseDist(
              fIn, outDir, cone, "QA_response_ptgen", objName, label,
              {collection, variant, "vs_ptgen", L2Name::PtGenKey(lo, hi)},
              kResponseVariantLabels[iv], halfWidth, minEntries, pb);
        }
        delete hJes;
      }
    }
  }

  // per-variant summary: 3 collections overlaid
  for (int iv = 0; iv < kNResponseVariants; iv++) {
    const TString variant = kResponseVariants[iv];
    DrawResponseSummary(fIn, outDir, cone, "JES", {variant, "vs_ptgen"},
                        kPtGenAxisTitle, kResponseVariantLabels[iv],
                        "JES_" + variant + "_vs_ptgen", pb);
    DrawResponseSummary(fIn, outDir, cone, "JER", {variant, "vs_ptgen"},
                        kPtGenAxisTitle, kResponseVariantLabels[iv],
                        "JER_" + variant + "_vs_ptgen", pb);
  }

  // variant comparison: corr/reco/raw overlaid, one per collection
  for (int ic = 0; ic < kNResponseCollections; ic++) {
    DrawVariantComparison(fIn, outDir, cone, kResponseCollections[ic], "JES",
                          pb);
    DrawVariantComparison(fIn, outDir, cone, kResponseCollections[ic], "JER",
                          pb);
  }

  // detector-region overlay: coarse |eta| slices, incl jets
  DrawEtaRangeComparison(fIn, outDir, cone, "JES", pb);
  DrawEtaRangeComparison(fIn, outDir, cone, "JER", pb);
}

#endif
