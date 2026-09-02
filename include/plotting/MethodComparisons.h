#ifndef L2RESIDUALS_PLOTTING_METHOD_COMPARISONS_H
#define L2RESIDUALS_PLOTTING_METHOD_COMPARISONS_H

#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"

#include "plotting/Style.h"
#include "plotting/Utilities.h"
#include "Binning.h"
#include "Naming.h"
#include "ProgressBar.h"

#include <algorithm>
#include <vector>

// method comparison (gauss vs trunc90 vs trunc95). One canvas per (cone,
// ptavg slice, eta type): top = all three overlaid, bottom = trunc/gauss ratios.

inline void PlotMethodComp(TFile *fIn, const TString &outDir,
                           const TString &cone, const BinningConfig &bins,
                           bool fullEta, ProgressBar &pb, bool useJer = false) {
  const TString etaMode = L2Name::EtaModeKey(fullEta);
  const double xMin =
      fullEta ? kEtaEdges.front() : (double)kAbsEtaEdges.front();
  const double xMax = fullEta ? kEtaEdges.back() : (double)kAbsEtaEdges.back();
  const TString xTitle = fullEta ? "#eta" : "|#eta|";
  const TString calibKey = useJer ? "jer" : "jec";

  for (int ip = 0; ip < (int)bins.ptavgSlices.size(); ip++) {
    if (!pb.ShouldKeep())
      continue;
    const auto &sl = bins.ptavgSlices[ip];

    std::vector<TH1D *> hists(kNMethods, nullptr);
    for (int m = 0; m < kNMethods; m++) {
      TString name =
          L2Name::ObjectName(cone, CalibKind("intercept", useJer),
                             {etaMode, L2Name::PtKey(sl)}, {kMethodKeys[m]});
      hists[m] = GetHAny(fIn, {cone + "/" + name});
    }

    if (!hists[0]) {
      for (auto *h : hists)
        delete h;
      continue;
    }

    // ratios to gauss for non-null trunc methods
    std::vector<TH1D *> ratios;
    std::vector<int> ratioIdx;
    for (int m = 1; m < kNMethods; m++) {
      if (!hists[m])
        continue;
      TString rname = Form("%s_mcomp%s%s_r%d", cone.Data(), sl.shortName.Data(),
                           etaMode.Data(), m);
      ratios.push_back(RatioH(hists[m], hists[0], rname));
      ratioIdx.push_back(m);
    }

    TString ptKey = L2Name::PtKey(sl);
    const TString cvName = Form("methods_%s_%s_%s_%s", cone.Data(),
                                calibKey.Data(), etaMode.Data(), ptKey.Data());
    TwoPad cv = MakeTwoPad(cvName);

    // main pad
    cv.main->cd();
    cv.main->SetGridx();
    cv.main->SetGridy();

    // y starts fixed at 0.9, tops out just above the real max instead of a
    // YRange() guess over only the first 3 methods
    double ylo = 0.9, yhi = -1e9;
    bool anyMax = false;
    for (auto *h : hists) {
      if (!h)
        continue;
      for (int i = 1; i <= h->GetNbinsX(); i++) {
        const double v = h->GetBinContent(i), e = h->GetBinError(i);
        if (v == 0 && e == 0)
          continue;
        yhi = std::max(yhi, v + e);
        anyMax = true;
      }
    }
    yhi = anyMax ? yhi + 0.10 * (yhi - ylo) : 1.3;

    auto [occXMin, occXMax] = OccupiedRangeWithMargin(hists);
    double plotXMin = xMin, plotXMax = xMax;
    if (occXMin < occXMax) {
      plotXMin = occXMin;
      plotXMax = occXMax;
    }

    bool first = true;

    // one legend (method markers + cone/pT/calib-tag info, |eta|/full-eta
    // dropped), upper-left, clear of the curve's low-|eta| flat region
    TLegend *leg = new TLegend(0.16, 0.60, 0.56, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.038);
    leg->AddEntry((TObject *)nullptr, cone, "");
    leg->AddEntry((TObject *)nullptr, sl.title, "");
    // this repo calls the mean-derived Step 2 output "JEC" everywhere else,
    // but the legend here is meant to read as the correction product itself
    leg->AddEntry((TObject *)nullptr, useJer ? "JER SF" : "L2Residual", "");

    for (int m = 0; m < kNMethods; m++) {
      if (!hists[m])
        continue;
      StyleH(hists[m], kMethodColors[m], kMethodStyles[m], 2.f);
      hists[m]->GetXaxis()->SetRangeUser(plotXMin, plotXMax);
      hists[m]->GetYaxis()->SetRangeUser(ylo, yhi);
      hists[m]->GetYaxis()->SetTitle(CalibYTitle(useJer));
      hists[m]->GetYaxis()->SetTitleSize(0.055);
      hists[m]->GetYaxis()->SetTitleOffset(1.10);
      hists[m]->GetYaxis()->SetLabelSize(0.040);
      hists[m]->GetXaxis()->SetLabelSize(0.0);
      hists[m]->GetXaxis()->SetTitle("");
      hists[m]->SetTitle("");
      hists[m]->Draw(first ? "E1" : "E1 same");
      first = false;
      leg->AddEntry(hists[m], kMethodLabels[m], "lp");
    }

    RefLine(cv.main, plotXMin, plotXMax, 1.0);
    leg->Draw();
    DrawCMSInternalHeader(0.13, 0.96, 0.915, 0.040);

    // ratio pad
    cv.ratio->cd();
    cv.ratio->SetGridx();
    cv.ratio->SetGridy();

    bool firstR = true;
    for (int k = 0; k < (int)ratios.size(); k++) {
      const int m = ratioIdx[k];
      StyleH(ratios[k], kMethodColors[m], kMethodStyles[m], 1.5f);
      ratios[k]->GetXaxis()->SetRangeUser(plotXMin, plotXMax);
      TuneRatio(ratios[k], xTitle, "/ Gauss", 0.993, 1.007);
      // main/ratio pads are 0.69/0.30 of the canvas height (MakeTwoPad)
      ratios[k]->GetYaxis()->SetLabelSize(0.092);
      ratios[k]->GetXaxis()->SetLabelSize(0.092);
      ratios[k]->Draw(firstR ? "E1" : "E1 same");
      firstR = false;
    }
    if (!firstR) {
      RefLine(cv.ratio, plotXMin, plotXMax, 1.0);
    }

    cv.c->cd();
    SavePlot(cv.c, outDir, cone, "methods", {calibKey, etaMode, ptKey}, cvName);
    pb.Update();

    // ratios and hists were drawn on pads, canvas cascade deletes them
    // rleg, leg, tex also drawn, also cascade-deleted
    delete cv.c;
  }
}

#endif
