#ifndef L2RESIDUALS_PLOTTING_NORM_COMPARISONS_H
#define L2RESIDUALS_PLOTTING_NORM_COMPARISONS_H

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

#include <vector>

// Direct vs kFSR-norm correction factor overlay. Two-panel canvas per
// (cone, method, eta mode): top = direct (filled) vs kFSR-norm (open) vs eta,
// all pT slices overlaid; bottom = norm/direct ratio per slice.

inline void PlotNormComp(TFile *fIn, const TString &outDir, const TString &cone,
                         const BinningConfig &bins, bool fullEta,
                         ProgressBar &pb, bool useJer = false) {
  const TString etaMode = L2Name::EtaModeKey(fullEta);
  const TString xTitle = fullEta ? "#eta" : "|#eta|";
  const double xMin =
      fullEta ? kEtaEdges.front() : (double)kAbsEtaEdges.front();
  const double xMax = fullEta ? kEtaEdges.back() : (double)kAbsEtaEdges.back();
  const TString calibKey = useJer ? "jer" : "jec";
  const int nPt = (int)bins.ptavgSlices.size();
  for (int m = 0; m < kNMethods; m++) {
    if (!pb.ShouldKeep())
      continue;
    std::vector<TH1D *> hDirect(nPt, nullptr);
    std::vector<TH1D *> hNorm(nPt, nullptr);

    for (int ip = 0; ip < nPt; ip++) {
      const auto &sl = bins.ptavgSlices[ip];
      TString ptKey = L2Name::PtKey(sl);
      TString nameDirect =
          L2Name::ObjectName(cone, CalibKind("intercept", useJer),
                             {etaMode, ptKey}, {kMethodKeys[m]});
      TString nameNorm = nameDirect + "_norm";
      hDirect[ip] = GetHAny(fIn, {cone + "/" + nameDirect});
      hNorm[ip] = GetHAny(fIn, {cone + "/" + nameNorm});
    }

    bool anyValid = false;
    for (int ip = 0; ip < nPt; ip++)
      if (hDirect[ip] && hNorm[ip]) {
        anyValid = true;
        break;
      }

    if (!anyValid) {
      for (auto *h : hDirect)
        delete h;
      for (auto *h : hNorm)
        delete h;
      continue;
    }

    // ratio: norm / direct
    std::vector<TH1D *> hRatios(nPt, nullptr);
    for (int ip = 0; ip < nPt; ip++) {
      if (!hDirect[ip] || !hNorm[ip])
        continue;
      TString rn = Form("normcomp_%s_%s_%s_r%d", cone.Data(), etaMode.Data(),
                        kMethodKeys[m], ip);
      hRatios[ip] = RatioH(hNorm[ip], hDirect[ip], rn);
    }

    const TString cvName =
        Form("normcomp_%s_%s_%s_%s", cone.Data(), calibKey.Data(),
             kMethodKeys[m], etaMode.Data());
    TwoPad cv = MakeTwoPad(cvName);

    // main pad
    cv.main->cd();
    cv.main->SetGridx();
    cv.main->SetGridy();

    // combined y-range over direct + norm
    std::vector<TH1D *> allMain;
    for (int ip = 0; ip < nPt; ip++) {
      if (hDirect[ip])
        allMain.push_back(hDirect[ip]);
      if (hNorm[ip])
        allMain.push_back(hNorm[ip]);
    }
    auto [ylo, yhi] = YRange(allMain);
    auto [occXMin, occXMax] = OccupiedRangeWithMargin(allMain);
    double plotXMin = xMin, plotXMax = xMax;
    if (occXMin < occXMax) {
      plotXMin = occXMin;
      plotXMax = occXMax;
    }

    // find first valid histogram for cloning dummy legend entries
    TH1D *firstD = nullptr;
    for (int ip = 0; ip < nPt; ip++)
      if (hDirect[ip]) {
        firstD = hDirect[ip];
        break;
      }

    TH1D *dummyD = (TH1D *)firstD->Clone("_nc_dd");
    dummyD->SetDirectory(0);
    TH1D *dummyN = (TH1D *)firstD->Clone("_nc_dn");
    dummyN->SetDirectory(0);
    StyleH(dummyD, kGray + 2, 20, 1.5f);
    StyleH(dummyN, kGray + 2, 25, 2.0f); // open square + dotted errors, distinct
    dummyN->SetLineStyle(3);           // from Direct's filled circle/solid line

    // one legend (cone/method/calib-tag info, Direct/kFSR-norm marker key,
    // pT-slice color key), upper-middle of the plot -- |eta| vs full-eta
    // dropped, it's already the x-axis title
    const int nLegEntries = 3 + 2 + nPt;
    const double legX1 = 0.32, legX2 = 0.68;
    const double legY1 = 0.55;
    const double legY2 = legY1 + 0.026 * (double)nLegEntries;
    TLegend *leg = new TLegend(legX1, legY1, legX2, legY2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.030);
    leg->AddEntry((TObject *)nullptr, cone, "");
    leg->AddEntry((TObject *)nullptr, kMethodLabels[m], "");
    leg->AddEntry((TObject *)nullptr, CalibTag(useJer), "");
    leg->AddEntry(dummyD, "Direct", "lp");
    leg->AddEntry(dummyN, "kFSR norm.", "lp");

    bool first = true;
    for (int ip = 0; ip < nPt; ip++) {
      if (!hDirect[ip])
        continue;
      const auto &ptSl = bins.ptavgSlices[ip];
      StyleH(hDirect[ip], ptSl.color, 20, 1.5f);
      hDirect[ip]->GetXaxis()->SetRangeUser(plotXMin, plotXMax);
      hDirect[ip]->GetYaxis()->SetRangeUser(ylo, yhi);
      hDirect[ip]->GetYaxis()->SetTitle(useJer ? "JER SF"
                                               : "L2Residual Corr. Factor");
      hDirect[ip]->GetYaxis()->SetTitleSize(0.048);
      hDirect[ip]->GetYaxis()->SetTitleOffset(1.10);
      hDirect[ip]->GetYaxis()->SetLabelSize(0.040);
      hDirect[ip]->GetXaxis()->SetLabelSize(0.0);
      hDirect[ip]->GetXaxis()->SetTitle("");
      hDirect[ip]->SetTitle("");
      // Direct stays a filled circle with solid error bars; kFSR norm gets
      // its own marker shape (open square) and dotted error bars/bin-width
      // ticks so it reads as a genuinely different series instead of
      // vanishing under a same-color, same-position filled point
      hDirect[ip]->Draw(first ? "E1" : "E1 same");
      first = false;
      if (hNorm[ip]) {
        StyleH(hNorm[ip], ptSl.color, 25, 2.0f);
        hNorm[ip]->SetLineStyle(3);
        hNorm[ip]->SetTitle("");
        hNorm[ip]->Draw("E1 same");
      }
      leg->AddEntry(hDirect[ip], ptSl.title, "lp");
    }
    RefLine(cv.main, plotXMin, plotXMax, 1.0);
    leg->Draw();
    DrawCMSInternalHeader(0.13, 0.96, 0.915, 0.040);

    // ratio pad
    cv.ratio->cd();
    cv.ratio->SetGridx();
    cv.ratio->SetGridy();

    auto [rlo, rhi] = YRange(hRatios);
    if (rlo > rhi) {
      rlo = 0.985;
      rhi = 1.015;
    }

    bool firstR = true;
    for (int ip = 0; ip < nPt; ip++) {
      if (!hRatios[ip])
        continue;
      StyleH(hRatios[ip], bins.ptavgSlices[ip].color, 20, 1.5f);
      hRatios[ip]->GetXaxis()->SetRangeUser(plotXMin, plotXMax);
      TuneRatio(hRatios[ip], xTitle, "Norm / Direct", rlo, rhi);
      // main/ratio pads are 0.69/0.30 of the canvas height (MakeTwoPad)
      hRatios[ip]->GetYaxis()->SetLabelSize(0.092);
      hRatios[ip]->GetXaxis()->SetLabelSize(0.092);
      hRatios[ip]->Draw(firstR ? "E1" : "E1 same");
      firstR = false;
    }
    if (!firstR)
      RefLine(cv.ratio, plotXMin, plotXMax, 1.0);

    cv.c->cd();
    SavePlot(cv.c, outDir, cone, "normcomp", {calibKey, etaMode}, cvName);
    pb.Update();

    delete dummyD;
    delete dummyN;
    // hDirect and hNorm entries where one is null were not drawn, delete explicitly
    for (int ip = 0; ip < nPt; ip++) {
      if (!hDirect[ip] && hNorm[ip])
        delete hNorm[ip];
    }
    // all drawn objects (hDirect, hNorm, hRatios, legends, tex) cascade-deleted with canvas
    delete cv.c;
  }
}

#endif
