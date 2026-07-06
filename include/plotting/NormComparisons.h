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
    StyleH(dummyN, kGray + 2, 24, 1.5f);

    TLegend *legStyle = new TLegend(0.16, 0.14, 0.46, 0.25);
    legStyle->SetBorderSize(0);
    legStyle->SetFillStyle(0);
    legStyle->SetTextSize(0.046);
    legStyle->AddEntry(dummyD, "Direct", "lp");
    legStyle->AddEntry(dummyN, "kFSR norm.", "lp");

    TLegend *legPt = new TLegend(0.60, 0.88 - 0.055 * nPt, 0.94, 0.88);
    legPt->SetBorderSize(0);
    legPt->SetFillStyle(0);
    legPt->SetTextSize(0.046);

    bool first = true;
    for (int ip = 0; ip < nPt; ip++) {
      if (!hDirect[ip])
        continue;
      const auto &ptSl = bins.ptavgSlices[ip];
      StyleH(hDirect[ip], ptSl.color, 20, 1.5f);
      hDirect[ip]->GetYaxis()->SetRangeUser(ylo, yhi);
      hDirect[ip]->GetYaxis()->SetTitle(CalibYTitle(useJer));
      hDirect[ip]->GetYaxis()->SetTitleSize(0.052);
      hDirect[ip]->GetYaxis()->SetTitleOffset(1.15);
      hDirect[ip]->GetYaxis()->SetLabelSize(0.048);
      hDirect[ip]->GetXaxis()->SetLabelSize(0.0);
      hDirect[ip]->GetXaxis()->SetTitle("");
      hDirect[ip]->SetTitle("");
      hDirect[ip]->Draw(first ? "E1" : "E1 same");
      first = false;
      if (hNorm[ip]) {
        StyleH(hNorm[ip], ptSl.color, 24, 1.5f);
        hNorm[ip]->SetTitle("");
        hNorm[ip]->Draw("E1 same");
      }
      legPt->AddEntry(hDirect[ip], ptSl.title, "lp");
    }
    RefLine(cv.main, xMin, xMax, 1.0);
    legStyle->Draw();
    legPt->Draw();

    // legStyle (bottom-left) and legPt (top-right) already occupy those
    // corners -- top-left is the free one
    DrawInfoLegend(0.16, 0.68, 0.50, 0.90,
                   {cone, kMethodLabels[m], xTitle, CalibTag(useJer)});

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
      TuneRatio(hRatios[ip], xTitle, "Norm / Direct", rlo, rhi);
      hRatios[ip]->Draw(firstR ? "E1" : "E1 same");
      firstR = false;
    }
    if (!firstR)
      RefLine(cv.ratio, xMin, xMax, 1.0);

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
