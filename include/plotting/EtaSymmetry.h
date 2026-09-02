#ifndef L2RESIDUALS_PLOTTING_ETA_SYMMETRY_H
#define L2RESIDUALS_PLOTTING_ETA_SYMMETRY_H

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

// Eta symmetry check. One canvas per (cone, method, ptavg slice): top =
// full-eta corrections + |eta| reflected, bottom = full/reflected ratio.

inline void PlotEtaSym(TFile *fIn, const TString &outDir, const TString &cone,
                       const BinningConfig &bins, ProgressBar &pb,
                       bool useJer = false) {
  const TString calibKey = useJer ? "jer" : "jec";
  for (int m = 0; m < kNMethods; m++) {
    for (int ip = 0; ip < (int)bins.ptavgSlices.size(); ip++) {
      if (!pb.ShouldKeep())
        continue;
      const auto &sl = bins.ptavgSlices[ip];
      TString ptKey = L2Name::PtKey(sl);

      const TString nameAbs = L2Name::ObjectName(
          cone, CalibKind("intercept", useJer),
          {L2Name::EtaModeKey(false), ptKey}, {kMethodKeys[m]});
      const TString nameFull = L2Name::ObjectName(
          cone, CalibKind("intercept", useJer),
          {L2Name::EtaModeKey(true), ptKey}, {kMethodKeys[m]});
      TH1D *hAbs = GetHAny(fIn, {cone + "/" + nameAbs});
      TH1D *hFull = GetHAny(fIn, {cone + "/" + nameFull});

      if (!hAbs || !hFull) {
        delete hAbs;
        delete hFull;
        continue;
      }

      // Reflect |eta| onto full-eta axis then compute ratio
      TH1D *hRefl = Reflect(hAbs, nameAbs + "_refl");
      TH1D *hRatio = RatioH(hFull, hRefl, nameAbs + "_rat");

      const TString cvName =
          Form("etasym_%s_%s_%s_%s", cone.Data(), calibKey.Data(),
               kMethodKeys[m], ptKey.Data());
      TwoPad cv = MakeTwoPad(cvName);

      // main pad
      cv.main->cd();
      cv.main->SetGridx();
      cv.main->SetGridy();

      // two series only: always kBlue/kRed
      StyleH(hFull, kBlue, 20, 2.f);
      StyleH(hRefl, kRed, 21, 2.f);
      hRefl->SetLineStyle(2);

      auto [ylo, yhi] = YRange({hFull, hRefl});
      auto [xMin, xMax] = OccupiedRangeWithMargin({hFull, hRefl});
      if (xMin >= xMax) {
        xMin = kEtaEdges.front();
        xMax = kEtaEdges.back();
      }
      hFull->GetXaxis()->SetRangeUser(xMin, xMax);
      hFull->GetYaxis()->SetRangeUser(ylo, yhi);
      hFull->GetYaxis()->SetTitle(CalibYTitle(useJer));
      hFull->GetYaxis()->SetTitleSize(0.055);
      hFull->GetYaxis()->SetTitleOffset(1.10);
      hFull->GetYaxis()->SetLabelSize(0.040);
      hFull->GetXaxis()->SetLabelSize(0.0);
      hFull->GetXaxis()->SetTitle("");
      hFull->SetTitle("");

      hFull->Draw("E1");
      hRefl->Draw("E1 same");
      RefLine(cv.main, xMin, xMax, 1.0);
      DrawCMSInternalHeader(0.13, 0.96, 0.915, 0.040); // slightly larger

      // one legend for both the Full/Reflected markers and the cone/
      // method/pT/calib-tag info, upper-center of the plot
      const double legX1 = 0.345, legX2 = 0.745;
      TLegend *leg = new TLegend(legX1, 0.68, legX2, 0.872);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.038);
      leg->AddEntry((TObject *)nullptr, cone, "");
      leg->AddEntry((TObject *)nullptr, kMethodLabels[m], "");
      leg->AddEntry((TObject *)nullptr, sl.title, "");
      leg->AddEntry((TObject *)nullptr, CalibTag(useJer), "");
      leg->AddEntry(hFull, "Full #eta", "lp");
      leg->AddEntry(hRefl, "|#eta| reflected", "lp");
      leg->Draw();

      // ratio pad
      cv.ratio->cd();
      cv.ratio->SetGridx();
      cv.ratio->SetGridy();

      StyleH(hRatio, kBlack, 20, 1.5f);
      TuneRatio(hRatio, "#eta", "Full / Refl.", 0.975, 1.025);
      // main/ratio pads are 0.69/0.30 of the canvas height (MakeTwoPad) --
      // same scaling roverlay uses to make label sizes read as equal
      hRatio->GetYaxis()->SetLabelSize(0.092);
      hRatio->GetXaxis()->SetLabelSize(0.092);
      hRatio->GetXaxis()->SetRangeUser(xMin, xMax);
      hRatio->Draw("E1");
      RefLine(cv.ratio, xMin, xMax, 1.0);

      cv.c->cd();
      SavePlot(cv.c, outDir, cone, "etasym", {calibKey, ptKey}, cvName);
      pb.Update();

      // hAbs was not drawn, delete manually
      delete hAbs;
      // hFull, hRefl, hRatio, leg, tex drawn on pads: canvas cascade deletes them
      delete cv.c;
    }
  }
}

#endif
