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

      StyleH(hFull, kColFull, 20, 2.f);
      StyleH(hRefl, kColRefl, 21, 2.f);
      hRefl->SetLineStyle(2);

      auto [ylo, yhi] = YRange({hFull, hRefl});
      hFull->GetYaxis()->SetRangeUser(ylo, yhi);
      hFull->GetYaxis()->SetTitle(CalibYTitle(useJer));
      hFull->GetYaxis()->SetTitleSize(0.055);
      hFull->GetYaxis()->SetTitleOffset(1.10);
      hFull->GetYaxis()->SetLabelSize(0.050);
      hFull->GetXaxis()->SetLabelSize(0.0);
      hFull->GetXaxis()->SetTitle("");
      hFull->SetTitle("");

      hFull->Draw("E1");
      hRefl->Draw("E1 same");
      RefLine(cv.main, kEtaEdges.front(), kEtaEdges.back(), 1.0);

      TLegend *leg = new TLegend(0.16, 0.14, 0.56, 0.28);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.048);
      leg->AddEntry(hFull, "Full #eta", "lp");
      leg->AddEntry(hRefl, "|#eta| reflected", "lp");
      leg->Draw();

      TLatex *tex = new TLatex();
      tex->SetNDC();
      tex->SetTextSize(0.051);
      tex->SetTextFont(62);
      tex->DrawLatex(0.16, 0.91,
                     Form("%s   |   %s   |   %s   |   %s", cone.Data(),
                          kMethodLabels[m], sl.title.Data(),
                          CalibTag(useJer).Data()));

      // ratio pad
      cv.ratio->cd();
      cv.ratio->SetGridx();
      cv.ratio->SetGridy();

      StyleH(hRatio, kBlack, 20, 1.5f);
      TuneRatio(hRatio, "#eta", "Full / Refl.", 0.975, 1.025);
      hRatio->Draw("E1");
      RefLine(cv.ratio, kEtaEdges.front(), kEtaEdges.back(), 1.0);

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
