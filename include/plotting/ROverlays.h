#ifndef L2RESIDUALS_PLOTTING_ROVERLAYS_H
#define L2RESIDUALS_PLOTTING_ROVERLAYS_H

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"

#include "plotting/Style.h"
#include "plotting/Utilities.h"
#include "Binning.h"
#include "Naming.h"
#include "ProgressBar.h"

// ============================================================
// Plot type 4: R_data and R_MC overlay with ratio panel
//
// For each (cone, method, alpha slice, pT slice): one two-panel canvas with
//   top panel  : R_data (blue) and R_MC (red) vs |eta|, reference at 1
//   bottom panel : R_data/R_MC vs |eta|, reference at 1
// ============================================================

inline void PlotROverlay(TFile *fIn, const TString &outDir, const TString &cone,
                         const BinningConfig &bins, ProgressBar &pb,
                         bool useJer = false) {
  TDirectory *dRvals = (TDirectory *)fIn->Get(cone + "/Rvals");

  const int nPt = (int)bins.ptavgSlices.size();
  const int nAlpha = (int)bins.alphaSlices.size();

  for (int m = 0; m < kNMethods; m++) {
    for (int ia = 0; ia < nAlpha; ia++) {
      const auto &aSl = bins.alphaSlices[ia];
      for (int ip = 0; ip < nPt; ip++) {
        const auto &ptSl = bins.ptavgSlices[ip];

        TString ptKey = L2Name::PtKey(ptSl);
        TString alphaKey = L2Name::AlphaKey(aSl);
        TString rdName = L2Name::ObjectName(
            cone, CalibKind("R_data", useJer),
            {L2Name::EtaModeKey(false), ptKey, alphaKey}, {kMethodKeys[m]});
        TString rmName = L2Name::ObjectName(
            cone, CalibKind("R_mc", useJer),
            {L2Name::EtaModeKey(false), ptKey, alphaKey}, {kMethodKeys[m]});

        TH1D *hRd = GetHAny(dRvals, {rdName});
        TH1D *hRm = GetHAny(dRvals, {rmName});

        if (!hRd || !hRm) {
          if (hRd)
            delete hRd;
          if (hRm)
            delete hRm;
          continue;
        }

        TH1D *hRdc = (TH1D *)hRd->Clone(rdName + "_c");
        hRdc->SetDirectory(0);
        TH1D *hRmc = (TH1D *)hRm->Clone(rmName + "_c");
        hRmc->SetDirectory(0);
        TH1D *hRat = RatioH(hRmc, hRdc, rdName + "_rat");

        const TString calibKey = useJer ? "jer" : "jec";
        const TString cvName =
            Form("roverlay_%s_%s_%s_%s_%s", cone.Data(), calibKey.Data(),
                 kMethodKeys[m], ptKey.Data(), alphaKey.Data());
        TwoPad cv = MakeTwoPad(cvName);

        // ---- main pad ----
        cv.main->cd();
        cv.main->SetGridx();
        cv.main->SetGridy();

        StyleH(hRdc, kBlue + 1, 20, 1.5f);
        StyleH(hRmc, kRed + 1, 21, 1.5f);

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
        RefLine(cv.main, (double)kAbsEtaEdges.front(),
                (double)kAbsEtaEdges.back(), 1.0);

        TLegend *leg = new TLegend(0.16, 0.14, 0.50, 0.28);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.048);
        leg->AddEntry(hRmc, "R_{MC}", "lp");
        leg->AddEntry(hRdc, "R_{data}", "lp");
        leg->Draw();

        TLatex *tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.050);
        tex->SetTextFont(62);
        tex->DrawLatex(0.16, 0.91,
                       Form("%s  |  %s  |  %s  |  %s  |  %s", cone.Data(),
                            kMethodLabels[m], ptSl.title.Data(),
                            aSl.title.Data(), CalibTag(useJer).Data()));

        // ---- ratio pad ----
        cv.ratio->cd();
        cv.ratio->SetGridx();
        cv.ratio->SetGridy();

        StyleH(hRat, kBlack, 20, 1.5f);
        TuneRatio(hRat, "|#eta|", "MC/Data", 0.93, 1.07);
        hRat->Draw("E1");
        RefLine(cv.ratio, (double)kAbsEtaEdges.front(),
                (double)kAbsEtaEdges.back(), 1.0);

        cv.c->cd();
        SavePlot(cv.c, outDir, cone, "roverlay", {calibKey, ptKey, alphaKey},
                 cvName);
        pb.Update();

        delete hRd;
        delete hRm;
        delete cv.c; // cascade-deletes hRdc, hRmc, hRat, leg, tex
      }
    }
  }
}

#endif
