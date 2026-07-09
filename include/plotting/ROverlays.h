#ifndef L2RESIDUALS_PLOTTING_ROVERLAYS_H
#define L2RESIDUALS_PLOTTING_ROVERLAYS_H

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"

#include "plotting/Style.h"
#include "plotting/Utilities.h"
#include "Binning.h"
#include "Naming.h"
#include "ProgressBar.h"

// R_data/R_MC overlay with ratio panel. Two-panel canvas per (cone, method,
// alpha slice, pT slice): top = R_data (blue) and R_MC (red) vs |eta|,
// bottom = R_MC/R_data (JEC) or R_data/R_MC (JER, matching the smearing-SF
// convention), both with a reference line at 1.

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
        if (!pb.ShouldKeep())
          continue;
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
        // JEC: R_MC/R_data (matches the correction convention). JER: inverted,
        // R_data/R_MC (matches the smearing-SF convention) -- see
        // CalibrationExtractor.cxx's ratio computation.
        TH1D *hRat = useJer ? RatioH(hRdc, hRmc, rdName + "_rat")
                            : RatioH(hRmc, hRdc, rdName + "_rat");

        const TString calibKey = useJer ? "jer" : "jec";
        const TString cvName =
            Form("roverlay_%s_%s_%s_%s_%s", cone.Data(), calibKey.Data(),
                 kMethodKeys[m], ptKey.Data(), alphaKey.Data());
        TwoPad cv = MakeTwoPad(cvName);

        // end the axis one bin past the last populated one, same as finals
        // (they share this |eta| binning)
        auto [xMin, xMax] = OccupiedRangeWithMargin({hRdc, hRmc});
        if (xMin >= xMax) {
          xMin = (double)kAbsEtaEdges.front();
          xMax = (double)kAbsEtaEdges.back();
        }

        // main pad
        cv.main->cd();
        cv.main->SetLeftMargin(0.16); // wider than TwoPad's default 0.13 --
                                      // the explicit fraction title needs it
        cv.main->SetGridx();
        cv.main->SetGridy();

        // two series only: always kBlue/kRed
        StyleH(hRdc, kBlue, 20, 1.5f);
        StyleH(hRmc, kRed, 21, 1.5f);

        auto [ylo, yhi] = YRange({hRdc, hRmc});
        hRmc->GetXaxis()->SetRangeUser(xMin, xMax);
        hRmc->GetYaxis()->SetRangeUser(ylo, yhi);
        hRmc->GetYaxis()->SetTitle(kAsymFractionTitle);
        hRmc->GetYaxis()->SetTitleSize(0.040);
        hRmc->GetYaxis()->SetTitleOffset(1.35);
        hRmc->GetYaxis()->SetLabelSize(0.055);
        hRmc->GetYaxis()->CenterTitle();
        hRmc->GetXaxis()->SetLabelSize(0.0);
        hRmc->GetXaxis()->SetTitle("");
        hRmc->SetTitle("");

        hRmc->Draw("E1");
        hRdc->Draw("E1 same");
        RefLine(cv.main, xMin, xMax, 1.0);

        // close to the left border/barrel, where R sits flat near 1 and
        // stays out of the markers' way
        TLegend *leg = new TLegend(0.18, 0.14, 0.40, 0.28);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.048);
        leg->AddEntry(hRmc, "R_{MC}", "lp");
        leg->AddEntry(hRdc, "R_{data}", "lp");
        leg->Draw();

        DrawInfoLegend(
            0.58, 0.60, 0.94, 0.90,
            {cone, CalibMethodLabel(m, useJer), ptSl.title, aSl.title});

        // ratio pad
        cv.ratio->cd();
        cv.ratio->SetLeftMargin(0.16); // match the main pad so the x-axes align
        cv.ratio->SetGridx();
        cv.ratio->SetGridy();

        StyleH(hRat, kBlack, 20, 1.5f);
        hRat->GetXaxis()->SetRangeUser(xMin, xMax);
        TuneRatio(hRat, "|#eta_{reco}|", useJer ? "Data/MC" : "MC/Data", 0.93,
                 1.07);
        hRat->Draw("E1");
        RefLine(cv.ratio, xMin, xMax, 1.0);

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
