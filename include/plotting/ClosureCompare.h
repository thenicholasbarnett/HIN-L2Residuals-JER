#ifndef L2RESIDUALS_PLOTTING_CLOSURE_COMPARE_H
#define L2RESIDUALS_PLOTTING_CLOSURE_COMPARE_H

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"

#include "plotting/FinalCorrections.h" // kFinalsYTitle/kFinalsYTitleJer
#include "plotting/Style.h"
#include "plotting/Utilities.h"
#include "Binning.h"
#include "Naming.h"
#include "ProgressBar.h"

#include <vector>

// Overlay a finals-style curve (JEC L2Residual or JER SF) from its original
// derivation against the same quantity from a closure-check re-derivation,
// with a ratio panel below -- "did the closure land back on the original."
// Reads two separate input files (fIn = original, fInClosure = closure
// check), the first plot family in this repo needing two.
inline void PlotClosureCompare(TFile *fIn, TFile *fInClosure,
                               const TString &outDir, const TString &cone,
                               const BinningConfig &bins, ProgressBar &pb,
                               bool useJer = false) {
  const TString calibKey = useJer ? "jer" : "jec";
  const int nPt = (int)bins.ptavgSlices.size();

  for (int m = 0; m < kNMethods; m++) {
    for (int ieta = 0; ieta < 2; ieta++) { // 0 = |eta|, 1 = full eta
      const bool fullEta = (ieta == 1);
      const TString etaMode = L2Name::EtaModeKey(fullEta);
      const TString xTitle = fullEta ? "#eta_{reco}" : "|#eta_{reco}|";

      for (int ip = 0; ip < nPt; ip++) {
        if (!pb.ShouldKeep())
          continue;
        const auto &ptSl = bins.ptavgSlices[ip];

        TString name = L2Name::ObjectName(cone, CalibKind("intercept", useJer),
                                          {etaMode, L2Name::PtKey(ptSl)},
                                          {kMethodKeys[m]});
        TH1D *hOrig = GetHAny(fIn, {cone + "/" + name});
        TH1D *hClos = GetHAny(fInClosure, {cone + "/" + name});

        // Step 3 grid fallback, independently per file -- either side may
        // be Step 2 (direct histogram) or Step 3 (grid-projected) output.
        if (!hOrig) {
          TString gridName = L2Name::ObjectName(
              cone, CalibKind("corrfinal", useJer), {etaMode}, {kMethodKeys[m]});
          TH2D *h2 = GetH2Any(fIn, {cone + "/" + gridName + "_norm",
                                    gridName + "_norm", cone + "/" + gridName,
                                    gridName});
          if (h2) {
            TDirectory::TContext nodir(nullptr);
            hOrig = h2->ProjectionX(Form("%s_orig_px%d", gridName.Data(), ip),
                                    ip + 1, ip + 1);
            hOrig->SetDirectory(0);
            delete h2;
          }
        }
        if (!hClos) {
          TString gridName = L2Name::ObjectName(
              cone, CalibKind("corrfinal", useJer), {etaMode}, {kMethodKeys[m]});
          TH2D *h2 = GetH2Any(fInClosure, {cone + "/" + gridName + "_norm",
                                           gridName + "_norm",
                                           cone + "/" + gridName, gridName});
          if (h2) {
            TDirectory::TContext nodir(nullptr);
            hClos = h2->ProjectionX(Form("%s_clos_px%d", gridName.Data(), ip),
                                    ip + 1, ip + 1);
            hClos->SetDirectory(0);
            delete h2;
          }
        }

        if (!hOrig || !hClos) {
          delete hOrig;
          delete hClos;
          continue;
        }

        TH1D *hRat = RatioH(hClos, hOrig, name + "_closratio");

        const TString cvName =
            Form("closurecompare_%s_%s_%s_%s_%s", cone.Data(), calibKey.Data(),
                 kMethodKeys[m], etaMode.Data(), L2Name::PtKey(ptSl).Data());
        TwoPad cv = MakeTwoPad(cvName);

        auto [xMin, xMax] = OccupiedRangeWithMargin({hOrig, hClos});
        if (xMin >= xMax) {
          xMin = fullEta ? kEtaEdges.front() : (double)kAbsEtaEdges.front();
          xMax = fullEta ? kEtaEdges.back() : (double)kAbsEtaEdges.back();
        }

        cv.main->cd();
        cv.main->SetLeftMargin(0.16);
        cv.main->SetGridx();
        cv.main->SetGridy();

        StyleH(hOrig, kBlue, 20, 1.5f);
        StyleH(hClos, kRed, 21, 1.5f);

        auto [ylo, yhi] = YRange({hOrig, hClos});
        hOrig->GetXaxis()->SetRangeUser(xMin, xMax);
        hOrig->GetYaxis()->SetRangeUser(ylo, yhi);
        hOrig->GetYaxis()->SetTitle(useJer ? kFinalsYTitleJer : kFinalsYTitle);
        hOrig->GetYaxis()->SetTitleSize(0.040);
        hOrig->GetYaxis()->SetTitleOffset(1.5);
        hOrig->GetYaxis()->SetLabelSize(0.055);
        hOrig->GetYaxis()->CenterTitle();
        hOrig->GetXaxis()->SetLabelSize(0.0);
        hOrig->GetXaxis()->SetTitle("");
        hOrig->SetTitle("");

        hOrig->Draw("E1");
        hClos->Draw("E1 same");
        RefLine(cv.main, xMin, xMax, 1.0);

        TLegend *leg = new TLegend(0.18, 0.14, 0.45, 0.28);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.048);
        leg->AddEntry(hOrig, "Original", "lp");
        leg->AddEntry(hClos, "Closure", "lp");
        leg->Draw();

        DrawInfoLegend(0.55, 0.60, 0.94, 0.90,
                       {cone, CalibMethodLabel(m, useJer), ptSl.title});

        cv.ratio->cd();
        cv.ratio->SetLeftMargin(0.16);
        cv.ratio->SetGridx();
        cv.ratio->SetGridy();

        StyleH(hRat, kBlack, 20, 1.5f);
        hRat->GetXaxis()->SetRangeUser(xMin, xMax);
        TuneRatio(hRat, xTitle, "Closure / Original", 0.93, 1.07);
        hRat->Draw("E1");
        RefLine(cv.ratio, xMin, xMax, 1.0);

        cv.c->cd();
        SavePlot(cv.c, outDir, cone, "closurecompare",
                {calibKey, etaMode, L2Name::PtKey(ptSl)}, cvName);
        pb.Update();

        delete cv.c; // cascade-deletes hOrig, hClos, hRat, leg (and the
                     // DrawInfoLegend/RefLine objects), same convention
                     // as PlotROverlay -- hOrig/hClos are drawn directly,
                     // not cloned, so no separate delete needed
      }
    }
  }
}

#endif
