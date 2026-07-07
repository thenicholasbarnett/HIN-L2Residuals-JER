#ifndef L2RESIDUALS_PLOTTING_FINAL_CORRECTIONS_H
#define L2RESIDUALS_PLOTTING_FINAL_CORRECTIONS_H

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"

#include "plotting/Style.h"
#include "plotting/Utilities.h"
#include "Binning.h"
#include "Naming.h"
#include "ProgressBar.h"

#include <vector>

static const TString kFinalsYTitle =
    "k_{FSR} #scale[0.35]{#bullet} (#frac{R_{MC}}{R_{Data}})_{#alpha=0.30}";

// Final extrapolated values, all pT slices overlaid.
// finals_{cone}_{method}_abseta/fulleta: kFSR*R_MC/R_data at alpha=0.30 vs eta_reco.

inline void PlotFinals(TFile *fIn, const TString &outDir, const TString &cone,
                       const BinningConfig &bins, ProgressBar &pb,
                       bool isClosure = false, bool useJer = false) {
  const TString calibKey = useJer ? "jer" : "jec";
  for (int m = 0; m < kNMethods; m++) {
    for (int ieta = 0; ieta < 2; ieta++) { // 0 = |eta|, 1 = full eta
      if (!pb.ShouldKeep())
        continue;
      const bool fullEta = (ieta == 1);
      const TString xTitle = fullEta ? "#eta_{reco}" : "|#eta_{reco}|";
      const double xFullMin =
          fullEta ? kEtaEdges.front() : (double)kAbsEtaEdges.front();
      const double xFullMax =
          fullEta ? kEtaEdges.back() : (double)kAbsEtaEdges.back();
      const TString etaMode = L2Name::EtaModeKey(fullEta);

      std::vector<TH1D *> hists;
      for (const auto &ptSl : bins.ptavgSlices) {
        TString name = L2Name::ObjectName(cone, CalibKind("intercept", useJer),
                                          {etaMode, L2Name::PtKey(ptSl)},
                                          {kMethodKeys[m]});
        hists.push_back(GetHAny(fIn, {cone + "/" + name}));
      }

      bool anyValid = false;
      for (auto *h : hists)
        if (h) {
          anyValid = true;
          break;
        }

      // Step 3 stores one TH2D grid per (cone, etaMode, method); fall back to
      // projecting each pT slice's row. JEC-only -- no JER SF grid exists yet.
      if (!anyValid && !useJer) {
        // try norm first (matches PlotAlphaFit), then direct
        TString gridName =
            L2Name::ObjectName(cone, "corrfinal", {etaMode}, {kMethodKeys[m]});
        TString gridNormName = gridName + "_norm";
        TH2D *h2 = GetH2Any(fIn, {cone + "/" + gridNormName, gridNormName,
                                  cone + "/" + gridName, gridName});
        if (h2) {
          for (int ip = 0; ip < (int)bins.ptavgSlices.size(); ip++) {
            TH1D *px;
            {
              TDirectory::TContext nodir(nullptr);
              px = h2->ProjectionX(Form("%s_px%d", gridName.Data(), ip), ip + 1,
                                   ip + 1);
            }
            px->SetDirectory(0);
            hists[ip] = px;
            anyValid = true;
          }
          delete h2;
        }
      }

      if (!anyValid) {
        for (auto *h : hists)
          delete h;
        continue;
      }

      const TString cvName =
          Form("finals_%s_%s_%s_%s", cone.Data(), calibKey.Data(),
               kMethodKeys[m], etaMode.Data());
      TCanvas *c = new TCanvas(cvName, "", 800, 800);
      RealAspectRatio(c);
      c->SetLeftMargin(0.20);
      c->SetBottomMargin(0.12);
      c->SetGridx();
      c->SetGridy();

      double ylo = 0.9, yhi = 1.6;
      // closure passes check R_MC/R_data ~= 1 to a tight tolerance
      if (isClosure) {
        ylo = 0.95;
        yhi = 1.05;
      }

      // end the axis one bin past the last populated one instead of running
      // out to the binning's full extent, which is often much wider than
      // where jets actually get reconstructed -- restricted to points that
      // actually land inside [ylo, yhi], so a closure's tight y-window
      // doesn't drag the x-axis out to cover a point that isn't even drawn
      auto [xMin, xMax] = OccupiedRangeWithMargin(hists, ylo, yhi);
      if (xMin >= xMax) {
        xMin = xFullMin;
        xMax = xFullMax;
      }

      // cone/method header stays top -- upper-left for |eta|, upper-middle
      // for full eta so it doesn't sit over the eta<0 half of the data;
      // pulled in almost flush with the frame, just enough buffer to clear
      // the inward-facing tick marks, not sitting on top of them
      const double legX1 = fullEta ? 0.34 : 0.215;
      const double legX2 = fullEta ? 0.67 : 0.54;
      TLegend *legInfo = new TLegend(legX1, 0.82, legX2, 0.895);
      legInfo->SetBorderSize(0);
      legInfo->SetFillStyle(0);
      legInfo->SetTextFont(42); // non-bold -- was inheriting tdrStyle's bold default
      legInfo->SetTextSize(0.028);
      legInfo->AddEntry((TObject *)nullptr, cone.Data(), "");
      legInfo->AddEntry((TObject *)nullptr, CalibMethodLabel(m, useJer).Data(), "");

      // pT-slice entries sit low in the same column instead of stacked under
      // the header -- for closures the data band hugs 1.0 so the bottom
      // border itself is clear, right down to the tick-mark buffer; the
      // normal wide 0.9-1.6 range has the opposite problem (1.0 sits close
      // to the *bottom* of that range), so it keeps the higher anchor that
      // already sat clear of the data. Box height scales with slice count.
      const double ptLegY1 = isClosure ? 0.135 : 0.29;
      const double ptLegY2 = ptLegY1 + 0.032 * (double)bins.ptavgSlices.size();
      TLegend *legPt = new TLegend(legX1, ptLegY1, legX2, ptLegY2);
      legPt->SetBorderSize(0);
      legPt->SetFillStyle(0);
      legPt->SetTextFont(42);
      legPt->SetTextSize(0.028);

      bool first = true;
      for (int ip = 0; ip < (int)bins.ptavgSlices.size(); ip++) {
        if (!hists[ip])
          continue;
        const auto &ptSl = bins.ptavgSlices[ip];
        StyleH(hists[ip], ptSl.color, kMethodStyles[m], 1.5f);
        hists[ip]->GetXaxis()->SetRangeUser(xMin, xMax);
        hists[ip]->GetYaxis()->SetRangeUser(ylo, yhi);
        hists[ip]->GetYaxis()->SetTitle(kFinalsYTitle);
        hists[ip]->GetYaxis()->SetTitleSize(0.048);
        hists[ip]->GetYaxis()->SetTitleOffset(2.0);
        hists[ip]->GetYaxis()->SetLabelSize(0.038);
        hists[ip]->GetXaxis()->SetTitle(xTitle);
        hists[ip]->GetXaxis()->SetTitleSize(0.052);
        hists[ip]->GetXaxis()->SetLabelSize(0.038);
        hists[ip]->GetXaxis()->CenterTitle();
        hists[ip]->GetYaxis()->CenterTitle();
        hists[ip]->SetTitle("");
        hists[ip]->Draw(first ? "E1" : "E1 same");
        first = false;
        TString label = ptSl.title;
        label.ReplaceAll(" GeV", "");
        legPt->AddEntry(hists[ip], label, "lp");
      }

      TLine *rl = new TLine(xMin, 1.0, xMax, 1.0);
      rl->SetLineStyle(2);
      rl->SetLineColor(kGray + 2);
      rl->SetLineWidth(1);
      rl->Draw();

      if (isClosure) {
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

      legInfo->Draw();
      legPt->Draw();
      DrawAsymHeader(0.20);

      SavePlot(c, outDir, cone, "finals", {calibKey, etaMode}, cvName);
      pb.Update();

      delete c; // cascade-deletes hists, legInfo, legPt, rl (and rl99/rl101 if drawn)
    }
  }
}

#endif
