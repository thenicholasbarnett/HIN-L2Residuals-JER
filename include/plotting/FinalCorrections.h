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
    "k_{FSR} #scale[0.35]{#bullet} #frac{R_{MC}}{R_{Data}} |_{#alpha=0.30}";

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
      c->SetLeftMargin(0.16);
      c->SetGridx();
      c->SetGridy();

      // end the axis one bin past the last populated one instead of running
      // out to the binning's full extent, which is often much wider than
      // where jets actually get reconstructed
      auto [xMin, xMax] = OccupiedRangeWithMargin(hists);
      if (xMin >= xMax) {
        xMin = xFullMin;
        xMax = xFullMax;
      }

      double ylo = 0.9, yhi = 1.6;
      // closure passes check R_MC/R_data ~= 1 to a tight tolerance
      if (isClosure) {
        ylo = 0.95;
        yhi = 1.05;
      }

      // one legend (cone/method + all pT slices) instead of two separate
      // boxes -- upper-left for |eta|, upper-middle for full eta so it
      // doesn't sit over the eta<0 half of the data
      const double legX1 = fullEta ? 0.34 : 0.19;
      const double legX2 = fullEta ? 0.67 : 0.52;
      TLegend *leg = new TLegend(legX1, 0.50, legX2, 0.88);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42); // non-bold -- was inheriting tdrStyle's bold default
      leg->SetTextSize(0.028);
      leg->AddEntry((TObject *)nullptr, cone.Data(), "");
      leg->AddEntry((TObject *)nullptr, CalibMethodLabel(m, useJer).Data(), "");

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
        hists[ip]->GetYaxis()->SetTitleOffset(1.35);
        hists[ip]->GetYaxis()->SetLabelSize(0.048);
        hists[ip]->GetXaxis()->SetTitle(xTitle);
        hists[ip]->GetXaxis()->SetTitleSize(0.052);
        hists[ip]->GetXaxis()->SetLabelSize(0.048);
        hists[ip]->GetXaxis()->CenterTitle();
        hists[ip]->GetYaxis()->CenterTitle();
        hists[ip]->SetTitle("");
        hists[ip]->Draw(first ? "E1" : "E1 same");
        first = false;
        TString label = ptSl.title;
        label.ReplaceAll(" GeV", "");
        leg->AddEntry(hists[ip], label, "lp");
      }

      TLine *rl = new TLine(xMin, 1.0, xMax, 1.0);
      rl->SetLineStyle(2);
      rl->SetLineColor(kGray + 2);
      rl->SetLineWidth(1);
      rl->Draw();

      if (isClosure) {
        // distinct accent color -- plain gray would be invisible against the grid
        TLine *rl99 = new TLine(xMin, 0.99, xMax, 0.99);
        rl99->SetLineStyle(3);
        rl99->SetLineColor(HiroshigeLightRed());
        rl99->SetLineWidth(2);
        rl99->Draw();

        TLine *rl101 = new TLine(xMin, 1.01, xMax, 1.01);
        rl101->SetLineStyle(3);
        rl101->SetLineColor(HiroshigeLightRed());
        rl101->SetLineWidth(2);
        rl101->Draw();
      }

      leg->Draw();
      DrawAsymHeader(0.16);

      SavePlot(c, outDir, cone, "finals", {calibKey, etaMode}, cvName);
      pb.Update();

      delete c; // cascade-deletes hists, leg, rl (and rl99/rl101 if drawn)
    }
  }
}

#endif
