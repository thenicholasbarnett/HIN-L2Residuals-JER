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
#include "AnalysisConfig.h"
#include "Binning.h"
#include "Colors.h"
#include "Naming.h"
#include "ProgressBar.h"

#include <vector>

static const TString kFinalsYTitle =
    "k_{FSR} #scale[0.35]{#bullet} (#frac{R_{MC}}{R_{Data}})_{#alpha=0.30}";
static const TString kFinalsYTitleJer = "JER SF_{#alpha=0.30}";

// Final extrapolated values, all pT slices overlaid.
// finals_{cone}_{method}_abseta/fulleta: kFSR*R_MC/R_data at alpha=0.30 vs eta_probe.

inline void PlotFinals(TFile *fIn, const TString &outDir, const TString &cone,
                       const BinningConfig &bins, ProgressBar &pb,
                       bool isClosure = false, bool useJer = false) {
  const TString calibKey = useJer ? "jer" : "jec";
  for (int m = 0; m < kNMethods; m++) {
    for (int ieta = 0; ieta < 2; ieta++) { // 0 = |eta|, 1 = full eta
      if (!pb.ShouldKeep())
        continue;
      const bool fullEta = (ieta == 1);
      const TString xTitle = fullEta ? "#eta^{probe}" : "|#eta^{probe}|";
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
      // projecting each pT slice's row.
      if (!anyValid) {
        // try norm first (matches PlotAlphaFit), then direct
        TString gridName = L2Name::ObjectName(
            cone, CalibKind("corrfinal", useJer), {etaMode}, {kMethodKeys[m]});
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
      c->SetLeftMargin(0.14); // matches finalscone's tightened margin --
                              // same tall fraction title, same fix
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

      // one legend (cone/method info + pT-slice entries) -- upper-middle
      // for full eta, upper-left for |eta| (its data only occupies x>=0,
      // so upper-middle would sit on top of real points)
      const double legX1 = fullEta ? 0.345 : 0.15;
      const double legX2 = fullEta ? 0.675 : 0.42;
      const double legY1 = fullEta ? 0.62 : 0.68;
      const double legY2 = legY1 + 0.026 * (double)(bins.ptavgSlices.size() + 2);
      TLegend *leg = new TLegend(legX1, legY1, legX2, legY2);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42); // non-bold -- was inheriting tdrStyle's bold default
      leg->SetTextSize(0.023);
      leg->SetTextAlign(22); // center entries instead of hugging the marker column
      leg->AddEntry((TObject *)nullptr, cone.Data(), "");
      leg->AddEntry((TObject *)nullptr, CalibMethodLabel(m, useJer).Data(), "");

      bool first = true;
      for (int ip = 0; ip < (int)bins.ptavgSlices.size(); ip++) {
        if (!hists[ip])
          continue;
        const auto &ptSl = bins.ptavgSlices[ip];
        StyleH(hists[ip], ptSl.color, 20, 1.5f); // circles for every pT
                                                 // slice -- kMethodStyles
                                                 // varies by method, but
                                                 // this plot fixes the
                                                 // method and varies pT
        hists[ip]->GetXaxis()->SetRangeUser(xMin, xMax);
        hists[ip]->GetYaxis()->SetRangeUser(ylo, yhi);
        // closure mode reads the fixed 0.95-1.05 band as "does this equal
        // 1", not "what's the correction shape" -- a plain closure-factor
        // label fits that better than the derivation-shaped kFSR fraction
        hists[ip]->GetYaxis()->SetTitle(
            isClosure ? (useJer ? TString("JER SF Closure Factor")
                                : TString("L2Residual Closure Factors"))
                     : (useJer ? kFinalsYTitleJer : kFinalsYTitle));
        hists[ip]->GetYaxis()->SetTitleSize(isClosure ? 0.040 : 0.030);
        hists[ip]->GetYaxis()->SetTitleOffset(isClosure ? 1.38 : 1.50);
        hists[ip]->GetYaxis()->SetLabelSize(0.032);
        hists[ip]->GetXaxis()->SetTitle(xTitle);
        hists[ip]->GetXaxis()->SetTitleSize(0.052);
        hists[ip]->GetXaxis()->SetLabelSize(0.032);
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

      leg->Draw();
      DrawAsymHeader(0.14);

      SavePlot(c, outDir, cone, "finals", {calibKey, etaMode}, cvName);
      pb.Update();

      delete c; // cascade-deletes hists, leg, rl (and rl99/rl101 if drawn)
    }
  }
}

// Same Step 3 corrfinal/corrfinal_jer grid PlotFinals slices into 1D pT-slice
// overlays above, drawn directly as a 2D eta-vs-pT_avg color map instead.
// Step 3 output only -- no Step 2 fallback, since Step 2 never stores this
// shape (per-slice histograms only).
inline void PlotFinalsGrid(TFile *fIn, const TString &outDir,
                           const TString &cone, ProgressBar &pb,
                           bool useJer = false) {
  const TString calibKey = useJer ? "jer" : "jec";
  for (int m = 0; m < kNMethods; m++) {
    for (int ieta = 0; ieta < 2; ieta++) { // 0 = |eta|, 1 = full eta
      if (!pb.ShouldKeep())
        continue;
      const bool fullEta = (ieta == 1);
      const TString xTitle = fullEta ? "#eta_{probe}" : "|#eta_{probe}|";
      const TString etaMode = L2Name::EtaModeKey(fullEta);

      TString gridName = L2Name::ObjectName(
          cone, CalibKind("corrfinal", useJer), {etaMode}, {kMethodKeys[m]});
      TH2D *h2 =
          GetH2Any(fIn, {cone + "/" + gridName + "_norm", gridName + "_norm",
                        cone + "/" + gridName, gridName});
      if (!h2) {
        continue;
      }

      const TString cvName =
          Form("finalsgrid_%s_%s_%s_%s", cone.Data(), calibKey.Data(),
               kMethodKeys[m], etaMode.Data());
      TCanvas *c = new TCanvas(cvName, "", 800, 600);
      RealAspectRatio(c);
      c->SetLeftMargin(0.13);
      c->SetRightMargin(0.17);
      c->SetBottomMargin(0.12);

      auto [zlo, zhi] = ZRange(h2);
      h2->SetTitle("");
      h2->GetXaxis()->SetTitle(xTitle);
      h2->GetXaxis()->CenterTitle();
      h2->GetYaxis()->CenterTitle();
      h2->GetZaxis()->SetTitle(useJer ? "JER SF" : "Correction factor");
      h2->GetZaxis()->SetRangeUser(zlo, zhi);
      h2->Draw("COLZ");

      DrawAsymHeader(0.13);

      SavePlot(c, outDir, cone, "finals", {calibKey, etaMode}, cvName);
      pb.Update();

      delete c; // cascade-deletes h2 (drawn on the canvas, same as PlotFinals)
    }
  }
}

// Cone-size overlay: for each (method, etaMode, pT slice), overlay all cone
// sizes across eta -- the transpose of PlotFinals' pT-slice overlay for one
// fixed cone. Reads across every cone's directory in fIn, so unlike
// PlotFinals/PlotFinalsGrid (called once per cone from runPlotting.C's
// per-cone loop) this is called exactly once.
inline void PlotFinalsByCone(TFile *fIn, const TString &outDir,
                             const BinningConfig &bins, ProgressBar &pb,
                             bool isClosure = false, bool useJer = false) {
  static Color_t (*const kConeColors[])() = {KlimtPink, KlimtRed, KlimtYellow,
                                             KlimtGreen, KlimtBlue};
  static constexpr int kNConeColors =
      sizeof(kConeColors) / sizeof(kConeColors[0]);

  const std::vector<TString> &cones = Config().coneLabels;
  const int nCones = (int)cones.size();
  const TString calibKey = useJer ? "jer" : "jec";

  for (int m = 0; m < kNMethods; m++) {
    for (int ieta = 0; ieta < 2; ieta++) { // 0 = |eta|, 1 = full eta
      const bool fullEta = (ieta == 1);
      const TString xTitle = fullEta ? "#eta^{probe}" : "|#eta^{probe}|";
      const double xFullMin =
          fullEta ? kEtaEdges.front() : (double)kAbsEtaEdges.front();
      const double xFullMax =
          fullEta ? kEtaEdges.back() : (double)kAbsEtaEdges.back();
      const TString etaMode = L2Name::EtaModeKey(fullEta);

      for (int ip = 0; ip < (int)bins.ptavgSlices.size(); ip++) {
        if (!pb.ShouldKeep())
          continue;
        const auto &ptSl = bins.ptavgSlices[ip];

        std::vector<TH1D *> hists(nCones, nullptr);
        for (int ic = 0; ic < nCones; ic++) {
          TString name = L2Name::ObjectName(
              cones[ic], CalibKind("intercept", useJer),
              {etaMode, L2Name::PtKey(ptSl)}, {kMethodKeys[m]});
          hists[ic] = GetHAny(fIn, {cones[ic] + "/" + name});
        }

        bool anyValid = false;
        for (auto *h : hists)
          if (h) {
            anyValid = true;
            break;
          }

        // Step 3 grid fallback per cone, same convention as PlotFinals --
        // project this pT slice's row out of each cone's own grid.
        if (!anyValid) {
          for (int ic = 0; ic < nCones; ic++) {
            TString gridName =
                L2Name::ObjectName(cones[ic], CalibKind("corrfinal", useJer),
                                   {etaMode}, {kMethodKeys[m]});
            TString gridNormName = gridName + "_norm";
            TH2D *h2 = GetH2Any(
                fIn, {cones[ic] + "/" + gridNormName, gridNormName,
                     cones[ic] + "/" + gridName, gridName});
            if (!h2) {
              continue;
            }
            TH1D *px;
            {
              TDirectory::TContext nodir(nullptr);
              px = h2->ProjectionX(Form("%s_px%d", gridName.Data(), ip),
                                   ip + 1, ip + 1);
            }
            px->SetDirectory(0);
            hists[ic] = px;
            anyValid = true;
            delete h2;
          }
        }

        if (!anyValid) {
          for (auto *h : hists)
            delete h;
          continue;
        }

        // find ak4PF for the ratio panel -- absent from this run period's
        // cone list, the ratio panel is simply skipped below
        int refIdx = -1;
        for (int ic = 0; ic < nCones; ic++) {
          if (cones[ic] == "ak4PF") {
            refIdx = ic;
            break;
          }
        }

        const TString cvName =
            Form("finalscone_%s_%s_%s_%s", calibKey.Data(), kMethodKeys[m],
                 etaMode.Data(), L2Name::PtKey(ptSl).Data());
        TwoPad cv = MakeTwoPad(cvName);

        // main pad
        cv.main->cd();
        cv.main->SetLeftMargin(0.14); // wider than TwoPad's default -- the
                                      // explicit fraction title needs it
        cv.main->SetGridx();
        cv.main->SetGridy();

        double ylo = 0.9, yhi = 1.6;
        if (isClosure) {
          ylo = 0.95;
          yhi = 1.05;
        }

        auto [xMin, xMax] = OccupiedRangeWithMargin(hists, ylo, yhi);
        if (xMin >= xMax) {
          xMin = xFullMin;
          xMax = xFullMax;
        }

        // non-closure yhi was a flat 1.6 regardless of what's actually
        // plotted -- rescale to the real max (plus a small pad) instead,
        // same padding convention as the alpha combined-overlay's dynamic
        // range. Closure mode keeps its fixed 0.95-1.05 window untouched.
        if (!isClosure) {
          double realMax = -1e9;
          bool any = false;
          for (auto *h : hists) {
            if (!h)
              continue;
            for (int i = 1; i <= h->GetNbinsX(); i++) {
              const double v = h->GetBinContent(i), e = h->GetBinError(i);
              if (v == 0 && e == 0)
                continue;
              const double x = h->GetBinCenter(i);
              if (x < xMin || x > xMax)
                continue;
              realMax = std::max(realMax, v + e);
              any = true;
            }
          }
          if (any)
            yhi = realMax + 0.10 * (realMax - ylo);
        }

        // one legend, upper-middle of the plot (clear of the curve's
        // central dip) instead of two separate boxes or dead center --
        // centered on the frame itself (same margins in both eta modes)
        // rather than data range, since the frame's NDC width is fixed
        const double legX1 = 0.405, legX2 = 0.735;
        const double legY1 = 0.62;
        const double legY2 = legY1 + 0.032 * (double)(nCones + 2);
        TLegend *leg = new TLegend(legX1, legY1, legX2, legY2);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.028);
        leg->SetTextAlign(22); // center entries instead of hugging the marker column
        leg->AddEntry((TObject *)nullptr, ptSl.title.Data(), "");
        leg->AddEntry((TObject *)nullptr, CalibMethodLabel(m, useJer).Data(),
                     "");

        bool first = true;
        for (int ic = 0; ic < nCones; ic++) {
          if (!hists[ic])
            continue;
          StyleH(hists[ic], kConeColors[ic % kNConeColors](), kMethodStyles[m],
                1.5f);
          hists[ic]->GetXaxis()->SetRangeUser(xMin, xMax);
          hists[ic]->GetYaxis()->SetRangeUser(ylo, yhi);
          hists[ic]->GetYaxis()->SetTitle(useJer ? kFinalsYTitleJer
                                                 : kFinalsYTitle);
          hists[ic]->GetYaxis()->SetTitleSize(0.048);
          hists[ic]->GetYaxis()->SetTitleOffset(1.20);
          hists[ic]->GetYaxis()->SetLabelSize(0.038);
          hists[ic]->GetXaxis()->SetLabelSize(0.0);
          hists[ic]->GetXaxis()->SetTitle("");
          hists[ic]->GetYaxis()->CenterTitle();
          hists[ic]->SetTitle("");
          hists[ic]->Draw(first ? "E1" : "E1 same");
          first = false;
          leg->AddEntry(hists[ic], cones[ic], "lp");
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

        leg->Draw();
        DrawCMSInternalHeader(0.14, 0.96, 0.915, 0.040); // matches etasym/methods

        // ratio pad: each cone's correction relative to ak4PF's, same idea
        // as the methods plot's "/ Gauss" panel
        cv.ratio->cd();
        cv.ratio->SetLeftMargin(0.14); // match the main pad so the x-axes align
        cv.ratio->SetGridx();
        cv.ratio->SetGridy();

        std::vector<TH1D *> ratios;
        std::vector<int> ratioConeIdx;
        const double rlo = 0.975, rhi = 1.025;
        if (refIdx >= 0 && hists[refIdx]) {
          for (int ic = 0; ic < nCones; ic++) {
            if (ic == refIdx || !hists[ic])
              continue;
            TH1D *r = RatioH(hists[ic], hists[refIdx],
                             Form("%s_finalscone_rat%d", cones[ic].Data(), ic));
            r->GetXaxis()->SetRangeUser(xMin, xMax);
            ratios.push_back(r);
            ratioConeIdx.push_back(ic);
          }
        }

        bool firstR = true;
        for (size_t k = 0; k < ratios.size(); k++) {
          const int ic = ratioConeIdx[k];
          StyleH(ratios[k], kConeColors[ic % kNConeColors](), kMethodStyles[m],
                1.5f);
          TuneRatio(ratios[k], xTitle, "/ ak4PF", rlo, rhi);
          ratios[k]->GetYaxis()->SetLabelSize(0.092);
          ratios[k]->GetXaxis()->SetLabelSize(0.092);
          ratios[k]->Draw(firstR ? "E1" : "E1 same");
          firstR = false;
        }
        if (!firstR) {
          RefLine(cv.ratio, xMin, xMax, 1.0);
        }

        cv.c->cd();
        SavePlot(cv.c, outDir, "finalscone", "finals",
                {calibKey, etaMode, L2Name::PtKey(ptSl)}, cvName);
        pb.Update();

        delete cv.c; // cascade-deletes hists, ratios, leg, rl (and rl99/rl101 if drawn)
      }
    }
  }
}

#endif
