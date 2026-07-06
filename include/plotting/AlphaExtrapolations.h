#ifndef L2RESIDUALS_PLOTTING_ALPHA_EXTRAPOLATIONS_H
#define L2RESIDUALS_PLOTTING_ALPHA_EXTRAPOLATIONS_H

#include "TFile.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"

#include "plotting/Style.h"
#include "plotting/Utilities.h"
#include "Binning.h"
#include "Naming.h"
#include "ProgressBar.h"

#include <algorithm>

static const TString kAlphaYTitle =
    "R_{MC}/R_{data}|_{#alpha} / R_{MC}/R_{data}|_{#alpha=0.30}";

// Alpha fit plots. One canvas per (cone, method, pT slice, eta bin): all
// alpha threshold points, fit line drawn only through [0, 0.30]. A second
// canvas per (cone, pT slice, eta bin), in the same directory, overlays all
// 4 methods together.

inline void PlotAlphaFit(TFile *fIn, const TString &outDir, const TString &cone,
                         const BinningConfig &bins, ProgressBar &pb,
                         bool useJer = false) {
  TDirectory *dGraphs = (TDirectory *)fIn->Get(cone + "/graphs");

  const int nPt = (int)bins.ptavgSlices.size();
  const int nEta = (int)kAbsEtaEdges.size() - 1;
  const TString calibKey = useJer ? "jer" : "jec";

  for (int m = 0; m < kNMethods; m++) {
    for (int ip = 0; ip < nPt; ip++) {
      const auto &ptSl = bins.ptavgSlices[ip];
      for (int ie = 0; ie < nEta; ie++) {
        TString etaKey = L2Name::EtaKey(ie, false);
        TString ptKey = L2Name::PtKey(ptSl);
        TString gname = L2Name::ObjectName(
            cone, CalibKind("R", useJer),
            {L2Name::EtaModeKey(false), etaKey, ptKey}, {kMethodKeys[m]});
        TString gnorm = gname + "_norm";

        TGraphErrors *gr = GetGraphAny(dGraphs, {gnorm, gname});

        if (!gr || gr->GetN() < 2) {
          continue;
        }
        TGraphErrors *gc = (TGraphErrors *)gr->Clone(gnorm + "_c");

        const TString cvName =
            Form("alphafit_%s_%s_%s_%s_%s", cone.Data(), calibKey.Data(),
                 kMethodKeys[m], etaKey.Data(), ptKey.Data());
        TCanvas *c = new TCanvas(cvName, "", 800, 800);
        RealAspectRatio(c);
        c->SetLeftMargin(0.18);
        c->SetGridx();
        c->SetGridy();

        // single series (markers + one fit line): always kRed/kBlue
        gc->SetMarkerStyle(20);
        gc->SetMarkerColor(kRed);
        gc->SetLineColor(kRed);
        gc->SetMarkerSize(0.9);

        gc->GetXaxis()->SetTitle("#alpha threshold");
        gc->GetYaxis()->SetTitle(kAlphaYTitle);
        gc->GetXaxis()->CenterTitle();
        gc->GetYaxis()->CenterTitle();
        gc->GetYaxis()->SetTitleOffset(1.9);
        gc->GetXaxis()->SetLimits(0.0, 0.50);
        gc->SetTitle("");

        gc->Draw("AP"); // embedded fit function draws automatically

        TF1 *fitFn = GetGraphFit(gc);
        if (fitFn) {
          fitFn->SetLineColor(kBlue);
          fitFn->SetLineWidth(3);
        }

        // vertical reference at x=0.30 to mark the fit boundary
        double ylo = gc->GetHistogram()->GetMinimum();
        double yhi = gc->GetHistogram()->GetMaximum();
        TLine *vl = new TLine(0.30, ylo, 0.30, yhi);
        vl->SetLineStyle(2);
        vl->SetLineColor(kGray + 2);
        vl->SetLineWidth(1);
        vl->Draw();

        // horizontal reference at y=1
        TLine *hl = new TLine(0.0, 1.0, 0.50, 1.0);
        hl->SetLineStyle(2);
        hl->SetLineColor(kGray + 2);
        hl->SetLineWidth(1);
        hl->Draw();

        DrawInfoLegend(0.62, 0.68, 0.92, 0.90,
                       {cone, kMethodLabels[m], ptSl.title,
                        Form("|#eta| bin %d", ie), CalibTag(useJer)});

        SavePlot(c, outDir, cone, "alpha", {calibKey, etaKey, ptKey}, cvName);
        pb.Update();

        delete gc;
        delete c;
      }
    }
  }

  // combined overlay: all 4 methods together, same (pT slice, eta bin), same dir
  for (int ip = 0; ip < nPt; ip++) {
    const auto &ptSl = bins.ptavgSlices[ip];
    for (int ie = 0; ie < nEta; ie++) {
      TString etaKey = L2Name::EtaKey(ie, false);
      TString ptKey = L2Name::PtKey(ptSl);

      TGraphErrors *gc[kNMethods] = {};
      bool any = false;
      for (int m = 0; m < kNMethods; m++) {
        TString gname = L2Name::ObjectName(
            cone, CalibKind("R", useJer),
            {L2Name::EtaModeKey(false), etaKey, ptKey}, {kMethodKeys[m]});
        TString gnorm = gname + "_norm";
        TGraphErrors *gr = GetGraphAny(dGraphs, {gnorm, gname});
        if (gr && gr->GetN() >= 2) {
          gc[m] = (TGraphErrors *)gr->Clone(Form("%s_ovc%d", gnorm.Data(), m));
          any = true;
        }
      }
      if (!any)
        continue;

      const TString cvName =
          Form("alphafit_allmethods_%s_%s_%s_%s", cone.Data(),
               calibKey.Data(), etaKey.Data(), ptKey.Data());
      TCanvas *c = new TCanvas(cvName, "", 800, 800);
      RealAspectRatio(c);
      c->SetLeftMargin(0.18);
      c->SetGridx();
      c->SetGridy();

      // combined y-range across whichever methods are present
      double ylo = 1e9, yhi = -1e9;
      for (int m = 0; m < kNMethods; m++) {
        if (!gc[m])
          continue;
        for (int i = 0; i < gc[m]->GetN(); i++) {
          double px, py;
          gc[m]->GetPoint(i, px, py);
          const double pey = gc[m]->GetErrorY(i);
          ylo = std::min(ylo, py - pey);
          yhi = std::max(yhi, py + pey);
        }
      }
      if (ylo > yhi) {
        ylo = 0.99;
        yhi = 1.01;
      }
      const double pad = 0.15 * (yhi - ylo);
      ylo -= pad;
      yhi += pad;

      TLegend *leg = new TLegend(0.16, 0.16, 0.50, 0.16 + 0.05 * kNMethods);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.030);

      bool first = true;
      for (int m = 0; m < kNMethods; m++) {
        if (!gc[m])
          continue;
        gc[m]->SetMarkerStyle(kMethodStyles[m]);
        gc[m]->SetMarkerColor(kMethodColors[m]);
        gc[m]->SetLineColor(kMethodColors[m]);
        gc[m]->SetMarkerSize(0.9);
        gc[m]->SetTitle("");

        if (first) {
          gc[m]->GetXaxis()->SetTitle("#alpha threshold");
          gc[m]->GetYaxis()->SetTitle(kAlphaYTitle);
          gc[m]->GetXaxis()->CenterTitle();
          gc[m]->GetYaxis()->CenterTitle();
          gc[m]->GetYaxis()->SetTitleOffset(1.9);
          gc[m]->GetXaxis()->SetLimits(0.0, 0.50);
          gc[m]->Draw("AP");
          gc[m]->GetHistogram()->SetMinimum(ylo);
          gc[m]->GetHistogram()->SetMaximum(yhi);
        } else {
          gc[m]->Draw("P");
        }

        TF1 *fitFn = GetGraphFit(gc[m]);
        if (fitFn) {
          fitFn->SetLineColor(kMethodColors[m]);
          fitFn->SetLineStyle(1);
          fitFn->SetLineWidth(2);
          fitFn->Draw("same");
        }

        leg->AddEntry(gc[m], kMethodLabels[m], "lp");
        first = false;
      }
      leg->Draw();

      // vertical reference at x=0.30 to mark the fit boundary
      TLine *vl = new TLine(0.30, ylo, 0.30, yhi);
      vl->SetLineStyle(2);
      vl->SetLineColor(kGray + 2);
      vl->SetLineWidth(1);
      vl->Draw();

      // horizontal reference at y=1
      TLine *hl = new TLine(0.0, 1.0, 0.50, 1.0);
      hl->SetLineStyle(2);
      hl->SetLineColor(kGray + 2);
      hl->SetLineWidth(1);
      hl->Draw();

      DrawInfoLegend(
          0.62, 0.75, 0.92, 0.90,
          {cone, ptSl.title, Form("|#eta| bin %d", ie), CalibTag(useJer)});

      SavePlot(c, outDir, cone, "alpha", {calibKey, etaKey, ptKey}, cvName);
      pb.Update();

      for (auto *g : gc)
        delete g;
      delete c;
    }
  }
}

#endif
