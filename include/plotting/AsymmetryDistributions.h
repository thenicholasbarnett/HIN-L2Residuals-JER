#ifndef L2RESIDUALS_PLOTTING_ASYMMETRY_DISTRIBUTIONS_H
#define L2RESIDUALS_PLOTTING_ASYMMETRY_DISTRIBUTIONS_H

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TString.h"

#include "plotting/Style.h"
#include "plotting/Utilities.h"
#include "Binning.h"
#include "Naming.h"
#include "ProgressBar.h"

#include <algorithm>
#include <utility>

// Return the x-values at the low and high truncation boundaries for `fraction` of the area.
inline std::pair<double, double> TruncBounds(TH1D *h, double fraction) {
  if (!h || h->Integral() <= 0)
    return {0, 0};
  double total = h->Integral();
  double tailN = 0.5 * (1.0 - fraction) * total;
  int nBins = h->GetNbinsX();
  double running = 0;
  int binLo = 1;
  for (int b = 1; b <= nBins; b++) {
    running += h->GetBinContent(b);
    if (running > tailN) {
      binLo = b;
      break;
    }
  }
  running = 0;
  int binHi = nBins;
  for (int b = nBins; b >= 1; b--) {
    running += h->GetBinContent(b);
    if (running > tailN) {
      binHi = b;
      break;
    }
  }
  return {h->GetBinLowEdge(binLo),
          h->GetBinLowEdge(binHi) + h->GetBinWidth(binHi)};
}

inline void NormalizeDensity(TH1D *h) {
  if (!h)
    return;
  const double norm = h->Integral("width");
  if (norm > 0)
    h->Scale(1.0 / norm);
}

inline void StyleGuideLine(TLine *l, Color_t col) {
  l->SetLineColorAlpha(col, 0.65);
  l->SetLineStyle(3);
  l->SetLineWidth(1);
}

inline void DrawTruncLines(TH1D *h, double fraction, Color_t col, double yMin,
                           double yMax, TLine *&loLine, TLine *&hiLine) {
  auto [xlo, xhi] = TruncBounds(h, fraction);
  loLine = new TLine(xlo, yMin, xlo, yMax);
  hiLine = new TLine(xhi, yMin, xhi, yMax);
  StyleGuideLine(loLine, col);
  StyleGuideLine(hiLine, col);
  loLine->Draw();
  hiLine->Draw();
}

inline TF1 *FitGaussianGuide(TH1D *h, const TString &name, Color_t col,
                             int minEntries) {
  if (!h || h->GetEntries() < minEntries)
    return nullptr;
  TF1 *fit = new TF1(name, "gaus", -0.5, 0.5);
  fit->SetParameter(0, h->GetMaximum());
  fit->SetParameter(1, h->GetMean());
  fit->SetParameter(2, std::max(h->GetRMS(), 1e-3));
  StyleFit(fit, col);
  h->Fit(fit, "NQSR");
  return fit;
}

// A = (p_T^probe - p_T^tag) / (p_T^probe + p_T^tag), see Dijet.h::MakeDijet
static const TString kAsymXTitle =
    "#frac{p_{T}^{probe} - p_{T}^{tag}}{p_{T}^{probe} + p_{T}^{tag}}";

inline void DrawAsymBase(TH1D *hData, TH1D *hMC, const TString &xTitle,
                         double yMin, double yMax) {
  // both open circles -- distinguished by color, not marker fill
  StyleH(hData, kBlue, 24, 1.5f);
  StyleH(hMC, kRed, 24, 1.5f);

  hData->SetTitle("");
  hData->GetXaxis()->SetTitle(xTitle);
  hData->GetXaxis()->CenterTitle();
  hData->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{dA}");
  hData->GetYaxis()->CenterTitle();
  hData->GetYaxis()->SetTitleOffset(1.5);
  hData->SetMinimum(yMin);
  hData->SetMaximum(yMax);

  hData->Draw("E1");
  hMC->Draw("E1 same");
}

// Asymmetry distributions. Three canvases per (cone, pT slice, alpha slice,
// eta bin), data (blue) vs MC (red), log-y: trunc90/trunc95 truncation
// bounds, gauss fit guides. Skips bins under minEntries -- same threshold
// ([cuts] min_entries_per_bin) extraction itself requires to fit a bin, so
// plotting never re-fits something extraction already wrote off as too thin.

inline void PlotAsymDist(TFile *fIn, const TString &outDir, const TString &cone,
                         const BinningConfig &bins, int minEntries,
                         ProgressBar &pb) {
  TDirectory *dData = (TDirectory *)fIn->Get(cone + "/QA_data");
  TDirectory *dMC = (TDirectory *)fIn->Get(cone + "/QA_mc");

  const int nPt = (int)bins.ptavgSlices.size();
  const int nAlpha = (int)bins.alphaSlices.size();
  const int nEta = (int)kAbsEtaEdges.size() - 1;

  for (int ip = 0; ip < nPt; ip++) {
    const auto &ptSl = bins.ptavgSlices[ip];
    for (int ia = 0; ia < nAlpha; ia++) {
      const auto &aSl = bins.alphaSlices[ia];
      for (int ie = 0; ie < nEta; ie++) {
        if (!pb.ShouldKeep())
          continue;
        TString etaKey = L2Name::EtaKey(ie, false);
        TString ptKey = L2Name::PtKey(ptSl);
        TString alphaKey = L2Name::AlphaKey(aSl);
        TString dname = L2Name::ObjectName(
            cone, "A_data",
            {L2Name::EtaModeKey(false), etaKey, ptKey, alphaKey});
        TString mname = L2Name::ObjectName(
            cone, "A_mc", {L2Name::EtaModeKey(false), etaKey, ptKey, alphaKey});

        TH1D *hd = GetHAny(dData, {dname});
        TH1D *hm = GetHAny(dMC, {mname});

        if (!hd || !hm || hd->GetEntries() < minEntries) {
          if (hd)
            delete hd;
          if (hm)
            delete hm;
          continue;
        }

        Long64_t nData = (Long64_t)hd->GetEntries();
        Long64_t nMC = (Long64_t)hm->GetEntries();

        TH1D *hdc = (TH1D *)hd->Clone(dname + "_c");
        hdc->SetDirectory(0);
        TH1D *hmc = (TH1D *)hm->Clone(mname + "_c");
        hmc->SetDirectory(0);

        NormalizeDensity(hdc);
        NormalizeDensity(hmc);

        const double ymax =
            std::max(hdc->GetMaximum(), hmc->GetMaximum()) * 5.0;
        const double ymin = 1e-3;

        double etalo = kAbsEtaEdges[ie];
        double etahi = kAbsEtaEdges[ie + 1];

        auto drawInfo = [&]() {
          TLatex *tex = new TLatex();
          tex->SetNDC();
          tex->SetTextSize(0.031);
          tex->SetTextFont(42);
          tex->SetTextAlign(31);
          tex->DrawLatex(0.93, 0.840, cone);
          tex->DrawLatex(0.93, 0.785,
                         Form("%.3f < |#eta^{probe}| < %.3f", etalo, etahi));
          tex->DrawLatex(0.93, 0.730, ptSl.title.Data());
          tex->DrawLatex(0.93, 0.675, aSl.title.Data());
        };

        auto makeLegend = [&]() {
          TLegend *leg = new TLegend(0.16, 0.66, 0.48, 0.86);
          leg->SetBorderSize(0);
          leg->SetFillColorAlpha(kWhite, 0.75);
          leg->SetFillStyle(1001);
          leg->SetTextAlign(12);
          leg->SetTextSize(0.031);
          leg->AddEntry(hdc, Form("Data (%s)", FormatEntriesText(nData).Data()),
                        "lp");
          leg->AddEntry(hmc, Form("MC (%s)", FormatEntriesText(nMC).Data()),
                        "lp");
          return leg;
        };

        auto drawTruncPlot = [&](double fraction, const TString &tag,
                                 const TString &label) {
          TString cvName =
              Form("adist_%s_%s_%s_%s_%s", cone.Data(), etaKey.Data(),
                   ptKey.Data(), alphaKey.Data(), tag.Data());
          TCanvas *c = new TCanvas(cvName, "", 800, 800);
          RealAspectRatio(c);
          c->SetLogy();
          c->SetLeftMargin(0.15);
          c->SetRightMargin(0.05);

          DrawAsymBase(hdc, hmc, kAsymXTitle, ymin, ymax);
          DrawAsymHeader();
          drawInfo();

          TLine *ldLo = nullptr, *ldHi = nullptr, *lmLo = nullptr,
                *lmHi = nullptr;
          DrawTruncLines(hdc, fraction, kBlue, ymin, ymax, ldLo, ldHi);
          DrawTruncLines(hmc, fraction, kRed, ymin, ymax, lmLo, lmHi);

          TLegend *leg = makeLegend();
          leg->AddEntry(ldLo, Form("Data %s", label.Data()), "l");
          leg->AddEntry(lmLo, Form("MC %s", label.Data()), "l");
          leg->Draw();

          SavePlot(c, outDir, cone, "adist", {etaKey, ptKey, alphaKey}, cvName);
          pb.Update();
          delete c;
        };

        drawTruncPlot(0.90, "trunc90", "trunc. 90%");
        drawTruncPlot(0.95, "trunc95", "trunc. 95%");

        {
          TString cvName = Form("adist_%s_%s_%s_%s_gauss", cone.Data(),
                                etaKey.Data(), ptKey.Data(), alphaKey.Data());
          TCanvas *c = new TCanvas(cvName, "", 800, 800);
          RealAspectRatio(c);
          c->SetLogy();
          c->SetLeftMargin(0.15);
          c->SetRightMargin(0.05);

          DrawAsymBase(hdc, hmc, kAsymXTitle, ymin, ymax);
          DrawAsymHeader();
          drawInfo();

          TF1 *fd =
              FitGaussianGuide(hdc, cvName + "_data_fit", kBlue, minEntries);
          TF1 *fm =
              FitGaussianGuide(hmc, cvName + "_mc_fit", kRed, minEntries);
          if (fd)
            fd->Draw("same");
          if (fm)
            fm->Draw("same");

          TLegend *leg = makeLegend();
          if (fd)
            leg->AddEntry(fd, "Data Gaussian fit", "l");
          if (fm)
            leg->AddEntry(fm, "MC Gaussian fit", "l");
          leg->Draw();

          SavePlot(c, outDir, cone, "adist", {etaKey, ptKey, alphaKey}, cvName);
          pb.Update();
          delete fd;
          delete fm;
          delete c;
        }

        delete hdc;
        delete hmc;
        delete hd;
        delete hm;
      }
    }
  }
}

#endif
