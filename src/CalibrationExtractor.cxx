#include "CalibrationExtractor.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "THnSparse.h"
#include "TDirectory.h"
#include "TString.h"
#include "TMath.h"

#include "Binning.h"
#include "Naming.h"
#include "Utilities.h"
#include "ProgressBar.h"
#include "AnalysisConfig.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

static constexpr int kMaxOutOfRangeMeanWarnings = 20;

// axes
static constexpr int kEtaAxis = 0;
static constexpr int kPtAvgAxis = 1;
static constexpr int kAlphaAxis = 2;
static constexpr int kAAxis = 3;

// methods to capture mean
static constexpr int kNMethods = 4;
static const char *const kMethodNames[] = {"gauss", "doubleGauss", "trunc90",
                                           "trunc95"};

struct GaussResult {
  double mean = 0, meanErr = 0, sigma = 0, sigmaErr = 0, chi2ndf = -1;
  bool valid = false;
};

struct TruncResult {
  double mean = 0, meanErr = 0, sigma = 0, sigmaErr = 0, nEff = 0;
  bool valid = false;
};

struct ResidualFitControls {
  int minEntries = 0;
  double gausFitHalfWidth = 0.0;
  double maxAbsA = 0.0;
  double alphaFitHi = 0.0;
  int outOfRangeMeanWarnings = 0;
};

// for alpha extrapolations
struct RPoint {
  double alpha, val, err;
};

// gauss fitting
static GaussResult FitGauss(TH1D *h, const ResidualFitControls &controls) {
  GaussResult r;
  if (!CanFit(h, controls.minEntries))
    return r;

  TF1 *g = new TF1(Form("_gf_%s", h->GetName()), "gaus",
                   -controls.gausFitHalfWidth, controls.gausFitHalfWidth);
  g->SetParameter(0, h->GetMaximum());
  g->SetParameter(1, h->GetMean());
  g->SetParameter(2, std::max(h->GetRMS(), 1e-3));

  TFitResultPtr res = h->Fit(g, "NQSR");
  if (res.Get() && res->IsValid()) {
    r.mean = res->Parameter(1);
    r.meanErr = res->ParError(1);
    r.sigma = res->Parameter(2);
    r.sigmaErr = res->ParError(2);
    r.chi2ndf = (res->Ndf() > 0) ? res->Chi2() / res->Ndf() : -1.0;
    r.valid = true;
  }
  delete g;
  return r;
}

// double gauss fitting
static GaussResult FitDoubleGauss(TH1D *h,
                                  const ResidualFitControls &controls) {
  GaussResult r;
  if (!CanFit(h, controls.minEntries))
    return r;

  TF1 *g = new TF1(Form("_dgf_%s", h->GetName()), "gaus(0)+gaus(3)",
                   -controls.gausFitHalfWidth, controls.gausFitHalfWidth);
  const double rms = std::max(h->GetRMS(), 1e-3);
  g->SetParameter(0, h->GetMaximum());  // core amplitude
  g->SetParameter(1, h->GetMean());  // core mean
  g->SetParameter(2, 0.5 * rms);  // core sigma
  g->SetParameter(3, 0.3 * h->GetMaximum());  // tail amplitude
  g->SetParameter(4, h->GetMean());  // tail mean
  g->SetParameter(5, 1.5 * rms);  // tail sigma 

  TFitResultPtr res = h->Fit(g, "NQSR");
  if (res.Get() && res->IsValid()) {
    const double amp1 = res->Parameter(0), mean1 = res->Parameter(1),
                 sigma1 = res->Parameter(2);
    const double amp2 = res->Parameter(3), mean2 = res->Parameter(4),
                 sigma2 = res->Parameter(5);
    const double w1 = std::fabs(amp1 * sigma1);
    const double w2 = std::fabs(amp2 * sigma2);
    const double wsum = w1 + w2;
    if (wsum > 0) {
      const double f1 = w1 / wsum, f2 = w2 / wsum;
      r.mean = f1 * mean1 + f2 * mean2;
      // errors summed in quadrature
      r.meanErr = std::sqrt(TMath::Power(f1 * res->ParError(1), 2) +
                            TMath::Power(f2 * res->ParError(4), 2));
      // weighted sigma
      r.sigma = std::sqrt(std::max(
          0.0, f1 * (sigma1 * sigma1 + mean1 * mean1) +
                   f2 * (sigma2 * sigma2 + mean2 * mean2) - r.mean * r.mean));
      r.sigmaErr = std::sqrt(TMath::Power(f1 * res->ParError(2), 2) +
                             TMath::Power(f2 * res->ParError(5), 2));
      r.chi2ndf = (res->Ndf() > 0) ? res->Chi2() / res->Ndf() : -1.0;
      r.valid = true;
    }
  }
  delete g;
  return r;
}

// truncation
static std::pair<int, int> FindTruncBins(TH1D *h, double fraction,
                                         const ResidualFitControls &controls) {
  if (!CanFit(h, controls.minEntries))
    return {1, 0};
  double total = h->Integral();
  if (total <= 0)
    return {1, 0};
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
  return {binLo, binHi};
}

// mean and width of h in bins [binLo, binHi]
static TruncResult TruncMeanInRange(TH1D *h, int binLo, int binHi) {
  TruncResult r;
  if (binLo > binHi)
    return r;
  double nEff = h->Integral(binLo, binHi);
  double hTotal = h->Integral();
  double nEffEntries =
      (hTotal > 1e-10) ? h->GetEntries() * nEff / hTotal : nEff;
  // no separate entries gate here: FindTruncBins already required
  // h's raw (pre-truncation) entries >= min_entries_per_bin
  h->GetXaxis()->SetRange(binLo, binHi);
  double mean = h->GetMean();
  double rms = h->GetRMS();
  h->GetXaxis()->SetRange(0, 0);
  r.mean = mean;
  r.meanErr = rms / TMath::Sqrt(nEffEntries);
  r.sigma = rms;
  // large-N (asymptotic normal) approximation for the SE of a sample std-dev estimate
  // same tradeoff as meanErr above, not a full jackknife/bootstrap estimate
  r.sigmaErr = rms / TMath::Sqrt(2.0 * TMath::Max(nEffEntries - 1.0, 1.0));
  r.nEff = nEffEntries;
  r.valid = true;
  return r;
}

// A (avg or sigma) → R
static double ToR(double A) { return (1.0 + A) / (1.0 - A); }
static double ToRErr(double A, double dA) {
  return dA * 2.0 / ((1.0 - A) * (1.0 - A));
}

static void ResetRange(THnSparse *h, int axis) {
  h->GetAxis(axis)->SetRange(0, 0);
}

static void
WarnOutOfRangeValue(ResidualFitControls &controls, const TString &cone,
                    const TString &etaMode, const TString &etaKey,
                    const TString &ptKey, const TString &alphaKey,
                    const char *method, const char *sample,
                    const char *quantity, // "mean A" (JEC) or "stddev A" (JER)
                    double value) {
  if (controls.outOfRangeMeanWarnings < kMaxOutOfRangeMeanWarnings) {
    std::cerr << "WARNING: residual " << method << " " << sample << " "
              << quantity << "=" << value
              << " exceeds " << controls.maxAbsA << " for "
              << cone << " " << etaMode << " " << etaKey << " " << ptKey << " "
              << alphaKey << "; dropping this point, but please investigate this! (unexpected behaviour and a red flag) \n";
  } else if (controls.outOfRangeMeanWarnings == kMaxOutOfRangeMeanWarnings) {
    std::cerr << "WARNING: more values exceed configured max_abs_a="
              << controls.maxAbsA
              << "; suppressing more of these warnings\n";
  }
  controls.outOfRangeMeanWarnings++;
}

// direct correction and kFSR extraction
// JEC (mean-derived R points)
// JER (fed stddev-derived R points)
// each (alpha, ptSlice, etaBin) fit returns mean and width of A

struct ExtrapResult {
  double c0 = 0, ec0 = 0;  // direct intercept
  double c0n = 0, ec0n = 0;  // kFSR (normalized alpha values by nominal alpha)
  bool valid = false;
  bool normValid = false;
};

static ExtrapResult FitAndExtrapolate(const std::vector<RPoint> &pts,
                                      const TString &gname, Color_t color,
                                      const ResidualFitControls &controls,
                                      TDirectory *dGraphs) {

  ExtrapResult out;
  const int n = (int)pts.size();
  std::vector<double> x(n), y(n), ex(n, 0.0), ey(n);
  for (int k = 0; k < n; k++) {
    x[k] = pts[k].alpha;
    y[k] = pts[k].val;
    ey[k] = pts[k].err;
  }

  TGraphErrors *gr =
      new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data());
  gr->SetName(gname);
  gr->SetTitle(";#alpha threshold;R_{MC}/R_{data}");
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(color);
  gr->SetLineColor(color);

  if (!CanFit(gr, {0.0, controls.alphaFitHi}, 2)) {
    delete gr;
    return out;
  }

  // only points within the fit range enter chi2
  TF1 *fitFn = new TF1(gname + "_fit", "[0]+[1]*x", 0.0, controls.alphaFitHi);
  fitFn->SetParameter(0, 1.0);
  fitFn->SetParameter(1, 0.0);
  fitFn->SetLineColor(color);
  gr->Fit(fitFn, "QR");

  out.c0 = fitFn->GetParameter(0);
  out.ec0 = fitFn->GetParError(0);
  out.valid = true;

  // kFSR
  // fitting values at each alpha normalized by the value at a nominal alpha
  // extrapolated value multiplied by value at nominal alpha gives correction
  double normVal = 0, normErr = 0;
  for (int k = n - 1; k >= 0; k--) {
    if (pts[k].alpha <= controls.alphaFitHi + 1e-4 && pts[k].val > 1e-6) {
      normVal = pts[k].val;
      normErr = pts[k].err;
      break;
    }
  }
  if (normVal > 1e-6) {
    int nfit = 0;
    for (int k = 0; k < n && pts[k].alpha <= controls.alphaFitHi + 1e-4; k++)
      nfit++;

    std::vector<double> xn(nfit), yn(nfit), exn(nfit, 0.0), eyn(nfit);
    bool bad = false;
    for (int k = 0; k < nfit; k++) {
      if (std::abs(pts[k].val) < 1e-6) {
        bad = true;
        break;
      }
      xn[k] = pts[k].alpha;
      yn[k] = pts[k].val / normVal;
      eyn[k] = yn[k] * TMath::Sqrt(TMath::Power(pts[k].err / pts[k].val, 2.0) +
                                   TMath::Power(normErr / normVal, 2.0));
    }
    if (!bad && nfit >= 2) {
      TString gnorm = gname + "_norm";
      TGraphErrors *grn =
          new TGraphErrors(nfit, xn.data(), yn.data(), exn.data(), eyn.data());
      grn->SetName(gnorm);
      grn->SetTitle(";#alpha threshold;R_{MC}/R_{data} (k_{FSR})");
      grn->SetMarkerStyle(20);
      grn->SetMarkerColor(color);
      grn->SetLineColor(color);
      TF1 *fn = new TF1(gnorm + "_fit", "[0]+[1]*x", 0.0, controls.alphaFitHi);
      fn->SetParameter(0, 1.0);
      fn->SetParameter(1, 0.0);
      fn->SetLineColor(color);
      grn->Fit(fn, "QR");

      const double c0n = fn->GetParameter(0);
      const double ec0n = fn->GetParError(0);
      out.c0n = c0n * normVal;
      out.ec0n =
          (std::abs(c0n) > 1e-9)
              ? out.c0n * TMath::Sqrt(TMath::Power(ec0n / c0n, 2.0) +
                                      TMath::Power(normErr / normVal, 2.0))
              : ec0n * normVal;
      out.normValid = true;

      dGraphs->cd();
      grn->Write();
      delete fn;
      delete grn;
    }
  }

  dGraphs->cd();
  gr->Write();  // fit function clone
  delete fitFn;  // delete original
  delete gr;
  return out;
}

// full extraction loop for JER SF or L2Residuals
// progress bar updated once per pT slice
// namespace: cone, object kind, eta mode/bin, ptavg, alpha, method
//
// Outputs: 
// dGraph (method, ptSlice, etaBin)
// TGraphErrors of R_MC/R_data vs alpha (R_jer)

static void ExtractAndFit(THnSparse *hData, THnSparse *hMC, const TString &cone,
                          const BinningConfig &bins, TDirectory *dQA_data,
                          TDirectory *dQA_mc, TDirectory *dGraphs,
                          TDirectory *dRvals, TDirectory *dOut,
                          const std::vector<Double_t> &etaEdges,
                          const TString &nameSuffix,
                          ResidualFitControls &controls, CalibrationMode mode,
                          ProgressBar &pb) {
  const bool doJEC = (mode == CalibrationMode::JEC);
  const bool doJER = (mode == CalibrationMode::JER);
  const int nPt = (int)bins.ptavgSlices.size();
  const int nAlpha = (int)bins.alphaSlices.size();
  const int nEta = hData->GetAxis(kEtaAxis)->GetNbins();
  const bool fullEta = !nameSuffix.IsNull();
  const TString etaMode = L2Name::EtaModeKey(fullEta);

  // rpts[method][ipt][ieta]: all nAlpha bins
  // rptsJer mirrors rpts, given stddev R points instead of mean-derived
  std::vector<std::vector<std::vector<std::vector<RPoint>>>> rpts(
      kNMethods, std::vector<std::vector<std::vector<RPoint>>>(
                     nPt, std::vector<std::vector<RPoint>>(nEta)));
  std::vector<std::vector<std::vector<std::vector<RPoint>>>> rptsJer(
      kNMethods, std::vector<std::vector<std::vector<RPoint>>>(
                     nPt, std::vector<std::vector<RPoint>>(nEta)));

  // R_data and R_mc vs eta, TH1Ds (method, ialpha, ipt)
  using RH = std::vector<std::vector<std::vector<TH1D *>>>;
  RH hRd(kNMethods, std::vector<std::vector<TH1D *>>(
                        nAlpha, std::vector<TH1D *>(nPt, nullptr)));
  RH hRm(kNMethods, std::vector<std::vector<TH1D *>>(
                        nAlpha, std::vector<TH1D *>(nPt, nullptr)));
  RH hRdJer(kNMethods, std::vector<std::vector<TH1D *>>(
                           nAlpha, std::vector<TH1D *>(nPt, nullptr)));
  RH hRmJer(kNMethods, std::vector<std::vector<TH1D *>>(
                           nAlpha, std::vector<TH1D *>(nPt, nullptr)));

  for (int m = 0; m < kNMethods; m++) {
    for (int ia = 0; ia < nAlpha; ia++) {
      for (int ip = 0; ip < nPt; ip++) {
        const std::vector<TString> keys = {
            etaMode, L2Name::PtKey(bins.ptavgSlices[ip]),
            L2Name::AlphaKey(bins.alphaSlices[ia])};
        hRd[m][ia][ip] = new TH1D(
            L2Name::ObjectName(cone, "R_data", keys, {kMethodNames[m]}), "",
            (int)etaEdges.size() - 1, etaEdges.data());
        hRm[m][ia][ip] =
            new TH1D(L2Name::ObjectName(cone, "R_mc", keys, {kMethodNames[m]}),
                     "", (int)etaEdges.size() - 1, etaEdges.data());
        hRdJer[m][ia][ip] = new TH1D(
            L2Name::ObjectName(cone, "R_data_jer", keys, {kMethodNames[m]}), "",
            (int)etaEdges.size() - 1, etaEdges.data());
        hRmJer[m][ia][ip] = new TH1D(
            L2Name::ObjectName(cone, "R_mc_jer", keys, {kMethodNames[m]}), "",
            (int)etaEdges.size() - 1, etaEdges.data());
      }
    }
  }

  // main extraction loop
  for (int ipt = 0; ipt < nPt; ipt++) {
    const auto &ptSlice = bins.ptavgSlices[ipt];

    for (int ialpha = 0; ialpha < nAlpha; ialpha++) {
      const auto &aSlice = bins.alphaSlices[ialpha];

      hData->GetAxis(kPtAvgAxis)->SetRangeUser(ptSlice.lo, ptSlice.hi);
      hMC->GetAxis(kPtAvgAxis)->SetRangeUser(ptSlice.lo, ptSlice.hi);
      hData->GetAxis(kAlphaAxis)->SetRangeUser(aSlice.lo, aSlice.hi);
      hMC->GetAxis(kAlphaAxis)->SetRangeUser(aSlice.lo, aSlice.hi);

      // 2D (eta, A) projection per (ipt, ialpha)
      // Projection(yDim, xDim) → TH2D with x=eta, y=A.
      TH2D *h2Data;
      TH2D *h2MC;
      {
        TDirectory::TContext nodir(nullptr);
        h2Data = hData->Projection(kAAxis, kEtaAxis);
        h2MC = hMC->Projection(kAAxis, kEtaAxis);
      }

      for (int ieta = 0; ieta < nEta; ieta++) {
        TString etaKey = L2Name::EtaKey(ieta, fullEta);
        TString ptKey = L2Name::PtKey(ptSlice);
        TString alphaKey = L2Name::AlphaKey(aSlice);

        TH1D *hAData;
        TH1D *hAMC;
        {
          TDirectory::TContext nodir(nullptr);
          hAData = h2Data->ProjectionY(
              L2Name::ObjectName(cone, "A_data",
                                 {etaMode, etaKey, ptKey, alphaKey}),
              ieta + 1, ieta + 1);
          hAMC = h2MC->ProjectionY(
              L2Name::ObjectName(cone, "A_mc",
                                 {etaMode, etaKey, ptKey, alphaKey}),
              ieta + 1, ieta + 1);
        }

        dQA_data->cd();
        hAData->Write();
        dQA_mc->cd();
        hAMC->Write();

        GaussResult gd = FitGauss(hAData, controls);
        GaussResult gm = FitGauss(hAMC, controls);
        GaussResult ddg = FitDoubleGauss(hAData, controls);
        GaussResult mdg = FitDoubleGauss(hAMC, controls);
        auto [dlo90, dhi90] = FindTruncBins(hAData, 0.90, controls);
        auto [mlo90, mhi90] = FindTruncBins(hAMC, 0.90, controls);
        auto [dlo95, dhi95] = FindTruncBins(hAData, 0.95, controls);
        auto [mlo95, mhi95] = FindTruncBins(hAMC, 0.95, controls);
        TruncResult td90 = TruncMeanInRange(hAData, dlo90, dhi90);
        TruncResult tm90 = TruncMeanInRange(hAMC, mlo90, mhi90);
        TruncResult td95 = TruncMeanInRange(hAData, dlo95, dhi95);
        TruncResult tm95 = TruncMeanInRange(hAMC, mlo95, mhi95);

        double alphaX = aSlice.hi;

        // accumulate alpha bins
        // L2Res/JEC=fed mean/meanErr, JER(SF)=fed sigma/sigmaErr
        // R=(1+x)/(1-x)
        auto accumulate =
            [&](std::vector<std::vector<std::vector<std::vector<RPoint>>>>
                    &rptsOut,
                RH &hRdOut, RH &hRmOut, const char *quantity, int method,
                double Ad, double eAd, double Am, double eAm, bool ok) {
              if (!ok)
                return;
              bool outsideConfiguredRange = false;
              if (std::abs(Ad) > controls.maxAbsA) {
                WarnOutOfRangeValue(controls, cone, etaMode, etaKey, ptKey,
                                    alphaKey, kMethodNames[method], "data",
                                    quantity, Ad);
                outsideConfiguredRange = true;
              }
              if (std::abs(Am) > controls.maxAbsA) {
                WarnOutOfRangeValue(controls, cone, etaMode, etaKey, ptKey,
                                    alphaKey, kMethodNames[method], "MC",
                                    quantity, Am);
                outsideConfiguredRange = true;
              }
              if (outsideConfiguredRange)
                return;
              double Rd = ToR(Ad), Rm = ToR(Am);
              double eRd = ToRErr(Ad, eAd), eRm = ToRErr(Am, eAm);
              if (std::abs(Rd) < 1e-6)
                return;
              double ratio = Rm / Rd;
              double eRatio = ratio * TMath::Sqrt((eRd / Rd) * (eRd / Rd) +
                                                  (eRm / Rm) * (eRm / Rm));
              rptsOut[method][ipt][ieta].push_back({alphaX, ratio, eRatio});
              hRdOut[method][ialpha][ipt]->SetBinContent(ieta + 1, Rd);
              hRdOut[method][ialpha][ipt]->SetBinError(ieta + 1, eRd);
              hRmOut[method][ialpha][ipt]->SetBinContent(ieta + 1, Rm);
              hRmOut[method][ialpha][ipt]->SetBinError(ieta + 1, eRm);
            };

        if (doJEC) {
          accumulate(rpts, hRd, hRm, "mean A", 0, gd.mean, gd.meanErr, gm.mean,
                     gm.meanErr, gd.valid && gm.valid);
          accumulate(rpts, hRd, hRm, "mean A", 1, ddg.mean, ddg.meanErr,
                     mdg.mean, mdg.meanErr, ddg.valid && mdg.valid);
          accumulate(rpts, hRd, hRm, "mean A", 2, td90.mean, td90.meanErr,
                     tm90.mean, tm90.meanErr, td90.valid && tm90.valid);
          accumulate(rpts, hRd, hRm, "mean A", 3, td95.mean, td95.meanErr,
                     tm95.mean, tm95.meanErr, td95.valid && tm95.valid);
        }
        if (doJER) {
          accumulate(rptsJer, hRdJer, hRmJer, "stddev A", 0, gd.sigma,
                     gd.sigmaErr, gm.sigma, gm.sigmaErr, gd.valid && gm.valid);
          accumulate(rptsJer, hRdJer, hRmJer, "stddev A", 1, ddg.sigma,
                     ddg.sigmaErr, mdg.sigma, mdg.sigmaErr,
                     ddg.valid && mdg.valid);
          accumulate(rptsJer, hRdJer, hRmJer, "stddev A", 2, td90.sigma,
                     td90.sigmaErr, tm90.sigma, tm90.sigmaErr,
                     td90.valid && tm90.valid);
          accumulate(rptsJer, hRdJer, hRmJer, "stddev A", 3, td95.sigma,
                     td95.sigmaErr, tm95.sigma, tm95.sigmaErr,
                     td95.valid && tm95.valid);
        }

        delete hAData;
        delete hAMC;
      }

      delete h2Data;
      delete h2MC;

      ResetRange(hData, kPtAvgAxis);
      ResetRange(hMC, kPtAvgAxis);
      ResetRange(hData, kAlphaAxis);
      ResetRange(hMC, kAlphaAxis);
    }

    pb.Update();
  }

  // write R_data and R_mc histograms
  dRvals->cd();
  for (int m = 0; m < kNMethods; m++) {
    for (int ia = 0; ia < nAlpha; ia++) {
      for (int ip = 0; ip < nPt; ip++) {
        if (doJEC) {
          hRd[m][ia][ip]->Write();
          hRm[m][ia][ip]->Write();
        }
        if (doJER) {
          hRdJer[m][ia][ip]->Write();
          hRmJer[m][ia][ip]->Write();
        }
        delete hRd[m][ia][ip];
        delete hRm[m][ia][ip];
        delete hRdJer[m][ia][ip];
        delete hRmJer[m][ia][ip];
      }
    }
  }

  // TGraphErrors and fit R_ratio vs alpha
  for (int method = 0; method < kNMethods; method++) {
    for (int ipt = 0; ipt < nPt; ipt++) {
      const auto &ptSlice = bins.ptavgSlices[ipt];

      TString ptKey = L2Name::PtKey(ptSlice);
      TString corrName = L2Name::ObjectName(cone, "intercept", {etaMode, ptKey},
                                            {kMethodNames[method]});
      TString corrNameJer = L2Name::ObjectName(
          cone, "intercept_jer", {etaMode, ptKey}, {kMethodNames[method]});

      const TString etaAxisTitle = nameSuffix.IsNull() ? "|#eta|" : "#eta";
      const int nEtaBins = (int)etaEdges.size() - 1;

      TH1D *hCorr =
          doJEC ? new TH1D(corrName, "", nEtaBins, etaEdges.data()) : nullptr;
      TH1D *hCorrNorm =
          doJEC ? new TH1D(corrName + "_norm", "", nEtaBins, etaEdges.data())
                : nullptr;
      TH1D *hCorrJer =
          doJER ? new TH1D(corrNameJer, "", nEtaBins, etaEdges.data())
                : nullptr;
      TH1D *hCorrNormJer =
          doJER ? new TH1D(corrNameJer + "_norm", "", nEtaBins, etaEdges.data())
                : nullptr;

      if (hCorr) {
        hCorr->GetXaxis()->SetTitle(etaAxisTitle);
        hCorr->GetYaxis()->SetTitle("R_{MC}/R_{data} at #alpha=0");
        hCorrNorm->GetXaxis()->SetTitle(etaAxisTitle);
        hCorrNorm->GetYaxis()->SetTitle(
            "k_{FSR} #cdot R_{MC}/R_{data}|_{fit range high edge}");
      }
      if (hCorrJer) {
        hCorrJer->GetXaxis()->SetTitle(etaAxisTitle);
        hCorrJer->GetYaxis()->SetTitle(
            "JER SF: R_{MC}/R_{data} (stddev A) at #alpha=0");
        hCorrNormJer->GetXaxis()->SetTitle(etaAxisTitle);
        hCorrNormJer->GetYaxis()->SetTitle(
            "k_{FSR} #cdot JER SF|_{fit range high edge}");
      }

      for (int ieta = 0; ieta < nEta; ieta++) {
        TString etaKey = L2Name::EtaKey(ieta, fullEta);

        if (doJEC) {
          TString gname = L2Name::ObjectName(
              cone, "R", {etaMode, etaKey, ptKey}, {kMethodNames[method]});
          ExtrapResult jes = FitAndExtrapolate(
              rpts[method][ipt][ieta], gname, ptSlice.color, controls, dGraphs);
          if (jes.valid) {
            hCorr->SetBinContent(ieta + 1, jes.c0);
            hCorr->SetBinError(ieta + 1, jes.ec0);
          }
          if (jes.normValid) {
            hCorrNorm->SetBinContent(ieta + 1, jes.c0n);
            hCorrNorm->SetBinError(ieta + 1, jes.ec0n);
          }
        }
        if (doJER) {
          TString gnameJer = L2Name::ObjectName(
              cone, "R_jer", {etaMode, etaKey, ptKey}, {kMethodNames[method]});
          ExtrapResult jer =
              FitAndExtrapolate(rptsJer[method][ipt][ieta], gnameJer,
                                ptSlice.color, controls, dGraphs);
          if (jer.valid) {
            hCorrJer->SetBinContent(ieta + 1, jer.c0);
            hCorrJer->SetBinError(ieta + 1, jer.ec0);
          }
          if (jer.normValid) {
            hCorrNormJer->SetBinContent(ieta + 1, jer.c0n);
            hCorrNormJer->SetBinError(ieta + 1, jer.ec0n);
          }
        }
      }

      dOut->cd();
      if (hCorr) {
        hCorr->Write();
        hCorrNorm->Write();
      }
      if (hCorrJer) {
        hCorrJer->Write();
        hCorrNormJer->Write();
      }
      delete hCorr;
      delete hCorrNorm;
      delete hCorrJer;
      delete hCorrNormJer;
    }
  }
}

// main
void runCalibration(TString dataFile, TString mcFile, TString outputFile,
                    CalibrationMode mode) {

  const AnalysisConfig &cfg = Config();
  PrintConfigSummary(cfg);
  ResidualFitControls controls = {cfg.minEntriesPerBin,
                                  cfg.residualGausFitHalfWidth, cfg.maxAbsA,
                                  cfg.residualAlphaFitHi, 0};

  TFile *fData = TFile::Open(dataFile, "read");
  TFile *fMC = TFile::Open(mcFile, "read");
  if (!fData || fData->IsZombie()) {
    std::cerr << "Cannot open " << dataFile << "\n";
    return;
  }
  if (!fMC || fMC->IsZombie()) {
    std::cerr << "Cannot open " << mcFile << "\n";
    return;
  }

  TFile *fOut = new TFile(outputFile, "recreate");

  BinningConfig bins;
  const int nPt = (int)bins.ptavgSlices.size();

  // |eta| and full eta ExtractAndFit calls per cone, each loops through nPt pT slices
  const int totalSteps = (int)cfg.coneLabels.size() * 2 * nPt;
  ProgressBar pb("Extracting:", totalSteps);

  for (const TString &cone : cfg.coneLabels) {

    TDirectory *coneDataDir = (TDirectory *)fData->Get(cone);
    TDirectory *coneMCDir = (TDirectory *)fMC->Get(cone);
    THnSparse *hRawData =
        coneDataDir ? (THnSparse *)coneDataDir->Get(cone + "_asym") : nullptr;
    THnSparse *hRawMC =
        coneMCDir ? (THnSparse *)coneMCDir->Get(cone + "_asym") : nullptr;
    if (!hRawData) {
      std::cerr << "Missing " << cone << "_asym in data\n";
      continue;
    }
    if (!hRawMC) {
      std::cerr << "Missing " << cone << "_asym in MC\n";
      continue;
    }

    THnSparse *hData =
        FoldEtaAxis(hRawData, kEtaAxis, cone + "_asym_data_abseta");
    THnSparse *hMC = FoldEtaAxis(hRawMC, kEtaAxis, cone + "_asym_mc_abseta");

    TDirectory *coneDir = fOut->mkdir(cone.Data());
    coneDir->cd();
    hData->Write();
    hMC->Write();

    TDirectory *dQA_data = coneDir->mkdir("QA_data");
    TDirectory *dQA_mc = coneDir->mkdir("QA_mc");
    TDirectory *dQA_data_full = coneDir->mkdir("QA_data_fulleta");
    TDirectory *dQA_mc_full = coneDir->mkdir("QA_mc_fulleta");
    TDirectory *dGraphs = coneDir->mkdir("graphs");
    TDirectory *dRvals = coneDir->mkdir("Rvals");
    TDirectory *dRvals_full = coneDir->mkdir("Rvals_fulleta");

    ExtractAndFit(hData, hMC, cone, bins, dQA_data, dQA_mc, dGraphs, dRvals,
                  coneDir, kAbsEtaEdges, "", controls, mode, pb);

    ExtractAndFit(hRawData, hRawMC, cone, bins, dQA_data_full, dQA_mc_full,
                  dGraphs, dRvals_full, coneDir, kEtaEdges, "_fulleta",
                  controls, mode, pb);

    delete hData;
    delete hMC;
  }

  pb.Finish();
  fOut->Close();
  fData->Close();
  fMC->Close();
}
