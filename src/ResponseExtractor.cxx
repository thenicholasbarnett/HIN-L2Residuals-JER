#include "ResponseExtractor.h"

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFitResult.h"
#include "THnSparse.h"
#include "TDirectory.h"
#include "TString.h"
#include "TMath.h"

#include "Binning.h"
#include "CrystalBall.h"
#include "Methods.h" // kNMethods, kMethodKeys
#include "Naming.h"
#include "Utilities.h"
#include "Truncation.h"
#include "ProgressBar.h"
#include "AnalysisConfig.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

// axis order {cone}_incl_resp/_tag_resp/_probe_resp
static constexpr int kRespEtaRecoAxis = 0;
static constexpr int kRespPtGenAxis = 1;
static constexpr int kRespCorrRAxis = 2;
static constexpr int kRespJtPtRAxis = 3;
static constexpr int kRespRawRAxis = 4;
// rho rides in the sparse but isn't sliced on here -- every extraction below
// projects the ratio axis, integrating over rho (inclusive-in-rho by default,
// mirroring how eta_reco was carried before it got its own sliced path)
static constexpr int kRespRhoAxis = 5;

static constexpr int kNCollections = 3;
static const char *const kCollectionKeys[kNCollections] = {"incl", "tag",
                                                           "probe"};
static const char *const kCollectionSuffixes[kNCollections] = {
    "_incl_resp", "_tag_resp", "_probe_resp"};

// Three response ratios in one sparse
struct ResponseVariant {
  int axis;
  const char *tag;
  const char *label;
};
static constexpr int kNVariants = 3;
static const ResponseVariant kVariants[kNVariants] = {
    {kRespCorrRAxis, "corr", "p_{T}^{corr}/p_{T}^{gen}"},
    {kRespJtPtRAxis, "reco", "p_{T}^{reco}/p_{T}^{gen}"},
    {kRespRawRAxis, "raw", "p_{T}^{raw}/p_{T}^{gen}"},
};

// Five parallel extraction methods -- Gauss/doubleGauss/trunc90/trunc95/
// crystalball -- matching Step 2's method set (Methods.h), each converting
// a fitted/truncated (mean, sigma) of the response ratio into JES (mean)
// and JER (sigma/mean).

struct ResponseFitResult {
  double jes = 0, jesErr = 0, jer = 0, jerErr = 0;
  bool valid = false;
};

static ResponseFitResult ToJesJer(double mean, double meanErr, double sigma,
                                  double sigmaErr) {
  ResponseFitResult r;
  if (std::abs(mean) < 1e-6) {
    return r;
  }
  r.jes = mean;
  r.jesErr = meanErr;
  r.jer = sigma / mean;
  r.jerErr = r.jer * TMath::Sqrt(TMath::Power(sigmaErr / sigma, 2.0) +
                                 TMath::Power(meanErr / mean, 2.0));
  r.valid = true;
  return r;
}

struct PeakLocation {
  double mean = 1.0, sigma = 0.1;
  bool valid = false;
};

// Pass 1 of the two-pass strategy already proven for this histogram shape in
// runPlotting's response_dist QA guide fit (ResponsePlots.h::FitResponseGuide):
// locate the real peak within [1-halfWidth, 1+halfWidth] before any of the
// Gauss/doubleGauss/crystalball methods below run their own, narrower
// (+-2 sigma) fit -- a single unconstrained pass chases the response axis's
// long non-Gaussian tail toward its hard [0,2] edges. trunc90/trunc95 don't
// need this, they read a fraction of the integral directly.
static PeakLocation LocateResponsePeak(TH1D *h, double halfWidth,
                                       int minEntries) {
  PeakLocation loc;
  if (!CanFit(h, minEntries)) {
    return loc;
  }
  TF1 pass1(Form("_resploc_%s", h->GetName()), "gaus", 1.0 - halfWidth,
           1.0 + halfWidth);
  pass1.SetParameter(0, h->GetMaximum());
  pass1.SetParameter(1, h->GetMean());
  pass1.SetParameter(2, std::max(h->GetRMS(), 1e-3));
  TFitResultPtr res = h->Fit(&pass1, "NQSR");
  if (!res.Get() || !res->IsValid()) {
    return loc;
  }
  loc.mean = res->Parameter(1);
  loc.sigma = std::max(res->Parameter(2), 1e-3);
  loc.valid = true;
  return loc;
}

static ResponseFitResult ExtractResponseGauss(TH1D *h, double halfWidth,
                                              int minEntries) {
  ResponseFitResult r;
  PeakLocation loc = LocateResponsePeak(h, halfWidth, minEntries);
  if (!loc.valid) {
    return r;
  }
  TF1 fit(Form("_respg_%s", h->GetName()), "gaus", loc.mean - 2.0 * loc.sigma,
         loc.mean + 2.0 * loc.sigma);
  fit.SetParameters(h->GetMaximum(), loc.mean, loc.sigma);
  TFitResultPtr res = h->Fit(&fit, "NQSR");
  if (res.Get() && res->IsValid()) {
    r = ToJesJer(res->Parameter(1), res->ParError(1), res->Parameter(2),
                res->ParError(2));
  }
  return r;
}

// mirrors CalibrationExtractor.cxx's FitDoubleGauss weighting scheme
// (amplitude-sigma-weighted mixture mean/sigma), windowed around the
// response ratio's own located peak instead of a fixed +-halfWidth about 0.
static ResponseFitResult ExtractResponseDoubleGauss(TH1D *h, double halfWidth,
                                                     int minEntries) {
  ResponseFitResult r;
  PeakLocation loc = LocateResponsePeak(h, halfWidth, minEntries);
  if (!loc.valid) {
    return r;
  }
  TF1 fit(Form("_respdg_%s", h->GetName()), "gaus(0)+gaus(3)",
         loc.mean - 2.0 * loc.sigma, loc.mean + 2.0 * loc.sigma);
  fit.SetParameter(0, h->GetMaximum());
  fit.SetParameter(1, loc.mean);
  fit.SetParameter(2, 0.5 * loc.sigma);
  fit.SetParameter(3, 0.3 * h->GetMaximum());
  fit.SetParameter(4, loc.mean);
  fit.SetParameter(5, 1.5 * loc.sigma);
  TFitResultPtr res = h->Fit(&fit, "NQSR");
  if (!res.Get() || !res->IsValid()) {
    return r;
  }
  const double amp1 = res->Parameter(0), mean1 = res->Parameter(1),
              sigma1 = res->Parameter(2);
  const double amp2 = res->Parameter(3), mean2 = res->Parameter(4),
              sigma2 = res->Parameter(5);
  const double w1 = std::fabs(amp1 * sigma1);
  const double w2 = std::fabs(amp2 * sigma2);
  const double wsum = w1 + w2;
  if (wsum <= 0) {
    return r;
  }
  const double f1 = w1 / wsum, f2 = w2 / wsum;
  const double mean = f1 * mean1 + f2 * mean2;
  const double meanErr = std::sqrt(TMath::Power(f1 * res->ParError(1), 2) +
                                   TMath::Power(f2 * res->ParError(4), 2));
  const double sigma = std::sqrt(std::max(
      0.0, f1 * (sigma1 * sigma1 + mean1 * mean1) +
               f2 * (sigma2 * sigma2 + mean2 * mean2) - mean * mean));
  const double sigmaErr = std::sqrt(TMath::Power(f1 * res->ParError(2), 2) +
                                    TMath::Power(f2 * res->ParError(5), 2));
  return ToJesJer(mean, meanErr, sigma, sigmaErr);
}

// Central-fraction truncated mean/RMS, same estimator as Step 2's
// trunc90/trunc95 (CalibrationExtractor.cxx), fed the response ratio
// (peaked at 1.0) instead of A (peaked at 0). fraction is hardcoded at the
// call site (0.90 or 0.95), matching Step 2's own hardcoded fractions --
// the whole point of trunc90/trunc95 side by side is a fixed comparison,
// not a tunable single fraction (this used to be config-driven via
// [cuts] response_trunc_fraction; dropped when this became one of five
// parallel methods instead of the sole estimator).
static ResponseFitResult ExtractResponseTrunc(TH1D *h, double fraction,
                                              int minEntries) {
  auto [binLo, binHi] = FindTruncBins(h, fraction, minEntries);
  TruncResult t = TruncMeanInRange(h, binLo, binHi);
  if (!t.valid) {
    return ResponseFitResult{};
  }
  return ToJesJer(t.mean, t.meanErr, t.sigma, t.sigmaErr);
}

// double-sided Crystal Ball (DoubleSidedCB, include/CrystalBall.h), same
// two-stage seeding as CalibrationExtractor.cxx's FitCrystalBall: LocateResponsePeak
// plays the role of Step 2's Gaussian pre-fit (seeds mean/sigma), and the CB
// itself fits the same fixed +-halfWidth window as the pre-fit (centered on
// 1.0, the response ratio's natural center, mirroring Step 2's fixed window
// centered on 0) -- not a narrower +-2 sigma window. A double-sided CB's four
// tail parameters are constrained by the distribution's tails, which a tight
// +-2 sigma window mostly excludes; that under-constraint was previously
// causing MINOS to fail validity on the large majority of bins.
static ResponseFitResult ExtractResponseCrystalBall(TH1D *h, double halfWidth,
                                                     int minEntries) {
  ResponseFitResult r;
  PeakLocation loc = LocateResponsePeak(h, halfWidth, minEntries);
  if (!loc.valid) {
    return r;
  }
  TF1 *cb = new TF1(Form("_respcb_%s", h->GetName()), DoubleSidedCB,
                    1.0 - halfWidth, 1.0 + halfWidth, 7);
  cb->SetParameter(0, h->GetMaximum());
  cb->SetParameter(1, loc.mean);
  cb->SetParameter(2, loc.sigma);
  cb->SetParameter(3, 2.0); // a1
  cb->SetParameter(4, 5.0); // p1
  cb->SetParameter(5, 2.0); // a2
  cb->SetParameter(6, 5.0); // p2
  cb->SetParLimits(3, 0.1, 10.0);
  cb->SetParLimits(4, 0.1, 100.0);
  cb->SetParLimits(5, 0.1, 10.0);
  cb->SetParLimits(6, 0.1, 100.0);
  TFitResultPtr res = h->Fit(cb, "NQSR");
  if (res.Get() && res->IsValid()) {
    r = ToJesJer(res->Parameter(1), res->ParError(1), res->Parameter(2),
                res->ParError(2));
  }
  delete cb;
  return r;
}

static ResponseFitResult ExtractResponseMethod(TH1D *h, int method,
                                                double halfWidth,
                                                int minEntries) {
  switch (method) {
  case 0:
    return ExtractResponseGauss(h, halfWidth, minEntries);
  case 1:
    return ExtractResponseDoubleGauss(h, halfWidth, minEntries);
  case 2:
    return ExtractResponseTrunc(h, 0.90, minEntries);
  case 3:
    return ExtractResponseTrunc(h, 0.95, minEntries);
  case 4:
    return ExtractResponseCrystalBall(h, halfWidth, minEntries);
  default:
    return ResponseFitResult{};
  }
}

// vs pt_gen -- one JES/JER histogram pair per method (5 parallel estimators,
// see the block above), all fit against the same per-bin projection.
static void ExtractVsPtGen(THnSparse *h, const TString &cone,
                           const TString &collection,
                           const ResponseVariant &variant,
                           const AxisBins &ptBins, TDirectory *dQA,
                           TDirectory *dOut, double halfWidth,
                           int minEntries) {

  const int nPt = ptBins.nBins;
  const double lo = ptBins.lo;
  const double hi = ptBins.hi;
  const double width = (hi - lo) / nPt;

  TH1D *hJES[kNMethods];
  TH1D *hJER[kNMethods];
  for (int m = 0; m < kNMethods; m++) {
    TString jesName =
        L2Name::ObjectName(cone, "JES", {variant.tag, "vs_ptgen"},
                           {collection, kMethodKeys[m]});
    TString jerName =
        L2Name::ObjectName(cone, "JER", {variant.tag, "vs_ptgen"},
                           {collection, kMethodKeys[m]});
    hJES[m] = new TH1D(jesName, "", nPt, lo, hi);
    hJER[m] = new TH1D(jerName, "", nPt, lo, hi);
    hJES[m]->GetXaxis()->SetTitle("p_{T}^{gen} [GeV]");
    hJES[m]->GetYaxis()->SetTitle(Form("JES = #LT %s #GT", variant.label));
    hJER[m]->GetXaxis()->SetTitle(hJES[m]->GetXaxis()->GetTitle());
    hJER[m]->GetYaxis()->SetTitle(
        Form("JER = #sigma / #LT %s #GT", variant.label));
  }

  for (int ip = 0; ip < nPt; ip++) {
    const double ptLo = lo + ip * width;
    const double ptHi = ptLo + width;

    TString distName =
        L2Name::ObjectName(cone, TString("response_") + variant.tag,
                           {L2Name::PtGenKey(ptLo, ptHi)}, {collection});

    TH1D *hProj =
        ProjectTHnSparse1D(h, variant.axis, {{kRespPtGenAxis, ptLo, ptHi}});
    hProj->SetName(distName);

    dQA->cd();
    hProj->Write();

    for (int m = 0; m < kNMethods; m++) {
      ResponseFitResult fr =
          ExtractResponseMethod(hProj, m, halfWidth, minEntries);
      if (fr.valid) {
        hJES[m]->SetBinContent(ip + 1, fr.jes);
        hJES[m]->SetBinError(ip + 1, fr.jesErr);
        hJER[m]->SetBinContent(ip + 1, fr.jer);
        hJER[m]->SetBinError(ip + 1, fr.jerErr);
      }
    }
    delete hProj;
  }

  dOut->cd();
  for (int m = 0; m < kNMethods; m++) {
    hJES[m]->Write();
    hJER[m]->Write();
    delete hJES[m];
    delete hJER[m];
  }
}

// Adjacent-pair (disjoint) bins from a flat sorted edge array, e.g.
// kAbsEtaEdges/kEtaRangeEdges/kEtaEdges -- the shape ExtractPerEtaVsPtGen
// always took before it needed to also support cumulative (overlapping,
// all-from-zero) ranges, which can't be expressed as one flat edge array.
static std::vector<std::pair<double, double>>
EdgesToPairs(const std::vector<Double_t> &edges) {
  std::vector<std::pair<double, double>> pairs;
  for (int i = 0; i + 1 < (int)edges.size(); i++) {
    pairs.push_back({edges[i], edges[i + 1]});
  }
  return pairs;
}

// RangeBin's lo/hi as (lo,hi) pairs -- BuildEtaRangeSlicesCumulative's
// ranges overlap (all start at 0), so they can't round-trip through
// EdgesToPairs' flat-edge-array assumption.
static std::vector<std::pair<double, double>>
EtaRangeSlicesToPairs(const std::vector<RangeBin> &slices) {
  std::vector<std::pair<double, double>> pairs;
  for (const auto &sl : slices) {
    pairs.push_back({sl.lo, sl.hi});
  }
  return pairs;
}

// JES/JER vs pT_gen for each (lo,hi) eta range in etaRanges. The caller
// supplies the scheme and the matching sparse: kAbsEtaEdges/kEtaRangeEdges
// (via EdgesToPairs) on a folded |eta| sparse, kEtaEdges on the raw signed
// one, or an explicit cumulative pair list (BuildEtaRangeSlicesCumulative).
// Feeds the pT-resolution text writer (fine schemes) and the detector-region
// overlay plots (coarse disjoint + cumulative schemes).
//
// restricts axis 0 to each eta range in turn, projects and extracts per pT_gen bin
// writes {cone}_JES_{variant}_vs_ptgen_{etaKey}_{collection} to dOut
static void
ExtractPerEtaVsPtGen(THnSparse *h, const TString &cone,
                     const TString &collection, const ResponseVariant &variant,
                     const AxisBins &ptBins, TDirectory *dOut, double halfWidth,
                     int minEntries,
                     const std::vector<std::pair<double, double>> &etaRanges) {

  const int nPt = ptBins.nBins;
  const double lo = ptBins.lo, hi = ptBins.hi;
  const double width = (hi - lo) / nPt;

  for (const auto &etaRange : etaRanges) {
    const TString etaKey = L2Name::EtaKey(etaRange.first, etaRange.second);

    h->GetAxis(kRespEtaRecoAxis)->SetRangeUser(etaRange.first, etaRange.second);

    TH1D *hJES[kNMethods];
    TH1D *hJER[kNMethods];
    for (int m = 0; m < kNMethods; m++) {
      TString jesName =
          L2Name::ObjectName(cone, "JES", {variant.tag, "vs_ptgen", etaKey},
                             {collection, kMethodKeys[m]});
      TString jerName =
          L2Name::ObjectName(cone, "JER", {variant.tag, "vs_ptgen", etaKey},
                             {collection, kMethodKeys[m]});
      hJES[m] = new TH1D(jesName, "", nPt, lo, hi);
      hJER[m] = new TH1D(jerName, "", nPt, lo, hi);
      hJES[m]->GetXaxis()->SetTitle("p_{T}^{gen} [GeV]");
      hJES[m]->GetYaxis()->SetTitle(Form("JES = #LT %s #GT", variant.label));
      hJER[m]->GetXaxis()->SetTitle(hJES[m]->GetXaxis()->GetTitle());
      hJER[m]->GetYaxis()->SetTitle(
          Form("JER = #sigma / #LT %s #GT", variant.label));
    }

    for (int ip = 0; ip < nPt; ip++) {
      const double ptLo = lo + ip * width;
      const double ptHi = ptLo + width;

      h->GetAxis(kRespPtGenAxis)->SetRangeUser(ptLo, ptHi);
      TH1D *hProj = ProjectTHnSparse1D(h, variant.axis, {});
      h->GetAxis(kRespPtGenAxis)->SetRange(0, 0);

      for (int m = 0; m < kNMethods; m++) {
        ResponseFitResult fr =
            ExtractResponseMethod(hProj, m, halfWidth, minEntries);
        if (fr.valid) {
          hJES[m]->SetBinContent(ip + 1, fr.jes);
          hJES[m]->SetBinError(ip + 1, fr.jesErr);
          hJER[m]->SetBinContent(ip + 1, fr.jer);
          hJER[m]->SetBinError(ip + 1, fr.jerErr);
        }
      }
      delete hProj;
    }

    h->GetAxis(kRespEtaRecoAxis)->SetRange(0, 0);

    dOut->cd();
    for (int m = 0; m < kNMethods; m++) {
      hJES[m]->Write();
      hJER[m]->Write();
      delete hJES[m];
      delete hJER[m];
    }
  }
}

// main
void runResponse(TString inputFile, TString outputFile) {

  const AnalysisConfig &cfg = Config();
  PrintConfigSummary(cfg);

  const int minEntries =
      (cfg.minEntriesPerBin > 0) ? cfg.minEntriesPerBin : 100;
  const double halfWidth = cfg.responseGausFitHalfWidth;

  TFile *fIn = TFile::Open(inputFile, "read");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "Cannot open " << inputFile << "\n";
    return;
  }

  // incl/tag/probe response sparses, ConeHistograms instantiation
  int totalSteps = 0;
  for (const TString &cone : cfg.coneLabels) {
    if (fIn->Get(cone + "/" + cone + "_incl_resp")) {
      totalSteps += kNCollections * kNVariants;
    }
  }
  if (totalSteps == 0) {
    std::cerr << "No MC response histograms found in " << inputFile
              << " -- are these MC asymmetries? (they need to be)\n";
    fIn->Close();
    return;
  }

  {
    Ssiz_t sl = outputFile.Last('/');
    if (sl != kNPOS) {
      gSystem->mkdir(TString(outputFile(0, sl)), kTRUE);
    }
  }
  TFile *fOut = new TFile(outputFile, "recreate");

  const bool wantAbsEta = (cfg.etaModeOutput != "eta");
  const bool wantFullEta = (cfg.etaModeOutput != "abseta");

  BinningConfig bins;
  ProgressBar pb("Extracting response:", totalSteps);

  for (const TString &cone : cfg.coneLabels) {

    TDirectory *coneDir = (TDirectory *)fIn->Get(cone);
    THnSparse *hSentinel =
        coneDir ? (THnSparse *)coneDir->Get(cone + "_incl_resp") : nullptr;
    if (!hSentinel) {
      continue;
    }

    TDirectory *coneDirOut = fOut->mkdir(cone.Data());
    TDirectory *dQA_ptgen = coneDirOut->mkdir("QA_response_ptgen");
    TDirectory *dPerAbsEta =
        wantAbsEta ? coneDirOut->mkdir("JER_per_abseta") : nullptr;
    TDirectory *dPerEta =
        wantFullEta ? coneDirOut->mkdir("JER_per_eta") : nullptr;
    // coarse |eta| detector-region slices for the JES/JER overlay plot
    TDirectory *dPerEtaRange =
        wantAbsEta ? coneDirOut->mkdir("JER_per_etarange") : nullptr;
    // cumulative (nested-from-zero) variant of the same coarse regions
    TDirectory *dPerEtaRangeCumulative =
        wantAbsEta ? coneDirOut->mkdir("JER_per_etarange_cumulative")
                  : nullptr;

    for (int ic = 0; ic < kNCollections; ic++) {
      const TString collection = kCollectionKeys[ic];
      THnSparse *hRaw =
          (THnSparse *)coneDir->Get(cone + kCollectionSuffixes[ic]);
      if (!hRaw) {
        std::cerr << "Missing " << cone << kCollectionSuffixes[ic] << " in "
                  << inputFile << "\n";
        continue;
      }

      for (int iv = 0; iv < kNVariants; iv++) {
        ExtractVsPtGen(hRaw, cone, collection, kVariants[iv], bins.pt,
                       dQA_ptgen, coneDirOut, halfWidth, minEntries);
        pb.Update();
      }

      // per |eta| corr-variant extraction, folded once: the fine kAbsEtaEdges
      // scheme feeds the CMS JER pT-resolution text files, the coarse
      // kEtaRangeEdges scheme feeds the detector-region JES/JER overlay plot
      if (wantAbsEta) {
        THnSparse *hFolded = FoldEtaAxis(
            hRaw, kRespEtaRecoAxis, cone + kCollectionSuffixes[ic] + "_abseta");
        ExtractPerEtaVsPtGen(hFolded, cone, collection, kVariants[0], bins.pt,
                             dPerAbsEta, halfWidth, minEntries,
                             EdgesToPairs(kAbsEtaEdges));
        ExtractPerEtaVsPtGen(hFolded, cone, collection, kVariants[0], bins.pt,
                             dPerEtaRange, halfWidth, minEntries,
                             EdgesToPairs(kEtaRangeEdges));
        ExtractPerEtaVsPtGen(hFolded, cone, collection, kVariants[0], bins.pt,
                             dPerEtaRangeCumulative, halfWidth, minEntries,
                             EtaRangeSlicesToPairs(BuildEtaRangeSlicesCumulative()));
        delete hFolded;
      }
      // per (signed) eta_reco extraction for the full-eta pT-resolution files
      if (wantFullEta) {
        ExtractPerEtaVsPtGen(hRaw, cone, collection, kVariants[0], bins.pt,
                             dPerEta, halfWidth, minEntries,
                             EdgesToPairs(kEtaEdges));
      }
    }
  }

  pb.Finish();
  fOut->Close();
  fIn->Close();
}
