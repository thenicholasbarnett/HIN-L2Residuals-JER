#include "ResponseExtractor.h"

#include "TFile.h"
#include "TH1D.h"
#include "THnSparse.h"
#include "TDirectory.h"
#include "TString.h"
#include "TMath.h"

#include "Binning.h"
#include "Naming.h"
#include "Utilities.h"
#include "Truncation.h"
#include "ProgressBar.h"
#include "AnalysisConfig.h"

#include <vector>
#include <iostream>
#include <cmath>

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

// central-95% truncated RMS -> JES (mean) and JER (sigma/mean). Same trunc95
// estimator as Step 2's A-distribution extraction (CalibrationExtractor's
// FindTruncBins/TruncMeanInRange), here on the response ratio peaked at 1.0
// rather than A peaked at 0 -- trims the long non-Gaussian tail toward the
// response axis's [0,2] edges without a fit window.
static constexpr double kResponseTruncFraction = 0.95;

struct ResponseFitResult {
  double jes = 0, jesErr = 0, jer = 0, jerErr = 0;
  bool valid = false;
};

static ResponseFitResult ExtractResponse(TH1D *h, int minEntries) {
  ResponseFitResult r;
  auto [binLo, binHi] = FindTruncBins(h, kResponseTruncFraction, minEntries);
  TruncResult t = TruncMeanInRange(h, binLo, binHi);
  if (!t.valid || std::abs(t.mean) < 1e-6) {
    return r;
  }
  r.jes = t.mean;
  r.jesErr = t.meanErr;
  r.jer = t.sigma / t.mean;
  r.jerErr = r.jer * TMath::Sqrt(TMath::Power(t.sigmaErr / t.sigma, 2.0) +
                                 TMath::Power(t.meanErr / t.mean, 2.0));
  r.valid = true;
  return r;
}

// vs pt_gen
static void ExtractVsPtGen(THnSparse *h, const TString &cone,
                           const TString &collection,
                           const ResponseVariant &variant,
                           const AxisBins &ptBins, TDirectory *dQA,
                           TDirectory *dOut, int minEntries) {

  const int nPt = ptBins.nBins;
  const double lo = ptBins.lo;
  const double hi = ptBins.hi;
  const double width = (hi - lo) / nPt;

  TString jesName =
      L2Name::ObjectName(cone, "JES", {variant.tag, "vs_ptgen"}, {collection});
  TString jerName =
      L2Name::ObjectName(cone, "JER", {variant.tag, "vs_ptgen"}, {collection});
  TH1D *hJES = new TH1D(jesName, "", nPt, lo, hi);
  TH1D *hJER = new TH1D(jerName, "", nPt, lo, hi);
  hJES->GetXaxis()->SetTitle("p_{T}^{gen} [GeV]");
  hJES->GetYaxis()->SetTitle(Form("JES = #LT %s #GT", variant.label));
  hJER->GetXaxis()->SetTitle(hJES->GetXaxis()->GetTitle());
  hJER->GetYaxis()->SetTitle(Form("JER = #sigma / #LT %s #GT", variant.label));

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

    ResponseFitResult fr = ExtractResponse(hProj, minEntries);
    if (fr.valid) {
      hJES->SetBinContent(ip + 1, fr.jes);
      hJES->SetBinError(ip + 1, fr.jesErr);
      hJER->SetBinContent(ip + 1, fr.jer);
      hJER->SetBinError(ip + 1, fr.jerErr);
    }
    delete hProj;
  }

  dOut->cd();
  hJES->Write();
  hJER->Write();
  delete hJES;
  delete hJER;
}

// per eta bin (|eta| if fullEta=false, signed eta_reco if fullEta=true) JER
// and JES vs pT_gen, used by runTextFilePtResolution
//
// restricts axis 0 to each eta bin in turn, projects and extracts per pT_gen bin
// writes {cone}_JES_{variant}_vs_ptgen_{etaKey}_{collection} to dOut
static void ExtractPerEtaVsPtGen(THnSparse *h, const TString &cone,
                                 const TString &collection,
                                 const ResponseVariant &variant,
                                 const AxisBins &ptBins, TDirectory *dOut,
                                 int minEntries, bool fullEta) {

  const std::vector<Double_t> &etaEdges = fullEta ? kEtaEdges : kAbsEtaEdges;
  const int nEta = (int)etaEdges.size() - 1;
  const int nPt = ptBins.nBins;
  const double lo = ptBins.lo, hi = ptBins.hi;
  const double width = (hi - lo) / nPt;

  for (int ieta = 0; ieta < nEta; ieta++) {
    const TString etaKey = L2Name::EtaKey(ieta, fullEta);

    h->GetAxis(kRespEtaRecoAxis)
        ->SetRangeUser(etaEdges[ieta], etaEdges[ieta + 1]);

    TString jesName = L2Name::ObjectName(
        cone, "JES", {variant.tag, "vs_ptgen", etaKey}, {collection});
    TString jerName = L2Name::ObjectName(
        cone, "JER", {variant.tag, "vs_ptgen", etaKey}, {collection});
    TH1D *hJES = new TH1D(jesName, "", nPt, lo, hi);
    TH1D *hJER = new TH1D(jerName, "", nPt, lo, hi);
    hJES->GetXaxis()->SetTitle("p_{T}^{gen} [GeV]");
    hJES->GetYaxis()->SetTitle(Form("JES = #LT %s #GT", variant.label));
    hJER->GetXaxis()->SetTitle(hJES->GetXaxis()->GetTitle());
    hJER->GetYaxis()->SetTitle(
        Form("JER = #sigma / #LT %s #GT", variant.label));

    for (int ip = 0; ip < nPt; ip++) {
      const double ptLo = lo + ip * width;
      const double ptHi = ptLo + width;

      h->GetAxis(kRespPtGenAxis)->SetRangeUser(ptLo, ptHi);
      TH1D *hProj = ProjectTHnSparse1D(h, variant.axis, {});
      h->GetAxis(kRespPtGenAxis)->SetRange(0, 0);

      ResponseFitResult fr = ExtractResponse(hProj, minEntries);
      if (fr.valid) {
        hJES->SetBinContent(ip + 1, fr.jes);
        hJES->SetBinError(ip + 1, fr.jesErr);
        hJER->SetBinContent(ip + 1, fr.jer);
        hJER->SetBinError(ip + 1, fr.jerErr);
      }
      delete hProj;
    }

    h->GetAxis(kRespEtaRecoAxis)->SetRange(0, 0);

    dOut->cd();
    hJES->Write();
    hJER->Write();
    delete hJES;
    delete hJER;
  }
}

// main
void runResponse(TString inputFile, TString outputFile) {

  const AnalysisConfig &cfg = Config();
  PrintConfigSummary(cfg);

  const int minEntries =
      (cfg.minEntriesPerBin > 0) ? cfg.minEntriesPerBin : 100;

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
                       dQA_ptgen, coneDirOut, minEntries);
        pb.Update();
      }

      // per |eta| and per (signed) eta_reco extraction for corr variant,
      // used to write CMS JER pT resolution text files
      if (wantAbsEta) {
        THnSparse *hFolded = FoldEtaAxis(
            hRaw, kRespEtaRecoAxis, cone + kCollectionSuffixes[ic] + "_abseta");
        ExtractPerEtaVsPtGen(hFolded, cone, collection, kVariants[0], bins.pt,
                             dPerAbsEta, minEntries, false);
        delete hFolded;
      }
      if (wantFullEta) {
        ExtractPerEtaVsPtGen(hRaw, cone, collection, kVariants[0], bins.pt,
                             dPerEta, minEntries, true);
      }
    }
  }

  pb.Finish();
  fOut->Close();
  fIn->Close();
}
