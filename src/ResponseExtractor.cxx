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
#include "Naming.h"
#include "Utilities.h"
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

// Gaussian fit around 1.0 -> JES (mean) and JER (sigma/mean)
struct ResponseFitResult {
  double jes = 0, jesErr = 0, jer = 0, jerErr = 0;
  bool valid = false;
};

static ResponseFitResult FitResponse(TH1D *h, double halfWidth,
                                     int minEntries) {
  ResponseFitResult r;
  if (!CanFit(h, minEntries))
    return r;

  TF1 *g = new TF1(Form("_rf_%s", h->GetName()), "gaus", 1.0 - halfWidth,
                   1.0 + halfWidth);
  g->SetParameter(0, h->GetMaximum());
  g->SetParameter(1, h->GetMean());
  g->SetParameter(2, std::max(h->GetRMS(), 1e-3));

  TFitResultPtr res = h->Fit(g, "NQSR");
  if (res.Get() && res->IsValid()) {
    const double mean = res->Parameter(1);
    const double meanErr = res->ParError(1);
    const double sigma = res->Parameter(2);
    const double sigmaErr = res->ParError(2);
    r.jes = mean;
    r.jesErr = meanErr;
    if (std::abs(mean) > 1e-6) {
      r.jer = sigma / mean;
      r.jerErr = r.jer * TMath::Sqrt(TMath::Power(sigmaErr / sigma, 2.0) +
                                     TMath::Power(meanErr / mean, 2.0));
      r.valid = true;
    }
  }
  delete g;
  return r;
}

// vs pt_gen
static void ExtractVsPtGen(THnSparse *h, const TString &cone,
                           const TString &collection,
                           const ResponseVariant &variant,
                           const AxisBins &ptBins, TDirectory *dQA,
                           TDirectory *dOut, double halfWidth, int minEntries) {

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

    ResponseFitResult fr = FitResponse(hProj, halfWidth, minEntries);
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

// per |eta| bin JER and JES vs pT_gen
// used by runTextFilePtResolution
//
// |eta| bin (restrict axis 0), projects and fit per pT_gen bin
// writes {cone}_JES_{variant}_vs_ptgen_{etaKey}_{collection} to dOut

static void ExtractPerAbsEtaVsPtGen(THnSparse *hFolded, const TString &cone,
                                    const TString &collection,
                                    const ResponseVariant &variant,
                                    const AxisBins &ptBins, TDirectory *dOut,
                                    double halfWidth, int minEntries) {

  const int nEta = (int)kAbsEtaEdges.size() - 1;
  const int nPt = ptBins.nBins;
  const double lo = ptBins.lo, hi = ptBins.hi;
  const double width = (hi - lo) / nPt;

  for (int ieta = 0; ieta < nEta; ieta++) {
    const TString etaKey = L2Name::EtaKey(ieta, false);

    hFolded->GetAxis(kRespEtaRecoAxis)
        ->SetRangeUser(kAbsEtaEdges[ieta], kAbsEtaEdges[ieta + 1]);

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

      hFolded->GetAxis(kRespPtGenAxis)->SetRangeUser(ptLo, ptHi);
      TH1D *hProj = ProjectTHnSparse1D(hFolded, variant.axis, {});
      hFolded->GetAxis(kRespPtGenAxis)->SetRange(0, 0);

      ResponseFitResult fr = FitResponse(hProj, halfWidth, minEntries);
      if (fr.valid) {
        hJES->SetBinContent(ip + 1, fr.jes);
        hJES->SetBinError(ip + 1, fr.jesErr);
        hJER->SetBinContent(ip + 1, fr.jer);
        hJER->SetBinError(ip + 1, fr.jerErr);
      }
      delete hProj;
    }

    hFolded->GetAxis(kRespEtaRecoAxis)->SetRange(0, 0);

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

  const double halfWidth = cfg.responseGausFitHalfWidth;
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
    if (fIn->Get(cone + "/" + cone + "_incl_resp"))
      totalSteps += kNCollections * kNVariants;
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

  BinningConfig bins;
  ProgressBar pb("Extracting response:", totalSteps);

  for (const TString &cone : cfg.coneLabels) {

    TDirectory *coneDir = (TDirectory *)fIn->Get(cone);
    THnSparse *hSentinel =
        coneDir ? (THnSparse *)coneDir->Get(cone + "_incl_resp") : nullptr;
    if (!hSentinel)
      continue;

    TDirectory *coneDirOut = fOut->mkdir(cone.Data());
    TDirectory *dQA_ptgen = coneDirOut->mkdir("QA_response_ptgen");
    TDirectory *dPerEta = coneDirOut->mkdir("JER_per_abseta");

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

      // per |eta| bin extraction for corr variant
      // used to write CMS JER pT resolution text files
      THnSparse *hFolded = FoldEtaAxis(
          hRaw, kRespEtaRecoAxis, cone + kCollectionSuffixes[ic] + "_abseta");
      ExtractPerAbsEtaVsPtGen(hFolded, cone, collection, kVariants[0], bins.pt,
                              dPerEta, halfWidth, minEntries);
      delete hFolded;
    }
  }

  pb.Finish();
  fOut->Close();
  fIn->Close();
}
