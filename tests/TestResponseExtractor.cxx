#include <iostream>
#include <cmath>
#include <cstdio>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"

#include "Binning.h"
#include "DijetHistograms.h"
#include "ResponseExtractor.h"

static int nPass = 0;
static int nFail = 0;
static constexpr double kEps = 1e-6;

void Check(bool cond, const char *msg) {
  if (cond) {
    std::cout << "  PASS  " << msg << std::endl;
    nPass++;
  } else {
    std::cout << "  FAIL  " << msg << std::endl;
    nFail++;
  }
}

// Response values symmetric about 1.0
static const double kResponseValues[] = {0.90, 0.93, 0.96, 0.99, 1.00,
                                         1.01, 1.04, 1.07, 1.10};
static constexpr int kNResponseValues = 9;
static constexpr int kRepeatsPerValue =
    20; // 180 total entries, well above min_entries_per_bin=100

int main() {
  gROOT->SetBatch(true);

  std::cout << "=== TestResponseExtractor ===" << std::endl;

  const char *inputPath = "/tmp/tre_input.root";
  const char *outputPath = "/tmp/tre_output.root";
  std::remove(inputPath);
  std::remove(outputPath);

  // synthetic MC asymmetry file (ak4PF, incl only)
  // corrPt, rawPt, and jtPt use same symmetric responses
  // JES(corr)=JES(reco)=JES(raw)~1.0 by construction -- 
  // this test extraction of all three variants, entry-count guard
  {
    BinningConfig bins;
    ConeHistograms h;
    h.Init("ak4PF", bins, true);

    // eta=0.13 -> eta_reco bin covering [0,0.261]
    // refPt=105 -> pt_gen bin [100,110)
    const float eta = 0.13f;
    const float refPt = 105.0f;
    for (int iv = 0; iv < kNResponseValues; iv++) {
      const float pt = (float)(kResponseValues[iv] * refPt);
      for (int r = 0; r < kRepeatsPerValue; r++) {
        h.FillInclJetResp(pt, pt, pt, eta, refPt, 1.0f);
      }
    }

    // deliberately under-filled pt_gen bin, only 5 entries
    // should NOT produce a valid fit.
    for (int r = 0; r < 5; r++) {
      h.FillInclJetResp(300.0f, 300.0f, 300.0f, eta, 300.0f, 1.0f);
    }

    TFile *fo = new TFile(inputPath, "recreate");
    TDirectory *dir = fo->mkdir("ak4PF");
    h.Write(dir);
    fo->Close();
    delete fo;
  }

  runResponse(inputPath, outputPath);

  TFile *fIn = TFile::Open(outputPath, "read");
  Check(fIn && !fIn->IsZombie(), "output file opens");
  if (!fIn || fIn->IsZombie()) {
    std::cout << "\n=== " << nPass << " passed, " << nFail
              << " failed ===" << std::endl;
    return 1;
  }

  std::cout << "\n[1] vs pt_gen, all three variants" << std::endl;
  {
    static const char *const kVariantTags[3] = {"corr", "reco", "raw"};
    for (const char *tag : kVariantTags) {
      TH1D *hJES =
          (TH1D *)fIn->Get(Form("ak4PF/ak4PF_JES_%s_vs_ptgen_incl", tag));
      TH1D *hJER =
          (TH1D *)fIn->Get(Form("ak4PF/ak4PF_JER_%s_vs_ptgen_incl", tag));
      Check(hJES != nullptr, Form("JES_%s_vs_ptgen_incl exists", tag));
      Check(hJER != nullptr, Form("JER_%s_vs_ptgen_incl exists", tag));
      if (hJES && hJER) {
        const int bin =
            hJES->GetXaxis()->FindBin(105.0); // pt_gen=105 -> [100,110) bin
        const double jes = hJES->GetBinContent(bin);
        const double jer = hJER->GetBinContent(bin);
        Check(std::fabs(jes - 1.0) < 0.02,
              Form("JES_%s(pt_gen~105) close to 1.0 (symmetric synthetic "
                   "response)",
                   tag));
        Check(jer > 0.01 && jer < 0.20,
              Form("JER_%s(pt_gen~105) in a sane positive range", tag));

        const int emptyBin =
            hJES->GetXaxis()->FindBin(305.0); // under-filled bin
        Check(std::fabs(hJES->GetBinContent(emptyBin)) < kEps,
              Form("under-filled pt_gen bin (5 entries < min_entries_per_bin) "
                   "stays at 0, no fit (%s)",
                   tag));
      }
    }
  }

  std::cout << "\n[2] tag/probe collections were never filled -> objects exist "
               "but bins stay zero"
            << std::endl;
  {
    // ConeHistograms::Init always constructs three response sparses
    // extraction writes JES/JER TH1Ds for unfilled hists still (see fallback to unity)
    TH1D *hJESTag = (TH1D *)fIn->Get("ak4PF/ak4PF_JES_corr_vs_ptgen_tag");
    Check(hJESTag != nullptr,
          "JES_corr_vs_ptgen_tag object still written (unfilled, not absent)");
    if (hJESTag) {
      bool allZero = true;
      for (int b = 1; b <= hJESTag->GetNbinsX(); b++) {
        if (std::fabs(hJESTag->GetBinContent(b)) > kEps) {
          allZero = false;
        }
      }
      Check(allZero, "tag collection (never filled) has no valid fits anywhere "
                     "-> all bins zero");
    }
  }

  std::cout << "\n[3] JER_per_abseta and JER_per_eta (default etaModeOutput "
               "= \"both\")"
            << std::endl;
  {
    // eta=0.13 lands in |eta| bin 0 ([0,0.261])
    // signed full-eta bin 18 ([0,0.261]) 
    // both should carry the real fit
    TH1D *hAbs = (TH1D *)fIn->Get(
        "ak4PF/JER_per_abseta/ak4PF_JER_corr_vs_ptgen_eta_0p0_0p261_incl");
    TH1D *hFull = (TH1D *)fIn->Get(
        "ak4PF/JER_per_eta/ak4PF_JER_corr_vs_ptgen_eta_0p0_0p261_incl");
    Check(hAbs != nullptr, "JER_per_abseta eta_0p0_0p261 exists");
    Check(hFull != nullptr, "JER_per_eta eta_0p0_0p261 exists");
    if (hAbs && hFull) {
      const int bin = hAbs->GetXaxis()->FindBin(105.0);
      const double jerAbs = hAbs->GetBinContent(bin);
      const double jerFull = hFull->GetBinContent(bin);
      Check(jerAbs > 0.01 && jerAbs < 0.20,
            "JER_per_abseta(pt_gen~105) in a sane positive range");
      Check(std::fabs(jerFull - jerAbs) < kEps,
            "JER_per_eta matches JER_per_abseta for the same physical eta bin");
    }

    // the mirror-image negative-eta bin no real entries 
    // should stay all zero, proving JER_per_eta isn't JER_per_abseta mirrored
    TH1D *hFullNeg = (TH1D *)fIn->Get(
        "ak4PF/JER_per_eta/ak4PF_JER_corr_vs_ptgen_eta_m0p261_0p0_incl");
    Check(hFullNeg != nullptr,
          "JER_per_eta eta_m0p261_0p0 (negative side) exists");
    if (hFullNeg) {
      bool allZero = true;
      for (int b = 1; b <= hFullNeg->GetNbinsX(); b++) {
        if (std::fabs(hFullNeg->GetBinContent(b)) > kEps) {
          allZero = false;
        }
      }
      Check(allZero, "JER_per_eta negative-eta bin has no fits -> all zero "
                     "(not mirrored from the positive side)");
    }
  }

  fIn->Close();
  std::remove(inputPath);
  std::remove(outputPath);

  std::cout << "\n=== " << nPass << " passed, " << nFail
            << " failed ===" << std::endl;
  return nFail > 0 ? 1 : 0;
}
