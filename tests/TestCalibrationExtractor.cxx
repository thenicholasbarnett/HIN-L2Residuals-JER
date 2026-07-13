#include <iostream>
#include <cmath>
#include <cstdio>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"

#include "Binning.h"
#include "Dijet.h"
#include "DijetHistograms.h"
#include "CalibrationExtractor.h"

static int nPass = 0;
static int nFail = 0;

void Check(bool cond, const char *msg) {
  if (cond) {
    std::cout << "  PASS  " << msg << std::endl;
    nPass++;
  } else {
    std::cout << "  FAIL  " << msg << std::endl;
    nFail++;
  }
}

// Deliberately right-skewed synthetic A distribution -- a Gaussian-ish core
// plus a longer, heavier positive tail than negative, a shape Gauss/
// double-Gauss can't represent (both symmetric or near-symmetric by
// construction) but Crystal Ball can, via independent per-side tail
// parameters (a1/p1 low side, a2/p2 high side). This is the point of this
// test: not just "does FitCrystalBall run," but "does it converge on data
// that actually needs its asymmetric tails."
struct AVal {
  double a;
  int count;
};
static const AVal kAValues[] = {
    {-0.10, 5},  {-0.07, 10}, {-0.05, 20}, {-0.03, 40}, {-0.01, 60},
    {0.00, 80},  {0.01, 60},  {0.03, 45},  {0.05, 35},  {0.07, 28},
    {0.09, 22},  {0.11, 18},  {0.13, 14},  {0.15, 11},  {0.18, 8},
    {0.21, 6},   {0.25, 4},   {0.30, 2},
};
static constexpr int kNAValues = sizeof(kAValues) / sizeof(kAValues[0]);

// fills the same asymmetric A distribution at a single low alpha (0.02),
// included in every cumulative alpha-threshold slice up to 0.31 (the fit
// range) -- enough distinct alpha points for FitAndExtrapolate's linear fit
// (needs >= 2), without needing a physically varying alpha dependence.
static void FillAsymmetricTail(ConeHistograms &h, float ptavg) {
  float pt[2] = {100.0f, 100.0f};
  float eta[2] = {0.13f, 0.13f}; // lands in |eta| bin [0, 0.261)
  float phi[2] = {0.0f, 1.0f};

  DijetResult d;
  d.valid = true;
  d.tagIdx = 0;
  d.probeIdx = 1;
  d.ptavg = ptavg;
  d.alpha = 0.02f;

  for (int iv = 0; iv < kNAValues; iv++) {
    d.A = (float)kAValues[iv].a;
    for (int r = 0; r < kAValues[iv].count; r++) {
      h.Fill(d, pt, eta, phi, 1.0f);
    }
  }
}

int main() {
  gROOT->SetBatch(true);

  std::cout << "=== TestCalibrationExtractor ===" << std::endl;

  const char *dataPath = "/tmp/tce_data.root";
  const char *mcPath = "/tmp/tce_mc.root";
  const char *outputPath = "/tmp/tce_output.root";
  std::remove(dataPath);
  std::remove(mcPath);
  std::remove(outputPath);

  // ptavg=105 -> pT_avg slice [100,175) -- a real slice from cfg's
  // ptavg_edges. Same asymmetric shape in both data and MC: this test only
  // checks the crystalball fit itself converges to a sane value, not the
  // physical accuracy of a data/MC ratio.
  {
    BinningConfig bins;
    ConeHistograms h;
    h.Init("ak4PF", bins, false);
    FillAsymmetricTail(h, 105.0f);

    TFile *fo = new TFile(dataPath, "recreate");
    TDirectory *dir = fo->mkdir("ak4PF");
    h.Write(dir);
    fo->Close();
    delete fo;
  }
  {
    BinningConfig bins;
    ConeHistograms h;
    h.Init("ak4PF", bins, false);
    FillAsymmetricTail(h, 105.0f);

    TFile *fo = new TFile(mcPath, "recreate");
    TDirectory *dir = fo->mkdir("ak4PF");
    h.Write(dir);
    fo->Close();
    delete fo;
  }

  std::cout << "\n[1] JEC mode: crystalball intercept exists and converges"
            << std::endl;
  {
    runCalibration(dataPath, mcPath, outputPath, CalibrationMode::JEC);

    TFile *fIn = TFile::Open(outputPath, "read");
    Check(fIn && !fIn->IsZombie(), "output file opens (JEC)");
    if (fIn && !fIn->IsZombie()) {
      TH1D *hGauss = (TH1D *)fIn->Get(
          "ak4PF/ak4PF_intercept_abseta_ptavg_100_175_gauss");
      TH1D *hCB = (TH1D *)fIn->Get(
          "ak4PF/ak4PF_intercept_abseta_ptavg_100_175_crystalball");
      Check(hGauss != nullptr, "gauss intercept object exists");
      Check(hCB != nullptr, "crystalball intercept object exists");
      if (hGauss && hCB) {
        // eta bin 0 ([0, 0.261)) is where the synthetic jets were filled
        const double gaussVal = hGauss->GetBinContent(1);
        const double cbVal = hCB->GetBinContent(1);
        Check(std::fabs(gaussVal) > 1e-9,
              "gauss intercept converged to a non-zero value");
        Check(std::fabs(cbVal) > 1e-9,
              "crystalball intercept converged to a non-zero value");
        Check(std::isfinite(cbVal), "crystalball intercept is finite");
        // same underlying data/MC A distributions -> R should extrapolate
        // close to 1 for both methods, and the two methods' fitted means
        // (both reading the same right-skewed distribution) should agree
        // to within a sane tolerance, not diverge wildly -- a real
        // convergence sanity check, not just "it produced some number"
        Check(std::fabs(cbVal - 1.0) < 0.5,
              "crystalball intercept is in a sane range around 1");
        Check(std::fabs(cbVal - gaussVal) < 0.3,
              "crystalball intercept tracks gauss reasonably closely on "
              "this synthetic sample (no wild divergence)");
      }
    }
    if (fIn) {
      fIn->Close();
    }
  }

  std::cout << "\n[2] JER mode: crystalball width-derived intercept exists"
            << std::endl;
  {
    std::remove(outputPath);
    runCalibration(dataPath, mcPath, outputPath, CalibrationMode::JER);

    TFile *fIn = TFile::Open(outputPath, "read");
    Check(fIn && !fIn->IsZombie(), "output file opens (JER)");
    if (fIn && !fIn->IsZombie()) {
      TH1D *hCB = (TH1D *)fIn->Get(
          "ak4PF/ak4PF_intercept_jer_abseta_ptavg_100_175_crystalball");
      Check(hCB != nullptr, "crystalball intercept_jer object exists");
      if (hCB) {
        const double cbVal = hCB->GetBinContent(1);
        Check(std::isfinite(cbVal) && cbVal > 0.0,
              "crystalball JER SF converged to a finite, positive value");
      }
    }
    if (fIn) {
      fIn->Close();
    }
  }

  std::remove(dataPath);
  std::remove(mcPath);
  std::remove(outputPath);

  std::cout << "\n=== " << nPass << " passed, " << nFail
            << " failed ===" << std::endl;
  return nFail > 0 ? 1 : 0;
}
