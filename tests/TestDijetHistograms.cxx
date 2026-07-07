#include <iostream>
#include <cmath>

#include "TROOT.h"
#include "THnSparse.h"

#include "Binning.h"
#include "DijetHistograms.h"
#include "Dijet.h"

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

int main() {
  gROOT->SetBatch(true);

  std::cout << "=== TestDijetHistograms ===" << std::endl;

  BinningConfig bins;

  std::cout << "\n[1] Data mode: no response histograms" << std::endl;
  {
    ConeHistograms h;
    h.Init("test_data", bins, false);
    Check(h.hInclJetResp == nullptr, "hInclJetResp null in data mode");
    Check(h.hTagJetResp == nullptr, "hTagJetResp null in data mode");
    Check(h.hProbeJetResp == nullptr, "hProbeJetResp null in data mode");
  }

  std::cout << "\n[2] MC mode: response histograms created" << std::endl;
  ConeHistograms h;
  h.Init("test_mc", bins, true);
  Check(h.hInclJetResp != nullptr, "hInclJetResp created in MC mode");
  Check(h.hTagJetResp != nullptr, "hTagJetResp created in MC mode");
  Check(h.hProbeJetResp != nullptr, "hProbeJetResp created in MC mode");

  // no float->double rounding ambiguity with FindBin
  std::cout << "\n[3] FillInclJetResp: matched jet" << std::endl;
  {
    // corrPt=105, rawPt=88, jtPt=100, eta=0.13, refPt=80
    // -> corr=1.3125, jtpt=1.25, raw=1.1
    h.FillInclJetResp(105.0f, 88.0f, 100.0f, 0.13f, 80.0f, 1.0f);
    Double_t x[5] = {0.13, 80.0, 1.3125, 1.25, 1.1};
    Int_t coords[5];
    for (Int_t i = 0; i < 5; i++)
      coords[i] = h.hInclJetResp->GetAxis(i)->FindBin(x[i]);
    Double_t content = h.hInclJetResp->GetBinContent(coords);
    Check(std::fabs(content - 1.0) < kEps,
          "matched jet fills all three response ratios correctly");
  }

  std::cout << "\n[4] FillInclJetResp: unmatched jet skipped" << std::endl;
  {
    Long64_t before = h.hInclJetResp->GetNbins();
    h.FillInclJetResp(50.0f, 45.0f, 48.0f, 1.0f, -1.0f, 1.0f);
    Long64_t after = h.hInclJetResp->GetNbins();
    Check(before == after, "unmatched jet (refPt<0) does not add a filled bin");
  }

  std::cout << "\n[5] FillResp: tag matched, probe unmatched" << std::endl;
  {
    float pt[] = {117.0f, 80.0f}; // corrPt
    float rawPt[] = {90.0f, 70.0f};
    float jtPt[] = {93.0f, 75.0f};
    float eta[] = {0.3f, 2.5f};
    float refPt[] = {96.0f, -1.0f}; // tag matched, probe unmatched

    DijetResult d;
    d.tagIdx = 0;
    d.probeIdx = 1;

    Long64_t beforeTag = h.hTagJetResp->GetNbins();
    Long64_t beforeProbe = h.hProbeJetResp->GetNbins();

    h.FillResp(d, pt, rawPt, jtPt, eta, refPt, 1.0f);

    Long64_t afterTag = h.hTagJetResp->GetNbins();
    Long64_t afterProbe = h.hProbeJetResp->GetNbins();

    Check(afterTag == beforeTag + 1, "tag jet (matched) fills hTagJetResp");
    Check(afterProbe == beforeProbe,
          "probe jet (unmatched) leaves hProbeJetResp untouched");

    Double_t x[5] = {eta[0], refPt[0], pt[0] / refPt[0], jtPt[0] / refPt[0],
                     rawPt[0] / refPt[0]};
    Int_t coords[5];
    for (Int_t i = 0; i < 5; i++)
      coords[i] = h.hTagJetResp->GetAxis(i)->FindBin(x[i]);
    Double_t content = h.hTagJetResp->GetBinContent(coords);
    Check(std::fabs(content - 1.0) < kEps,
          "tag response ratios = corrPt/refPt=1.21875, jtPt/refPt=0.96875, "
          "rawPt/refPt=0.9375");
  }

  std::cout << "\n=== " << nPass << " passed, " << nFail
            << " failed ===" << std::endl;
  return nFail > 0 ? 1 : 0;
}
