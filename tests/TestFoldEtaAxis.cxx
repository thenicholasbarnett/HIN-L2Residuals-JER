#include <iostream>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "THnSparse.h"
#include "TMath.h"

#include "Binning.h"
#include "Utilities.h"

static int nPass = 0;
static int nFail = 0;
static constexpr double kEps = 1e-9;

void Check(bool cond, const char *msg) {
  if (cond) {
    std::cout << "  PASS  " << msg << std::endl;
    nPass++;
  } else {
    std::cout << "  FAIL  " << msg << std::endl;
    nFail++;
  }
}

// Build a 3-axis THnSparse: full eta (36 CMS JEC bins) x asymmetry (10 bins) x pT_avg (5 bins)
THnSparse *MakeTestSparse(const char *name) {
  BinningConfig bins;
  AxisBins asymBins = {10, -1.0, 1.0, "A"};
  AxisBins ptBins = {5, 40.0, 1000.0, "p_{T,avg}"};
  THnSparse *h =
      MakeTHnSparse<THnSparseD>(name, "test", {bins.eta, asymBins, ptBins});
  SetEtaBins(h, 0);
  h->Sumw2();
  return h;
}

int main() {
  gROOT->SetBatch(true);

  std::cout << "=== TestFoldEtaAxis ===" << std::endl;

  // ----------------------------------------------------------------
  // Test 1: input axis properties before folding
  // ----------------------------------------------------------------
  std::cout << "\n[1] Input axis properties" << std::endl;
  {
    THnSparse *h = MakeTestSparse("h1");
    Check(h->GetAxis(0)->GetNbins() == 36, "eta axis has 36 bins");
    Check(std::abs(h->GetAxis(0)->GetXmin() - (-5.191)) < kEps,
          "eta axis lo = -5.191");
    Check(std::abs(h->GetAxis(0)->GetXmax() - 5.191) < kEps,
          "eta axis hi =  5.191");
    Check(std::abs(h->GetAxis(0)->GetBinLowEdge(1) - kEtaEdges[0]) < kEps,
          "bin 1 low edge matches kEtaEdges[0]");
    Check(std::abs(h->GetAxis(0)->GetBinUpEdge(36) - kEtaEdges[36]) < kEps,
          "bin 36 up edge matches kEtaEdges[36]");
    delete h;
  }

  // ----------------------------------------------------------------
  // Test 2: folded axis properties
  // ----------------------------------------------------------------
  std::cout << "\n[2] Folded axis properties" << std::endl;
  {
    THnSparse *h = MakeTestSparse("h2");
    THnSparse *hf = FoldEtaAxis(h, 0);

    Check(hf->GetAxis(0)->GetNbins() == 18, "|eta| axis has 18 bins");
    Check(std::abs(hf->GetAxis(0)->GetXmin() - 0.0) < kEps,
          "|eta| axis lo = 0");
    Check(std::abs(hf->GetAxis(0)->GetXmax() - 5.191) < kEps,
          "|eta| axis hi = 5.191");
    Check(std::abs(hf->GetAxis(0)->GetBinLowEdge(1) - kAbsEtaEdges[0]) < kEps,
          "folded bin 1 low edge matches kAbsEtaEdges[0]");
    Check(std::abs(hf->GetAxis(0)->GetBinUpEdge(1) - kAbsEtaEdges[1]) < kEps,
          "folded bin 1 up edge = 0.261");
    Check(std::abs(hf->GetAxis(0)->GetBinUpEdge(18) - kAbsEtaEdges[18]) < kEps,
          "folded bin 18 up edge = 5.191");

    // non-eta axes unchanged
    Check(hf->GetAxis(1)->GetNbins() == 10, "asymmetry axis: still 10 bins");
    Check(hf->GetAxis(2)->GetNbins() == 5, "pT_avg axis: still 5 bins");
    Check(std::abs(hf->GetAxis(1)->GetXmin() - (-1.0)) < kEps,
          "asymmetry axis lo unchanged");
    Check(std::abs(hf->GetAxis(2)->GetXmax() - 1000.0) < kEps,
          "pT_avg axis hi unchanged");

    delete h;
    delete hf;
  }

  // ----------------------------------------------------------------
  // Test 3: symmetric fill — both sides fold into the same |eta| bin
  //   fill eta = +0.13 and eta = -0.13 with weight 1.0 each
  //   both land in |eta| bin 1 [0, 0.261]
  //   expect: content = 2.0, error = sqrt(1^2 + 1^2) = sqrt(2)
  // ----------------------------------------------------------------
  std::cout << "\n[3] Symmetric fill" << std::endl;
  {
    THnSparse *h = MakeTestSparse("h3");
    Double_t xpos[3] = {0.13, 0.0, 100.0};
    Double_t xneg[3] = {-0.13, 0.0, 100.0};
    h->Fill(xpos, 1.0);
    h->Fill(xneg, 1.0);

    THnSparse *hf = FoldEtaAxis(h, 0);

    Int_t absEtaBin = hf->GetAxis(0)->FindBin(0.13); // should be 1
    Int_t aBin = hf->GetAxis(1)->FindBin(0.0);
    Int_t ptBin = hf->GetAxis(2)->FindBin(100.0);
    Int_t c[3] = {absEtaBin, aBin, ptBin};

    double content = hf->GetBinContent(c);
    double err = hf->GetBinError(c);

    Check(absEtaBin == 1, "eta=0.13 lands in |eta| bin 1");
    Check(std::abs(content - 2.0) < kEps, "symmetric fill: content = 2.0");
    Check(std::abs(err - TMath::Sqrt(2.0)) < kEps,
          "symmetric fill: error = sqrt(2)");

    delete h;
    delete hf;
  }

  // ----------------------------------------------------------------
  // Test 4: positive-only fill
  //   fill eta = +2.0 (|eta| bin 9, [1.930, 2.172]) with weight 3.0
  //   expect: content = 3.0, error = 3.0
  // ----------------------------------------------------------------
  std::cout << "\n[4] Positive-only fill" << std::endl;
  {
    THnSparse *h = MakeTestSparse("h4");
    Double_t x[3] = {2.0, 0.0, 100.0};
    h->Fill(x, 3.0);

    THnSparse *hf = FoldEtaAxis(h, 0);

    Int_t absEtaBin = hf->GetAxis(0)->FindBin(2.0);
    Int_t aBin = hf->GetAxis(1)->FindBin(0.0);
    Int_t ptBin = hf->GetAxis(2)->FindBin(100.0);
    Int_t c[3] = {absEtaBin, aBin, ptBin};

    Check(std::abs(hf->GetBinContent(c) - 3.0) < kEps,
          "positive fill: content = 3.0");
    Check(std::abs(hf->GetBinError(c) - 3.0) < kEps,
          "positive fill: error = 3.0");

    delete h;
    delete hf;
  }

  // ----------------------------------------------------------------
  // Test 5: negative-only fill maps to the same |eta| bin as positive
  //   fill eta = -2.0 with weight 5.0
  //   should map to same |eta| bin as eta = +2.0
  // ----------------------------------------------------------------
  std::cout << "\n[5] Negative-only fill maps to correct |eta| bin"
            << std::endl;
  {
    THnSparse *h = MakeTestSparse("h5");
    Double_t x[3] = {-2.0, 0.0, 100.0};
    h->Fill(x, 5.0);

    THnSparse *hf = FoldEtaAxis(h, 0);

    Int_t absEtaBin = hf->GetAxis(0)->FindBin(2.0);
    Int_t aBin = hf->GetAxis(1)->FindBin(0.0);
    Int_t ptBin = hf->GetAxis(2)->FindBin(100.0);
    Int_t c[3] = {absEtaBin, aBin, ptBin};

    Check(std::abs(hf->GetBinContent(c) - 5.0) < kEps,
          "negative fill: content = 5.0");
    Check(std::abs(hf->GetBinError(c) - 5.0) < kEps,
          "negative fill: error = 5.0");

    // also confirm no content anywhere else on the |eta| axis
    bool onlyOneBinFilled = true;
    for (Int_t ib = 1; ib <= 18; ib++) {
      if (ib == absEtaBin)
        continue;
      Int_t cc[3] = {ib, aBin, ptBin};
      if (hf->GetBinContent(cc) != 0.0) {
        onlyOneBinFilled = false;
        break;
      }
    }
    Check(onlyOneBinFilled,
          "negative fill: only the correct |eta| bin is filled");

    delete h;
    delete hf;
  }

  // ----------------------------------------------------------------
  // Test 6: negative and positive fills at different |eta| values
  //   fill eta = +0.5 (|eta| bin 2) and eta = -2.6 (|eta| bin 12, interior)
  //   these should NOT add together — they land in different |eta| bins.
  //   Note: -2.5 is exactly on edge kEtaEdges[7] and would go to the
  //   adjacent bin due to ROOT's bin-edge convention, so -2.6 is used instead.
  // ----------------------------------------------------------------
  std::cout << "\n[6] Fills at different |eta| values stay separate"
            << std::endl;
  {
    THnSparse *h = MakeTestSparse("h6");
    Double_t xA[3] = {0.5, 0.0, 100.0};
    Double_t xB[3] = {-2.6, 0.0, 100.0};
    h->Fill(xA, 2.0);
    h->Fill(xB, 7.0);

    THnSparse *hf = FoldEtaAxis(h, 0);

    Int_t binA = hf->GetAxis(0)->FindBin(0.5);
    Int_t binB = hf->GetAxis(0)->FindBin(2.6);
    Int_t aBin = hf->GetAxis(1)->FindBin(0.0);
    Int_t ptBin = hf->GetAxis(2)->FindBin(100.0);
    Int_t cA[3] = {binA, aBin, ptBin};
    Int_t cB[3] = {binB, aBin, ptBin};

    Check(binA != binB, "eta=0.5 and eta=-2.6 map to different |eta| bins");
    Check(std::abs(hf->GetBinContent(cA) - 2.0) < kEps,
          "eta=+0.5 bin has content 2.0");
    Check(std::abs(hf->GetBinContent(cB) - 7.0) < kEps,
          "eta=-2.6 bin has content 7.0");

    delete h;
    delete hf;
  }

  // ----------------------------------------------------------------
  // Test 7: eta axis index is not hardcoded — fold works on axis 1
  //   same symmetric fill, but eta is axis 1 instead of axis 0
  // ----------------------------------------------------------------
  std::cout << "\n[7] FoldEtaAxis works on arbitrary axis index" << std::endl;
  {
    BinningConfig bins;
    AxisBins asymBins = {10, -1.0, 1.0, "A"};
    AxisBins ptBins = {5, 40.0, 1000.0, "p_{T,avg}"};
    // axis order: asymmetry, eta, pT_avg
    THnSparse *h =
        MakeTHnSparse<THnSparseD>("h7", "test", {asymBins, bins.eta, ptBins});
    SetEtaBins(h, 1);
    h->Sumw2();

    Double_t xpos[3] = {0.0, 0.13, 100.0};
    Double_t xneg[3] = {0.0, -0.13, 100.0};
    h->Fill(xpos, 1.0);
    h->Fill(xneg, 1.0);

    THnSparse *hf = FoldEtaAxis(h, 1);

    Check(hf->GetAxis(1)->GetNbins() == 18,
          "folding axis 1: result has 18 |eta| bins");
    Check(hf->GetAxis(0)->GetNbins() == 10,
          "folding axis 1: axis 0 (asymmetry) unchanged");
    Check(hf->GetAxis(2)->GetNbins() == 5,
          "folding axis 1: axis 2 (pT_avg) unchanged");

    Int_t absEtaBin = hf->GetAxis(1)->FindBin(0.13);
    Int_t aBin = hf->GetAxis(0)->FindBin(0.0);
    Int_t ptBin = hf->GetAxis(2)->FindBin(100.0);
    Int_t c[3] = {aBin, absEtaBin, ptBin};

    Check(std::abs(hf->GetBinContent(c) - 2.0) < kEps,
          "arbitrary axis fold: content = 2.0");
    Check(std::abs(hf->GetBinError(c) - TMath::Sqrt(2.0)) < kEps,
          "arbitrary axis fold: error = sqrt(2)");

    delete h;
    delete hf;
  }

  // ----------------------------------------------------------------
  std::cout << "\n=== " << nPass << " passed, " << nFail
            << " failed ===" << std::endl;
  return nFail > 0 ? 1 : 0;
}
