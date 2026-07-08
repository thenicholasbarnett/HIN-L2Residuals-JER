#ifndef TRUNCATION_H
#define TRUNCATION_H

#include "TH1D.h"
#include "TMath.h"

#include "Utilities.h" // CanFit

#include <utility>

// Central-fraction truncated mean/RMS, shared by Step 2's A-distribution
// extraction (CalibrationExtractor) and runResponse's JES/JER extraction.

struct TruncResult {
  double mean = 0, meanErr = 0, sigma = 0, sigmaErr = 0, nEff = 0;
  bool valid = false;
};

// Bin range holding the central `fraction` of the integral, trimming
// half the excluded tail off each side. {1, 0} (empty) if h can't be fit.
inline std::pair<int, int> FindTruncBins(TH1D *h, double fraction,
                                         int minEntries) {
  if (!CanFit(h, minEntries)) {
    return {1, 0};
  }
  double total = h->Integral();
  if (total <= 0) {
    return {1, 0};
  }
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
inline TruncResult TruncMeanInRange(TH1D *h, int binLo, int binHi) {
  TruncResult r;
  if (binLo > binHi) {
    return r;
  }
  double nEff = h->Integral(binLo, binHi);
  double hTotal = h->Integral();
  double nEffEntries =
      (hTotal > 1e-10) ? h->GetEntries() * nEff / hTotal : nEff;
  // no separate entries gate here: FindTruncBins already required
  // h's raw (pre-truncation) entries >= minEntries
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

#endif
