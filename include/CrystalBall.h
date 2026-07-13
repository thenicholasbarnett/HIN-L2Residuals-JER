#ifndef CRYSTAL_BALL_H
#define CRYSTAL_BALL_H

#include "TMath.h"

// Double-sided Crystal Ball: a Gaussian core with independent power-law
// tails on each side, ported from the real fit used by CMS's official
// JRA-based JER tooling (JetMETAnalysis's jet_response_fitter_x.cc,
// fnc_dscb) -- ported here as an alternative to this repo's own Gauss/
// double-Gauss/truncated-mean estimators, shared by Step 2's A-distribution
// fitting (CalibrationExtractor.cxx) and runResponse's JES/JER extraction
// (ResponseExtractor.cxx), both of which window it around their own peak
// (0 for A, 1 for the response ratio) via TF1's fit range, not a parameter
// here.
//
// Deliberately just the raw function -- the reference implementation's
// TSpectrum peak-finding, histogram-name-driven fit-range parsing, and
// 4-iteration progressive parameter-fixing loop are all CMS-ntuple-specific
// or a bigger deviation from this repo's single-pass TF1::Fit(f, "NQSR")
// convention than warranted for a first pass; add iteration later only if
// convergence proves poor in practice.
//
// Parameters: [0] N (norm), [1] mu (mean), [2] sigma, [3] a1/[4] p1 (low-side
// tail), [5] a2/[6] p2 (high-side tail) -- independent per side, genuinely
// asymmetric.
inline double DoubleSidedCB(double *xx, double *pp) {
  double x = xx[0];
  double N = pp[0];
  double mu = pp[1];
  double sig = pp[2];
  double a1 = pp[3];
  double p1 = pp[4];
  double a2 = pp[5];
  double p2 = pp[6];

  double u = (x - mu) / sig;
  double A1 = TMath::Power(p1 / TMath::Abs(a1), p1) * TMath::Exp(-a1 * a1 / 2);
  double A2 = TMath::Power(p2 / TMath::Abs(a2), p2) * TMath::Exp(-a2 * a2 / 2);
  double B1 = p1 / TMath::Abs(a1) - TMath::Abs(a1);
  double B2 = p2 / TMath::Abs(a2) - TMath::Abs(a2);

  double result = N;
  if (u < -a1) {
    result *= A1 * TMath::Power(B1 - u, -p1);
  } else if (u < a2) {
    result *= TMath::Exp(-u * u / 2);
  } else {
    result *= A2 * TMath::Power(B2 + u, -p2);
  }
  return result;
}

#endif
