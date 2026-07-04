#ifndef RESPONSEEXTRACTOR_H
#define RESPONSEEXTRACTOR_H

#include "TString.h"

// Reads a hadded Step 1 MC file (runAsymmetry MODE=mc) and, for each cone
// and each of the three matched jet collections (incl/tag/probe), extracts
// JES/JER for three response ratios read off the 5D (eta_reco, pT_gen,
// corrPt/pT_gen, jtPt/pT_gen, rawPt/pT_gen) response THnSparse (see
// DijetHistograms.h):
//   corr — this framework's L2Relative+L2Residual-corrected pT
//   reco — the ntuple's own baked-in "jtpt" correction
//   raw  — the ntuple's uncorrected "rawpt"
// JES = mean of a Gaussian fit to the ratio, JER = sigma/mean of the same
// fit (fractional resolution, standard CMS convention), as a function of
// pT_gen (marginal over eta_reco and the other two ratio axes). Binning is
// on the GEN pT, never reco: conditioning on truth directly is what makes
// this an unbiased resolution measurement rather than a data-applicable
// correction curve (binning by reco pT would bias the result via
// falling-spectrum migration). eta_reco is deliberately NOT gen-binned here
// (Nicky's explicit call, 2026-07-02) -- the correction is applied to a jet
// by its reconstructed eta in data, so the axis is kept in the sparse for a
// future eta_reco-binned extraction pass, not marginalized away, but this
// pass doesn't slice on it yet. No data-mode input exists -- the response
// THnSparses this reads are MC-only.
void runResponse(TString inputFile, TString outputFile);

#endif
