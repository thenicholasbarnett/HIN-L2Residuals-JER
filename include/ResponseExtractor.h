#ifndef RESPONSEEXTRACTOR_H
#define RESPONSEEXTRACTOR_H

#include "TString.h"

// Reads a hadded Step 1 MC file (runAsymmetry MODE=mc). For each cone and
// each matched jet collection (incl/tag/probe), extracts JES/JER for three
// response ratios off the 5D response THnSparse (see DijetHistograms.h):
//   corr : this framework's L2Relative+L2Residual-corrected pT
//   reco : the ntuple's own baked-in "jtpt" correction
//   raw  : the ntuple's uncorrected "rawpt"
// JES = Gaussian fit mean, JER = sigma/mean, vs pT_gen (never reco pT, to
// avoid falling-spectrum migration bias). eta_reco stays in the sparse for
// a future binned pass, not sliced on yet. MC-only.
void runResponse(TString inputFile, TString outputFile);

#endif
