#ifndef RESPONSEEXTRACTOR_H
#define RESPONSEEXTRACTOR_H

#include "TString.h"

// Reads a hadded Step 1 MC file (runAsymmetry MODE=mc) and, for each cone
// and each of the three matched jet collections (incl/tag/probe), extracts:
//   JES = mean of a Gaussian fit to the pT_reco/pT_gen response
//   JER = sigma/mean of the same fit (fractional resolution, standard CMS
//         convention)
// as a function of eta_gen (both |eta_gen| and full eta_gen) and, separately,
// as a function of pT_gen. Binning is always on the GEN quantity, never
// reco: conditioning on truth directly is what makes this an unbiased
// resolution measurement rather than a data-applicable correction curve
// (binning by reco pT would bias the result via falling-spectrum migration;
// jointly binning both pT_reco and pT_gen at once collapses the response
// spread to just the bin width). No data-mode input exists -- the response
// THnSparses this reads are MC-only (see DijetHistograms.h).
void runResponse( TString inputFile, TString outputFile );

#endif
