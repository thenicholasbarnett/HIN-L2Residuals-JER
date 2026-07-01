#ifndef TEXTFILEWRITER_H
#define TEXTFILEWRITER_H

#include "TString.h"

// Step 3 entry point. Processes every cone in cfg.coneLabels in one call.
//
// hpResidualsFile, zbResidualsFile: Step 2 outputs for HardProbes and ZeroBias.
//   Per pT_avg slice, the slice is taken from hpResidualsFile if
//   slice.lo >= cfg.hltJ80Thresh, otherwise from zbResidualsFile — HP is
//   trigger-biased below its efficiency plateau, ZB fills in the rest.
// outputRootFile:   per-cone TDirectories, each with a corrfinal_{etaMode}_{method}
//                    TH2D (eta/|eta| vs pT_avg, z = final correction) and a
//                    graphs/ dir of per-eta-bin pT-dependence fit TGraphErrors.
// outputTextPrefix: base path; writes "<prefix>_<cone>_abseta.txt" (mirrored
//                    |eta| fit) and "<prefix>_<cone>_eta.txt" (independent
//                    full-eta fit, no mirroring) per cone. When useNorm is
//                    true, both get a "_norm" suffix before ".txt".
// method:           "gauss" | "doubleGauss" | "trunc90" | "trunc95"  (default: "gauss")
// useNorm:          use the kFSR-normalized intercepts instead of the direct
//                    ones (default: false). Affects the corrfinal grid and
//                    ptcorr graph names too (both get a "_norm" suffix).
void runTextFile( TString hpResidualsFile, TString zbResidualsFile,
                  TString outputRootFile, TString outputTextPrefix,
                  TString method = "gauss", bool useNorm = false );

#endif
