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
// useNorm:          use the kFSR-normalized intercepts (default: true — this
//                    is the standard method). Pass false for the direct,
//                    non-normalized variant. Affects the corrfinal grid and
//                    ptcorr graph names too (both get a "_norm" suffix when true).
void runTextFile( TString hpResidualsFile, TString zbResidualsFile,
                  TString outputRootFile, TString outputTextPrefix,
                  TString method = "gauss", bool useNorm = true );

// Single-dataset convenience overload, for systems with no HP/ZB split (e.g.
// a single min-bias or single-trigger sample — everything above threshold
// is unbiased, everything isn't, or there's simply only one dataset to use).
// Equivalent to calling the two-file version with the same file for both
// hpResidualsFile and zbResidualsFile: every pT slice reads from this one
// file regardless of cfg.hltJ80Thresh, so the merge becomes a no-op.
void runTextFile( TString residualsFile, TString outputRootFile,
                  TString outputTextPrefix,
                  TString method = "gauss", bool useNorm = true );

#endif
