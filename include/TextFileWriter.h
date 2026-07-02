#ifndef TEXTFILEWRITER_H
#define TEXTFILEWRITER_H

#include "TString.h"

// Step 3 entry point. Processes every cone in cfg.coneLabels in one call.
//
// triggeredResidualsFile, nonTriggeredResidualsFile: Step 2 outputs for a
//   triggered dataset (e.g. HardProbes) and a non-triggered dataset (e.g.
//   ZeroBias or MinBias). Per pT_avg slice, the slice is taken from
//   triggeredResidualsFile if slice.lo >= cfg.hltJ80Thresh, otherwise from
//   nonTriggeredResidualsFile — the triggered dataset is trigger-biased below
//   its efficiency plateau, the non-triggered dataset fills in the rest.
// outputRootFile:   per-cone TDirectories, each with a corrfinal_{etaMode}_{method}
//                    TH2D (eta/|eta| vs pT_avg, z = final correction) and a
//                    graphs/ dir of per-eta-bin pT-dependence fit TGraphErrors.
//                    Path is caller-specified, same as any other Step output.
// outputTextPrefix: a plain filename prefix, NOT a path — must not contain
//                    '/' (rejected with an error otherwise). Empty/omitted
//                    defaults to "L2Residual". The correction text files
//                    themselves always go to a fixed location,
//                    <repo root>/data/jec/preliminary/ (created automatically
//                    if missing) — that directory is gitignored, meant for
//                    locally-generated/preliminary corrections, not a
//                    caller-chosen path. Writes
//                    "data/jec/preliminary/<prefix>_<cone>_abseta.txt"
//                    (mirrored |eta| fit) and/or
//                    "data/jec/preliminary/<prefix>_<cone>_eta.txt"
//                    (independent full-eta fit, no mirroring) per cone,
//                    depending on [step3] eta_mode in the TOML ("both" |
//                    "abseta" | "eta"; default "both"). When useNorm is true,
//                    both get a "_norm" suffix before ".txt".
// method:           "gauss" | "doubleGauss" | "trunc90" | "trunc95"  (default
//                    parameter here is "gauss"; macros/runTextFile.C instead
//                    falls back to [step3] default_method from the TOML)
// useNorm:          use the kFSR-normalized intercepts (default: true — this
//                    is the standard method). Pass false for the direct,
//                    non-normalized variant. Affects the corrfinal grid and
//                    ptcorr graph names too (both get a "_norm" suffix when true).
void runTextFile( TString triggeredResidualsFile, TString nonTriggeredResidualsFile,
                  TString outputRootFile, TString outputTextPrefix = "",
                  TString method = "gauss", bool useNorm = true );

// Single-dataset convenience overload, for systems with no triggered/non-triggered
// split (e.g. a single min-bias or single-trigger sample — everything above
// threshold is unbiased, everything isn't, or there's simply only one dataset
// to use). Equivalent to calling the two-file version with the same file for
// both triggeredResidualsFile and nonTriggeredResidualsFile: every pT slice
// reads from this one file regardless of cfg.hltJ80Thresh, so the merge
// becomes a no-op.
void runTextFile( TString residualsFile, TString outputRootFile,
                  TString outputTextPrefix = "",
                  TString method = "gauss", bool useNorm = true );

#endif
