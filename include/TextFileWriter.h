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
// outputTag:        a plain filename prefix, NOT a path — must not contain
//                    '/' (rejected with an error otherwise). Empty/omitted
//                    defaults to "L2Residual". The correction text files
//                    themselves always go to a fixed location,
//                    <repo root>/data/jec/preliminary/ (created automatically
//                    if missing) — that directory is gitignored, meant for
//                    locally-generated/preliminary corrections, not a
//                    caller-chosen path. Writes
//                    "data/jec/preliminary/<tag>_<cone>_abseta.txt"
//                    (mirrored |eta| fit) and/or
//                    "data/jec/preliminary/<tag>_<cone>_eta.txt"
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
//
// JER SF text output (added 2026-07-03): alongside the JEC text files above,
// also writes "data/jec/preliminary/<tag>_<cone>_abseta_jer[_norm].txt"
// and "..._<cone>_eta_jer[_norm].txt", using the same intercept_jer_* input
// this repo's JER scale factors already come from (see CalibrationExtractor.cxx)
// and the same eta_mode/useNorm/tag rules as the JEC files. Unlike the JEC
// writer, this is a direct binned grid (one flat value per eta-bin/pT_avg-slice
// cell, no pT-dependence fit) written in the real, standard CMS JER text
// format via the vendored JME::JetResolutionObject
// (external/jetmet_jer/, from CMSSW CondFormats/JetMETObjects) rather than a
// hand-rolled format.
void runTextFile( TString triggeredResidualsFile, TString nonTriggeredResidualsFile,
                  TString outputRootFile, TString outputTag = "",
                  TString method = "gauss", bool useNorm = true );

// A scoped enum, not a bool: a bare bool in this parameter slot collides with
// the two-file overload above under C++ overload resolution -- a TString/
// const char* argument implicitly converts to bool (pointer-to-bool, a
// built-in conversion), which ranks *better* than the const char*->TString
// user-defined conversion the two-file overload needs at that same
// position. That silently sent every two-file call to this overload instead
// (verified: it's not a hypothetical, it actually happened). A scoped enum
// has no implicit conversion from a string type at all, so it can only ever
// bind to this overload -- the ambiguity is a compile error, not a runtime
// footgun, if it's ever reintroduced.
enum class SingleDatasetKind { Triggered, NonTriggered };

// Single-dataset mode: only one dataset available, no triggered/non-triggered
// merge. kind selects which behavior applies -- the two cases are NOT
// interchangeable, since only one of them is trigger-biased:
//   Triggered    -- residualsFile is itself a trigger-biased dataset (e.g. a
//            HardProbes-only sample with no ZeroBias/MinBias companion). pT_avg
//            slices below cfg.hltJ80Thresh are dropped entirely (left missing
//            in the output, same as any other unfilled slice) -- there is no
//            non-triggered fallback to fill them in, and using the triggered
//            dataset there anyway would silently bake trigger-efficiency bias
//            into the correction.
//   NonTriggered -- residualsFile is not trigger-biased (e.g. ZeroBias/MinBias-only,
//            or any single min-bias-like sample). Every pT_avg slice is used
//            unconditionally, no threshold cut applied.
void runTextFile( TString residualsFile, SingleDatasetKind kind, TString outputRootFile,
                  TString outputTag = "",
                  TString method = "gauss", bool useNorm = true );

// Reads a runResponse output file (produced by runResponse on a Step 1 MC
// file) and writes a CMS JER pT resolution text file per cone. The resolution
// values are JER = sigma/mean from the per-|eta|-bin Gaussian fits stored in
// the "JER_per_abseta/" subdirectory of each cone directory. Uses the "corr"
// variant (this framework's L2Relative+L2Residual-corrected pT) and the
// "incl" collection. Output format is the same flat-grid format as the JER SF
// file ({1 JetEta 1 JetPt [0] Resolution}), with resolution values in place
// of scale factors. Files written to data/jec/preliminary/ as
// "<tag>_<cone>_abseta_ptresolution.txt". eta_mode from the TOML controls
// whether the abseta (mirrored) variant, the independent full-eta variant, or
// both are written -- same rule as the JEC/JER SF writers.
void runTextFilePtResolution( TString responseFile, TString outputTag = "" );

#endif
