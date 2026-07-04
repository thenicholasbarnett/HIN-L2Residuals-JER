#ifndef TEXTFILEWRITER_H
#define TEXTFILEWRITER_H

#include "TString.h"

// Step 3 entry point. Processes every cone in cfg.coneLabels in one call.
//
// triggeredResidualsFile, nonTriggeredResidualsFile: Step 2 outputs to merge.
//   Per pT_avg slice, triggeredResidualsFile is used if slice.lo >=
//   cfg.hltJ80Thresh, else nonTriggeredResidualsFile -- triggered data is
//   trigger-biased below its efficiency plateau.
// outputRootFile: per-cone corrfinal TH2D + ptcorr fit graphs, caller path.
// outputTag: plain filename prefix (no '/'), default "L2Residual". Text
//   files always land in <repo root>/data/jec/preliminary/ (gitignored,
//   created automatically), never at outputRootFile's path --
//   "<tag>_<cone>_abseta[_norm].txt" and/or "<tag>_<cone>_eta[_norm].txt"
//   per [step3] eta_mode ("both" | "abseta" | "eta").
// method: "gauss" | "doubleGauss" | "trunc90" | "trunc95".
// useNorm: kFSR-normalized intercepts (default true).
//
// Also writes JER SF text files alongside the JEC ones
// ("<tag>_<cone>_abseta_jer[_norm].txt", etc.) from the same intercept_jer_*
// histograms JER scale factors come from (see CalibrationExtractor.cxx) --
// a flat value per (eta bin, pT_avg slice) cell, real CMS JER text format,
// via the vendored JME::JetResolutionObject (external/jetmet_jer/).
void runTextFile(TString triggeredResidualsFile,
                 TString nonTriggeredResidualsFile, TString outputRootFile,
                 TString outputTag = "", TString method = "gauss",
                 bool useNorm = true);

// Scoped enum, not a bool -- a bool here once collided with the two-file
// overload above under C++ overload resolution (a TString/const char*
// argument implicitly converts to bool, which outranked the
// const char*->TString conversion the two-file overload needed at the same
// position, silently sending every two-file call to this overload instead).
// A scoped enum has no implicit string conversion, so that's a compile
// error now, not a runtime footgun.
enum class SingleDatasetKind { Triggered, NonTriggered };

// Single-dataset mode, no merge -- the two cases aren't interchangeable,
// since only one is trigger-biased:
//   Triggered    -- residualsFile is itself trigger-biased (e.g. a
//     HardProbes-only sample). pT_avg slices below cfg.hltJ80Thresh are
//     dropped entirely -- no non-triggered fallback exists to fill them in.
//   NonTriggered -- residualsFile isn't trigger-biased (e.g. ZeroBias/MinBias-
//     only). Every pT_avg slice is used unconditionally.
void runTextFile(TString residualsFile, SingleDatasetKind kind,
                 TString outputRootFile, TString outputTag = "",
                 TString method = "gauss", bool useNorm = true);

// Reads a runResponse output file and writes a CMS JER pT-resolution text
// file per cone -- JER = sigma/mean from the per-|eta|-bin Gaussian fits in
// each cone's JER_per_abseta/ subdirectory ("corr" variant, "incl"
// collection). Same flat-grid format as the JER SF file
// ({1 JetEta 1 JetPt [0] Resolution}), written to
// "<tag>_<cone>_abseta_ptresolution.txt" (and the eta variant per eta_mode).
void runTextFilePtResolution(TString responseFile, TString outputTag = "");

#endif
