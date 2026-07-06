#ifndef TEXTFILEWRITER_H
#define TEXTFILEWRITER_H

#include "TString.h"

#include "CalibrationMode.h"

// Step 3 entry point, merges Step 2 outputs per cone in cfg.coneLabels.
// mode mirrors runCalibration's own mandatory mode -- a Step 2 file only
// ever has one calibration kind's objects, so Step 3 must be told which one
// to read: JEC writes the L2Residual correction text files (corrfinal grid +
// pT-dependence fit, from intercept_*); JER writes the JER SF text files
// (flat per-slice grid, from intercept_jer_*, real CMS JER text format via
// the vendored JME::JetResolutionObject). No default -- required, same as
// Step 2.
// Per pT_avg slice: triggeredResidualsFile if slice.lo >= cfg.hltJ80Thresh,
// else nonTriggeredResidualsFile.
// outputTag: filename prefix (no '/'), default "L2Residual". Text files land
// in <repo root>/data/jec/preliminary/, independent of outputRootFile's path.
// method: "gauss" | "doubleGauss" | "trunc90" | "trunc95".
void runTextFile(TString triggeredResidualsFile,
                 TString nonTriggeredResidualsFile, TString outputRootFile,
                 CalibrationMode mode, TString outputTag = "",
                 TString method = "gauss", bool useNorm = true);

// scoped enum, not a bool -- a bool here once silently collided with the
// two-file overload above under C++ overload resolution
enum class SingleDatasetKind { Triggered, NonTriggered };

// Single-dataset mode, no merge -- the two kinds aren't interchangeable:
// Triggered residualsFile is itself trigger-biased, so slices below
// cfg.hltJ80Thresh are dropped (no fallback); NonTriggered isn't
// trigger-biased, so every slice is used.
void runTextFile(TString residualsFile, SingleDatasetKind kind,
                 TString outputRootFile, CalibrationMode mode,
                 TString outputTag = "", TString method = "gauss",
                 bool useNorm = true);

// Reads a runResponse output file, writes a CMS JER pT-resolution text file
// per cone (JER = sigma/mean, "corr" variant, "incl" collection). Same
// flat-grid format as the JER SF file.
void runTextFilePtResolution(TString responseFile, TString outputTag = "");

#endif
