#ifndef CALIBRATIONEXTRACTOR_H
#define CALIBRATIONEXTRACTOR_H

#include "TString.h"

// Reads hadded Step 1 output (one data file, one MC file) and, for each
// (alpha threshold, pT_avg slice, eta bin), extracts both:
//   JEC    -- the mean of the asymmetry distribution, R = (1+<A>)/(1-<A>),
//             ratio R_MC/R_data, linear fit vs alpha, extrapolate to alpha=0
//             (the L2Residual correction; not to be confused with runResponse's
//             MC-only JES, a different quantity).
//   JER SF -- the same procedure, fed the *stddev* of the asymmetry
//             distribution instead of the mean: R = (1+stddev_A)/(1-stddev_A),
//             same ratio/fit/extrapolation machinery, same 4 methods
//             (Gauss, double-Gauss, trunc90, trunc95), same kFSR-normalized
//             variant. Both ride the same per-bin A-distribution projection
//             and fits -- JER SF is not a separate pass over the input.
void runCalibration( TString dataFile, TString mcFile, TString outputFile );

#endif
