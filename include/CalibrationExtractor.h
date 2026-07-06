#ifndef CALIBRATIONEXTRACTOR_H
#define CALIBRATIONEXTRACTOR_H

#include "TString.h"

#include "CalibrationMode.h"

// Both modes run the same per-bin A-distribution projection and fitting
// (Gaussian fit yields mean and sigma together), but only the requested
// quantity's R(alpha) points, R_data/R_mc histograms, and intercept TH1Ds
// are accumulated and written. Mode is required -- there is no default.

// Reads hadded Step 1 output (one data file, one MC file). For each
// (alpha threshold, pT_avg slice, eta bin), runs the full extraction loop
// for the selected mode and writes intercept TH1Ds + R(alpha) graphs.
void runCalibration(TString dataFile, TString mcFile, TString outputFile,
                    CalibrationMode mode);

#endif
