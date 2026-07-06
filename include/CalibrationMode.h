#ifndef CALIBRATIONMODE_H
#define CALIBRATIONMODE_H

// Shared between runCalibration (Step 2) and runTextFile (Step 3): which
// calibration quantity a pass produces (Step 2) or expects to find already
// written (Step 3).
//   JEC -- mean of the asymmetry distribution -> L2Residual correction.
//   JER -- stddev of the asymmetry distribution -> JER scale factor.
// A Step 2 output file only ever has one mode's objects; Step 3 must be
// told which one it's reading, mirroring Step 2's own mandatory choice.
enum class CalibrationMode { JEC, JER };

#endif
