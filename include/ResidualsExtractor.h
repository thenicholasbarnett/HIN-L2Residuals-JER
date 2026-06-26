#ifndef RESIDUALSEXTRACTOR_H
#define RESIDUALSEXTRACTOR_H

#include "TString.h"

// Step 2 entry point.
// dataFile / mcFile: hadded ROOT files from Step 1 (runAsymmetry output)
// outputFile:        TGraphErrors + fit functions + correction TH1Ds
void runResiduals(TString dataFile, TString mcFile, TString outputFile);

#endif
