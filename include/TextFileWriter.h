#ifndef TEXTFILEWRITER_H
#define TEXTFILEWRITER_H

#include "TString.h"

// Step 3 entry point.
// residualsFile: output of runResiduals (Step 2)
// outputFile:    destination JEC text file (standard CMS format, readable by JetCorrector)
// method:        "gauss" | "trunc90" | "trunc95"  (default: "gauss")
// cone:          cone label from kConeLabels, e.g. "ak4PF"  (default: "ak4PF")
void runTextFile(TString residualsFile, TString outputFile,
                 TString method = "gauss", TString cone = "ak4PF");

#endif
