#ifndef RUNASYMMETRY_H
#define RUNASYMMETRY_H

#include "TString.h"
#include "Rtypes.h"

// jerClosure: JER-smear MC jets before histogramming (hybrid method,
// JetSmearer.h), for the JER SF closure check -- only meaningful when
// modeFlag == "mc"; caller (macros/runAsymmetry.C) gates this on
// -calibration jer -closure true.
void runAsymmetry(TString input, TString output, TString modeFlag = "triggered",
                  Long64_t maxEvents = -1, bool jerClosure = false);

#endif
