#ifndef RUNASYMMETRY_H
#define RUNASYMMETRY_H

#include "TString.h"
#include "Rtypes.h"

// declaration header

void runAsymmetry(TString input, TString output, TString modeFlag = "triggered",
                  Long64_t maxEvents = -1);

#endif
