#ifndef RUNPOPULATION_H
#define RUNPOPULATION_H

#include "TString.h"
#include "Rtypes.h"

// one 7D THnSparse of every jet of the first configured cone:
// (eta, phi, corrected pT, passes jet ID, outside veto map, event not vetoed,
// leading jet)
void runPopulation(TString input, TString output, TString modeFlag = "triggered",
                   Long64_t maxEvents = -1);

#endif
