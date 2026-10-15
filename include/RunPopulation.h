#ifndef RUNPOPULATION_H
#define RUNPOPULATION_H

#include "TString.h"
#include "Rtypes.h"

// jet populations of the first configured cone (PF, calo, gen in MC), its
// correction stages and the MC response per stage; MC pT-hat spikes dropped,
// optional leading / subleading jet cuts from the config's [population]
void runPopulation(TString input, TString output, TString modeFlag = "triggered",
                   Long64_t maxEvents = -1);

#endif
