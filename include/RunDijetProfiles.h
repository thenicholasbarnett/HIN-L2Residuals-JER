#ifndef RUNDIJETPROFILES_H
#define RUNDIJETPROFILES_H

#include "TString.h"
#include "Rtypes.h"

// Same forest, event selection and jet corrections as runAsymmetry
// (ForestEventLoop.h), writing inputs for JME's dijet tools instead of
// asymmetry sparses: l2l3res-multijet profiles (DijetProfiles.h) and
// DiJetJERC sparses (DiJetJERCHistograms.h). jerClosure as in runAsymmetry.
void runDijetProfiles(TString input, TString output,
                      TString modeFlag = "triggered", Long64_t maxEvents = -1,
                      bool jerClosure = false);

#endif
