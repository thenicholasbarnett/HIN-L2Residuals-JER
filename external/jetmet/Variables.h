// Vendored verbatim from JetMETAnalysis/JetUtilities/interface/Variables.hh
// (https://github.com/lmartika/JetMETAnalysis) into L2Residuals-2024ppref on
// 2026-07-03. Zero CMSSW coupling -- just ROOT basics (TROOT/TSystem/TString).
// Used for JME-standard axis titles wherever this repo's own plot quantity has
// a clean, honest equivalent in JetMET's single-jet pt/eta vocabulary (gen-jet
// pT in runResponse's JES/JER-vs-pT_gen plots -> refpt; reco jet pT/eta in
// Kinematics.h's incl/tag/probe plots -> jtpt/jteta). Quantities with no JetMET
// equivalent -- dijet-averaged pT_avg, alpha, the asymmetry A, response-ratio
// axes, phi (not in this enum at all) -- keep this repo's own local titles.

#ifndef VARIABLES_HH
#define VARIABLES_HH

//ROOT libraries
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"

//C++ libraries
#include <iostream>
#include <string>
#include <map>
#include <vector>

namespace VARIABLES {

   enum Variable {none, refpt, ptclpt, refeta, jtpt, jteta, mu, npv, npu, rho, area, other, unknown};
	static const unsigned int nVariable = 12;


	// A routine that returns a string given the VARIABLE
	std::string getVariableString(Variable var);

	// A routine that returns a string fit for a title
	std::string getVariableTitleString(Variable var);

	// A routine that returns a string fit for an axis title
	std::string getVariableAxisTitleString(Variable var, bool unit=true);

	// A routine that returns a VARIABLE given a string
	Variable getVariable(std::string str);

}

#endif
