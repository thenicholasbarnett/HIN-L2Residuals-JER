#ifndef JETSTRUCT_H
#define JETSTRUCT_H

#include "TString.h"
#include "Rtypes.h"

#include <vector>
#include <utility>

template <Int_t MAXNREF> struct JetStruct {

  // upper bound on jets per event across all input ttrees
  static constexpr Int_t maxnref = MAXNREF;

  // reco jets
  struct RecoMomenta {
    // number in event
    Int_t nref;
    // momenta
    Float_t rawpt[MAXNREF];
    Float_t pt[MAXNREF];
    Float_t eta[MAXNREF];
    Float_t phi[MAXNREF];
    Float_t area[MAXNREF];
    struct PFInfo {
      // PF fractions
      Float_t CHF[MAXNREF];
      Float_t NHF[MAXNREF];
      Float_t CEF[MAXNREF];
      Float_t NEF[MAXNREF];
      Float_t MUF[MAXNREF];
      // PF multiplicities
      Int_t CHM[MAXNREF];
      Int_t NHM[MAXNREF];
      Int_t CEM[MAXNREF];
      Int_t NEM[MAXNREF];
      Int_t MUM[MAXNREF];
    } pf;
  } reco;

  // ref jets
  struct RefMomenta {
    Float_t pt[MAXNREF];
    Float_t eta[MAXNREF];
    Float_t phi[MAXNREF];
  } ref;

  // gen jets
  struct GenMomenta {
    Float_t pt[MAXNREF];
    Float_t eta[MAXNREF];
    Float_t phi[MAXNREF];
  } gen;

  std::vector<std::pair<TString, void *>> BranchMap(bool isMC) {
    std::vector<std::pair<TString, void *>> map = {
        {"nref", &reco.nref},     {"rawpt", reco.rawpt},
        {"jtpt", reco.pt},        {"jteta", reco.eta},
        {"jtphi", reco.phi},      {"jtarea", reco.area},
        {"jtPfCHF", reco.pf.CHF},
        {"jtPfNHF", reco.pf.NHF}, {"jtPfCEF", reco.pf.CEF},
        {"jtPfNEF", reco.pf.NEF}, {"jtPfMUF", reco.pf.MUF},
        {"jtPfCHM", reco.pf.CHM}, {"jtPfNHM", reco.pf.NHM},
        {"jtPfCEM", reco.pf.CEM}, {"jtPfNEM", reco.pf.NEM},
        {"jtPfMUM", reco.pf.MUM},
    };

    // including gen and ref jet collections iff MC
    if (isMC) {
      map.insert(map.end(), {
                                {"refpt", ref.pt},
                                {"refeta", ref.eta},
                                {"refphi", ref.phi},
                                {"genpt", gen.pt},
                                {"geneta", gen.eta},
                                {"genphi", gen.phi},
                            });
    }
    return map;
  }
};

#endif
