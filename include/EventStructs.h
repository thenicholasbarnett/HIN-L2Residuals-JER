#ifndef EVENTSTRUCTS_H
#define EVENTSTRUCTS_H

#include "TString.h"
#include "Rtypes.h"

#include <vector>
#include <utility>

struct EventStruct {

  // weight is 1 for data and is mapped to a branch for MC
  Float_t w = 1.0f;
  Float_t vz;
  UInt_t run;
  ULong64_t event;
  UInt_t lumi;

  std::vector<std::pair<TString, void *>> BranchMap(bool isMC) {
    std::vector<std::pair<TString, void *>> branches = {{"vz", &vz},
                                                        {"evt", &event}};
    if (isMC) {
      branches.push_back({"weight", &w});
    } else {
      branches.insert(branches.end(), {{"run", &run}, {"lumi", &lumi}});
    }
    return branches;
  }
};

struct FiltersStruct {

  // primary vertex filter
  Int_t ppvF;
  std::vector<std::pair<TString, void *>>
  BranchMap(const TString &filterBranch) {
    return {{filterBranch, &ppvF}};
  }
};

#endif
