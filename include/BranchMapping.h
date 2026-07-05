#ifndef BRANCHMAPPING_H
#define BRANCHMAPPING_H

#include "TTree.h"
#include "TString.h"

#include <stdexcept>
#include <vector>
#include <utility>

void SetBranches(TTree *t,
                 const std::vector<std::pair<TString, void *>> &branches) {
  t->SetBranchStatus("*", 0);

  for (const auto &branch : branches) {
    TString sBranchName = branch.first;
    void *address = branch.second;

    if (!t->GetBranch(sBranchName)) {
      throw std::runtime_error(
          Form("ERROR: Branch '%s' not found in TTree '%s'", sBranchName.Data(),
               t->GetName()));
    }

    t->SetBranchStatus(sBranchName, 1);
    t->SetBranchAddress(sBranchName, address);
  }
}

#endif
