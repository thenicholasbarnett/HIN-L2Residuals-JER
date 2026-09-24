#ifndef FORESTEVENTLOOP_H
#define FORESTEVENTLOOP_H

#include "TDirectory.h"
#include "TFile.h"
#include "TString.h"
#include "Rtypes.h"

#include "Dijet.h"
#include "EventStructs.h"
#include "JetStruct.h"

#include <memory>
#include <vector>

enum class RunMode { MC, NonTriggered, Triggered };

// "mc" | "non-triggered" | "triggered", throws otherwise
RunMode ParseRunMode(const TString &modeFlag);

// Shared HiForest skeleton for the per-event binaries (runAsymmetry,
// runDijetProfiles): opens the forest, binds branches, and applies the
// event selection
//   |vz| < 15 cm, PV filter (data), golden JSON (data),
//   trigger-cone nref >= 2, HLT bit + plateau cut (triggered)
// then corrects every jet (JEC chain, plus JER smearing for the MC closure)
// and pT-orders each cone. Next() only stops on events that passed.
class ForestEventLoop {
public:
  static constexpr Int_t kNRefMax = 200;
  using Jets = JetStruct<kNRefMax>;

  // throws std::runtime_error on bad config or input
  ForestEventLoop(const TString &input, RunMode mode, Long64_t maxEvents = -1,
                  bool jerClosure = false);
  ~ForestEventLoop();

  bool Next();

  RunMode Mode() const { return mode_; }
  bool IsMC() const { return mode_ == RunMode::MC; }
  size_t NCones() const { return jets_.size(); }

  const EventStruct &Event() const { return event_; }
  float Weight() const { return event_.w; }
  const Jets &JetsIn(size_t c) const { return jets_[c]; }
  const float *CorrPt(size_t c) const { return corrPt_[c].data(); }
  const SortedJets &Sorted(size_t c) const { return sorted_[c]; }

  // tight jet ID and outside the veto map
  bool GoodJet(size_t c, int j) const;
  // JME Run 3 event veto (JetSelector::VetoEvent) on cone c's corrected jets:
  // any jet with pT > 15, tight ID, EMF < 0.9 inside the veto map
  bool EventVetoed(size_t c) const;

  // vz, filter, trigger QA hists, at the top of the output file
  void WriteEventHists(TDirectory *dir) const;

private:
  bool SelectEntry(Long64_t i);

  struct Impl;
  std::unique_ptr<Impl> impl_;

  RunMode mode_;
  bool jerClosure_;
  Long64_t entry_ = -1;
  Long64_t nLoop_ = 0;

  std::vector<Jets> jets_;
  EventStruct event_;
  FiltersStruct filters_;
  Int_t hltJ80_ = 0;
  std::vector<std::vector<float>> corrPt_;
  std::vector<SortedJets> sorted_;
};

// output ROOT file, parent directory created as needed
TFile *CreateOutputFile(const TString &path);

#endif
