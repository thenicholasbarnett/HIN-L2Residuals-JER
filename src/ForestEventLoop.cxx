#include "ForestEventLoop.h"

#include "TH1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TTree.h"

#include "jetcorrector/JetCorrector.h"
#include "JetSelector.h"
#include "json_handler/JSON_handler.h"
#include "JetSmearer.h"

#include "BranchMapping.h"
#include "AnalysisConfig.h"

#include <cstdint>
#include <random>
#include <stdexcept>
#include <string>

static constexpr float kVzCut = 15.0f;

// matches SmearedJetProducerT.h's default "seed" fillDescriptions value --
// no reason to diverge, this is purely for reproducible smearing draws
static constexpr std::uint32_t kJerClosureSeed = 37428479;

RunMode ParseRunMode(const TString &modeFlag) {
  if (modeFlag == "mc") {
    return RunMode::MC;
  }
  if (modeFlag == "non-triggered") {
    return RunMode::NonTriggered;
  }
  if (modeFlag == "triggered") {
    return RunMode::Triggered;
  }
  throw std::invalid_argument(
      Form("ERROR: invalid mode '%s'. Expected mc, non-triggered, or "
           "triggered.",
           modeFlag.Data()));
}

TFile *CreateOutputFile(const TString &path) {
  Ssiz_t sl = path.Last('/');
  if (sl != kNPOS) {
    gSystem->mkdir(TString(path(0, sl)), kTRUE);
  }
  return new TFile(path, "recreate");
}

struct ForestEventLoop::Impl {
  std::vector<JetCorrector> jecs;
  std::vector<JetSmearerJME::JetResolution> jerResolution;
  std::vector<JetSmearerJME::JetResolutionScaleFactor> jerScaleFactor;
  std::mt19937 jerRng{kJerClosureSeed};
  JetSmearing::Method jerMethod = JetSmearing::Method::Hybrid;
  std::unique_ptr<JetSelector> js;
  std::unique_ptr<JSON_handler> dcs;

  TFile *fi = nullptr;
  std::vector<TTree *> jetTrees;
  TTree *evtTree = nullptr;
  TTree *ggTree = nullptr;
  TTree *skimTree = nullptr;
  TTree *trigTree = nullptr;

  size_t trigConeIdx = 0;
  float hltJ80Thresh = 0.0f;

  TH1D *hvz_all = nullptr;
  TH1D *hvz = nullptr;
  TH1I *hfilt = nullptr;
  TH1I *h_j80 = nullptr;
};

ForestEventLoop::ForestEventLoop(const TString &input, RunMode mode,
                                 Long64_t maxEvents, bool jerClosure)
    : impl_(std::make_unique<Impl>()), mode_(mode), jerClosure_(jerClosure) {

  const AnalysisConfig &cfg = Config();
  Impl &m = *impl_;

  if (jerClosure && mode != RunMode::MC) {
    throw std::runtime_error("ERROR: -calibration jer -closure true only "
                             "applies to -mode mc (JER smearing needs "
                             "gen-matched MC jets)");
  }

  // validate cone config
  const size_t nCones = cfg.coneLabels.size();
  if (cfg.jecFilesPerCone.size() != nCones ||
      cfg.jetTreePaths.size() != nCones) {
    throw std::runtime_error("ERROR: cone labels, JEC files, and jet tree "
                             "paths must be the same length");
  }

  // find trigger cone
  m.trigConeIdx = nCones;
  for (size_t c = 0; c < nCones; c++) {
    if (cfg.coneLabels[c] == cfg.trigCone) {
      m.trigConeIdx = c;
      break;
    }
  }
  if (m.trigConeIdx == nCones) {
    throw std::runtime_error(Form("ERROR: trigger cone '%s' not found in "
                                  "cone labels",
                                  cfg.trigCone.Data()));
  }
  m.hltJ80Thresh = cfg.hltJ80Thresh;

  // L2Residuals applied to data only
  m.jecs.reserve(nCones);
  for (size_t c = 0; c < nCones; c++) {
    std::vector<std::string> chain = cfg.jecFilesPerCone[c];
    if (mode != RunMode::MC && !cfg.residualFilesPerCone.empty()) {
      for (const auto &f : cfg.residualFilesPerCone[c]) {
        chain.push_back(f);
      }
    }
    m.jecs.emplace_back(chain);
  }

  // JER SF closure: MC jets get smeared (JetSmearer.h) using previously
  // derived per-cone resolution + scale-factor text files
  if (jerClosure && (cfg.jerResolutionFilesPerCone.size() != nCones ||
                     cfg.jerScaleFactorFilesPerCone.size() != nCones)) {
    throw std::runtime_error("ERROR: -calibration jer -closure true requires "
                             "jer_closure.resolution_files and "
                             ".scale_factor_files, one per cone, in the "
                             "config");
  }
  if (jerClosure) {
    m.jerMethod = JetSmearing::MethodFromString(cfg.jerMethod);
    m.jerResolution.reserve(nCones);
    m.jerScaleFactor.reserve(nCones);
    for (size_t c = 0; c < nCones; c++) {
      m.jerResolution.emplace_back(cfg.jerResolutionFilesPerCone[c]);
      m.jerScaleFactor.emplace_back(cfg.jerScaleFactorFilesPerCone[c]);
    }
  }

  // jet ID, veto map, golden json
  m.js = std::make_unique<JetSelector>(
      cfg.jetSystem == "ion" ? JetSelector::System::Ion
                             : JetSelector::System::pp,
      cfg.jetIdPath, cfg.vetoMapPath,
      cfg.jetPurpose == "analysis" ? JetSelector::Purpose::Analysis
                                   : JetSelector::Purpose::Calibration);
  if (mode != RunMode::MC) {
    m.dcs = std::make_unique<JSON_handler>(cfg.jsonPath.Data());
  }

  // structures
  jets_.resize(nCones);
  corrPt_.assign(nCones, std::vector<float>(kNRefMax, 0.0f));
  sorted_.resize(nCones);

  // input
  m.fi = TFile::Open(input, "read");
  if (!m.fi || m.fi->IsZombie()) {
    throw std::runtime_error(Form("Cannot open %s", input.Data()));
  }

  // ttrees: jet collections, event info, pileup density, filter, trigger
  m.jetTrees.assign(nCones, nullptr);
  for (size_t c = 0; c < nCones; c++) {
    m.jetTrees[c] = (TTree *)m.fi->Get(cfg.jetTreePaths[c]);
    if (!m.jetTrees[c]) {
      throw std::runtime_error(Form("Missing jet TTree %s in %s",
                                    cfg.jetTreePaths[c].Data(), input.Data()));
    }
  }
  m.evtTree = (TTree *)m.fi->Get(cfg.hiTreePath);
  if (!m.evtTree) {
    throw std::runtime_error(
        Form("Missing HiTree %s in %s", cfg.hiTreePath.Data(), input.Data()));
  }
  m.ggTree = (TTree *)m.fi->Get(cfg.ggTreePath);
  if (!m.ggTree) {
    throw std::runtime_error(Form("Missing ggHiNtuplizer TTree %s in %s",
                                  cfg.ggTreePath.Data(), input.Data()));
  }
  if (mode != RunMode::MC) {
    m.skimTree = (TTree *)m.fi->Get(cfg.skimTreePath);
    if (!m.skimTree) {
      throw std::runtime_error(Form("Missing skim TTree %s in %s\n(check "
                                    "cfg/2024ppRef.toml)",
                                    cfg.skimTreePath.Data(), input.Data()));
    }
  }
  if (mode == RunMode::Triggered) {
    m.trigTree = (TTree *)m.fi->Get(cfg.trigTreePath);
    if (!m.trigTree) {
      throw std::runtime_error(Form("Missing HLT TTree %s in %s\n(check "
                                    "cfg/2024ppRef.toml)",
                                    cfg.trigTreePath.Data(), input.Data()));
    }
  }

  // branch mapping
  const bool isMC = IsMC();
  SetBranches(m.evtTree, event_.BranchMap(isMC));
  SetBranches(m.ggTree, event_.RhoBranchMap());
  for (size_t c = 0; c < nCones; c++) {
    SetBranches(m.jetTrees[c], jets_[c].BranchMap(isMC));
  }
  if (mode != RunMode::MC) {
    SetBranches(m.skimTree, filters_.BranchMap(cfg.filterBranch));
  }
  if (mode == RunMode::Triggered) {
    SetBranches(m.trigTree, {{cfg.hltJ80Branch, &hltJ80_}});
  }

  // event hists
  m.hvz_all = new TH1D("hvz_all", "all events;v_{z} (cm);N", 40, -20, 20);
  m.hvz = new TH1D("hvz", "after vz+filter;v_{z} (cm);N", 40, -20, 20);
  if (mode != RunMode::MC) {
    m.hfilt = new TH1I("hfilt", "ppvF;filter;N", 2, 0, 2);
  }
  if (mode == RunMode::Triggered) {
    m.h_j80 = new TH1I("h_hlt_j80", "HLT_AK4PFJet80;bit;N", 2, 0, 2);
  }
  for (TH1 *h :
       {(TH1 *)m.hvz_all, (TH1 *)m.hvz, (TH1 *)m.hfilt, (TH1 *)m.h_j80}) {
    if (h) {
      h->SetDirectory(nullptr);
    }
  }

  const Long64_t nEvents = m.evtTree->GetEntries();
  nLoop_ = (maxEvents < 0) ? nEvents : std::min(maxEvents, nEvents);
}

ForestEventLoop::~ForestEventLoop() {
  delete impl_->hvz_all;
  delete impl_->hvz;
  delete impl_->hfilt;
  delete impl_->h_j80;
  if (impl_->fi) {
    impl_->fi->Close();
  }
}

bool ForestEventLoop::Next() {
  while (++entry_ < nLoop_) {
    if (SelectEntry(entry_)) {
      return true;
    }
  }
  return false;
}

bool ForestEventLoop::SelectEntry(Long64_t i) {
  Impl &m = *impl_;
  const size_t nCones = jets_.size();

  // vz, rho
  m.evtTree->GetEntry(i);
  m.ggTree->GetEntry(i);
  m.hvz_all->Fill(event_.vz);
  if (TMath::Abs(event_.vz) > kVzCut) {
    return false;
  }

  // primary vertex filter
  if (mode_ != RunMode::MC) {
    m.skimTree->GetEntry(i);
    m.hfilt->Fill(filters_.ppvF);
    if (filters_.ppvF == 0) {
      return false;
    }
  }
  m.hvz->Fill(event_.vz);

  // jet trigger hist
  if (mode_ == RunMode::Triggered) {
    m.trigTree->GetEntry(i);
    m.h_j80->Fill(hltJ80_);
  }

  // golden JSON
  if (m.dcs && !m.dcs->isGood(event_.run, event_.lumi)) {
    return false;
  }

  // jet trees
  for (size_t c = 0; c < nCones; c++) {
    m.jetTrees[c]->GetEntry(i);
  }
  if (jets_[m.trigConeIdx].reco.nref < 2) {
    return false;
  }

  // applying JEC
  for (size_t c = 0; c < nCones; c++) {
    for (int j = 0; j < jets_[c].reco.nref; j++) {
      m.jecs[c].SetJetPT(jets_[c].reco.rawpt[j]);
      m.jecs[c].SetJetEta(jets_[c].reco.eta[j]);
      m.jecs[c].SetJetPhi(jets_[c].reco.phi[j]);
      m.jecs[c].SetJetArea(jets_[c].reco.area[j]);
      m.jecs[c].SetRho(event_.rho);
      corrPt_[c][j] = (float)m.jecs[c].GetCorrectedPT();
    }
  }

  // JER SF closure: smear MC jets ahead of any downstream jet ordering or
  // selection, mirroring how L2Residual reprocessing sits on corrPt above
  if (jerClosure_) {
    for (size_t c = 0; c < nCones; c++) {
      for (int j = 0; j < jets_[c].reco.nref; j++) {
        JetSmearing::Result sm = JetSmearing::ComputeSmearFactor(
            corrPt_[c][j], jets_[c].reco.eta[j], event_.rho, jets_[c].ref.pt[j],
            m.jerResolution[c], m.jerScaleFactor[c], m.jerRng,
            Variation::NOMINAL, "", 3.0, m.jerMethod);
        corrPt_[c][j] =
            (float)JetSmearing::SmearedPt(corrPt_[c][j], sm.smearFactor);
      }
    }
  }

  // sorted structure of jet indices per cone
  for (size_t c = 0; c < nCones; c++) {
    sorted_[c] = FindLeadingJets(corrPt_[c].data(), jets_[c].reco.nref);
  }

  // cutting on trigger pT threshold near full efficiency
  if (mode_ == RunMode::Triggered) {
    if (hltJ80_ == 0) {
      return false;
    }
    const size_t t = m.trigConeIdx;
    if (hltJ80_ == 1 && corrPt_[t][sorted_[t].lead] <= m.hltJ80Thresh) {
      return false;
    }
  }

  return true;
}

bool ForestEventLoop::GoodJet(size_t c, int j) const {
  const auto &r = jets_[c].reco;
  return impl_->js->JetSelection(r.eta[j], r.phi[j], r.pf.CHF[j], r.pf.NHF[j],
                                 r.pf.CEF[j], r.pf.NEF[j], r.pf.MUF[j],
                                 r.pf.CHM[j], r.pf.NHM[j], r.pf.CEM[j],
                                 r.pf.NEM[j], r.pf.MUM[j]);
}

bool ForestEventLoop::EventVetoed(size_t c) const {
  const auto &r = jets_[c].reco;
  return impl_->js->VetoEvent(r.nref, corrPt_[c].data(), r.eta, r.phi, r.pf.CHF,
                              r.pf.NHF, r.pf.CEF, r.pf.NEF, r.pf.MUF, r.pf.CHM,
                              r.pf.NHM, r.pf.CEM, r.pf.NEM, r.pf.MUM);
}

void ForestEventLoop::WriteEventHists(TDirectory *dir) const {
  TDirectory::TContext ctx(dir);
  impl_->hvz_all->Write();
  impl_->hvz->Write();
  if (impl_->hfilt) {
    impl_->hfilt->Write();
  }
  if (impl_->h_j80) {
    impl_->h_j80->Write();
  }
}
