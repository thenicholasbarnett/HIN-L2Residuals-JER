#include "RunAsymmetry.h"

#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TString.h"
#include "Rtypes.h"

#include "Binning.h"
#include "Dijet.h"
#include "DijetHistograms.h"
#include "ForestEventLoop.h"

#include "AnalysisConfig.h"

#include <vector>

void runAsymmetry(TString input, TString output, TString modeFlag,
                  Long64_t maxEvents, bool jerClosure) {

  const AnalysisConfig &cfg = Config();
  PrintConfigSummary(cfg);

  // forest, event selection, JEC (+ JER smearing), pT ordering
  ForestEventLoop loop(input, ParseRunMode(modeFlag), maxEvents, jerClosure);
  const bool isMC = loop.IsMC();
  const size_t nCones = loop.NCones();

  // JME Run 3 event veto on the configured cone's jets (AK4, as the maps
  // are derived), whole event dropped for every cone; "" = off
  int vetoCone = -1;
  for (size_t c = 0; c < nCones; c++) {
    if (cfg.coneLabels[c] == cfg.jetEventVetoCone.c_str()) {
      vetoCone = (int)c;
    }
  }
  TH1I *hEventVeto = nullptr;
  if (vetoCone >= 0) {
    hEventVeto = new TH1I("h_eventveto", "JME event veto;vetoed;N", 2, 0, 2);
    hEventVeto->SetDirectory(nullptr);
  }

  // jet hists
  BinningConfig bins;
  std::vector<ConeHistograms> cones(nCones);
  for (size_t c = 0; c < nCones; c++) {
    cones[c].Init(cfg.coneLabels[c], bins, isMC);
  }

  // event loop
  while (loop.Next()) {

    const EventStruct &event = loop.Event();
    const float weight = loop.Weight();

    if (hEventVeto) {
      const bool vetoed = loop.EventVetoed(vetoCone);
      hEventVeto->Fill(vetoed);
      if (vetoed) {
        continue;
      }
    }

    // jet ID -> dijet -> A -> fill
    for (size_t c = 0; c < nCones; c++) {
      const auto &reco = loop.JetsIn(c).reco;
      const auto &ref = loop.JetsIn(c).ref;
      const float *corrPt = loop.CorrPt(c);
      const SortedJets &sorted = loop.Sorted(c);

      // inclusive jets
      for (int j = 0; j < reco.nref; j++) {
        if (corrPt[j] >= cfg.minJetPt) {
          cones[c].FillInclJet(corrPt[j], reco.eta[j], reco.phi[j], weight);
          if (isMC) {
            cones[c].FillInclJetResp(corrPt[j], reco.rawpt[j], reco.pt[j],
                                     reco.eta[j], ref.pt[j], event.rho, weight);
          }

          // same tight jet ID + veto map already applied to lead/sublead/
          // third below, but with no dijet/trigger-jet-pairing requirement
          // -- every jet in the event gets checked independently
          if (loop.GoodJet(c, j)) {
            cones[c].FillInclJetID(corrPt[j], reco.eta[j], reco.phi[j], weight);
          }
        }
      }

      if (sorted.sublead == -1) {
        continue;
      }

      // jet ID & veto map applied to all jets used in this analysis
      if (!loop.GoodJet(c, sorted.lead) || !loop.GoodJet(c, sorted.sublead)) {
        continue;
      }
      const bool hasThird = sorted.third != -1 && loop.GoodJet(c, sorted.third);

      // dijet logic (see include/Dijet.h)
      DijetResult dijet =
          MakeDijet(sorted, hasThird, corrPt, reco.eta, reco.phi, event.event,
                    cfg.minJetPt, cfg.dphiCut);
      if (!dijet.valid) {
        continue;
      }

      // |A| cut
      if (TMath::Abs(dijet.A) > cfg.maxAbsA) {
        continue;
      }

      // fill histograms
      cones[c].Fill(dijet, corrPt, reco.eta, reco.phi, weight);
      if (isMC) {
        cones[c].FillResp(dijet, corrPt, reco.rawpt, reco.pt, reco.eta, ref.pt,
                          event.rho, weight);
      }
    }
  }

  // output
  TFile *fo = CreateOutputFile(output);
  loop.WriteEventHists(fo);
  if (hEventVeto) {
    fo->cd();
    hEventVeto->Write();
  }
  for (size_t c = 0; c < nCones; c++) {
    TDirectory *dir = fo->mkdir(cfg.coneLabels[c].Data());
    cones[c].Write(dir);
    fo->cd();
  }
  fo->Close();
}
