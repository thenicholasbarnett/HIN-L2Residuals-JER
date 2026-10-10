#include "RunPopulation.h"

#include "TFile.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TString.h"
#include "Rtypes.h"

#include "Binning.h"
#include "ForestEventLoop.h"
#include "Utilities.h"

#include "AnalysisConfig.h"

#include <vector>

void runPopulation(TString input, TString output, TString modeFlag,
                   Long64_t maxEvents) {

  const AnalysisConfig &cfg = Config();
  PrintConfigSummary(cfg);

  // forest, event selection, JEC, pT ordering -- no nref requirement
  ForestEventLoop loop(input, ParseRunMode(modeFlag), maxEvents, false, 0);

  // first configured cone only (ak4PF in cfg/2024ppRef_population.toml)
  const size_t c = 0;

  // log pT bins
  const int nPt = 60;
  const double ptLo = 5.0;
  const double ptHi = 1000.0;
  std::vector<Double_t> ptEdges(nPt + 1);
  for (int i = 0; i <= nPt; i++) {
    ptEdges[i] = ptLo * TMath::Power(ptHi / ptLo, (double)i / nPt);
  }

  const AxisBins flag = {2, 0.0, 2.0, ""};
  THnSparse *h = MakeTHnSparse<THnSparseD>(
      cfg.coneLabels[c] + "_population", "",
      {{104, -5.2, 5.2, "#eta"},
       {72, -(Float_t)TMath::Pi(), (Float_t)TMath::Pi(), "#phi (rad)"},
       {nPt, (Float_t)ptLo, (Float_t)ptHi, "corrected p_{T} (GeV)"},
       flag,
       flag,
       flag,
       flag});
  h->GetAxis(2)->Set(nPt, ptEdges.data());
  h->GetAxis(3)->SetTitle("passes jet ID");
  h->GetAxis(4)->SetTitle("outside veto map");
  h->GetAxis(5)->SetTitle("event not vetoed");
  h->GetAxis(6)->SetTitle("leading jet");
  h->Sumw2();

  // event loop
  while (loop.Next()) {

    const float weight = loop.Weight();
    const auto &reco = loop.JetsIn(c).reco;
    const float *corrPt = loop.CorrPt(c);
    const int lead = loop.Sorted(c).lead;
    const bool eventPass = !loop.EventVetoed(c);

    for (int j = 0; j < reco.nref; j++) {
      if (corrPt[j] < cfg.minJetPt) {
        continue;
      }
      const Double_t x[7] = {reco.eta[j],
                             reco.phi[j],
                             corrPt[j],
                             (Double_t)loop.PassesJetID(c, j),
                             (Double_t)!loop.InVetoRegion(c, j),
                             (Double_t)eventPass,
                             (Double_t)(j == lead)};
      h->Fill(x, weight);
    }
  }

  // output
  TFile *fo = CreateOutputFile(output);
  fo->cd();
  h->Write();
  fo->Close();
}
