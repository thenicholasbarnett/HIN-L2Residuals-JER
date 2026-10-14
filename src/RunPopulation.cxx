#include "RunPopulation.h"

#include "TFile.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TString.h"
#include "Rtypes.h"

#include "Binning.h"
#include "ForestEventLoop.h"
#include "JetSmearer.h"
#include "Utilities.h"

#include "AnalysisConfig.h"

#include <vector>

namespace {

// cone of the population collections (ak4PF, its gen jets, AK4 calo jets)
constexpr double kConeR = 0.4;
constexpr float kConeArea = (float)(TMath::Pi() * kConeR * kConeR);
// gen and calo jets take the jet ID decision of the PF jet within
// kMatchDRFraction * R, closest pairs first, no PF jet used twice
constexpr double kMatchDRFraction = 0.5;

// (eta, phi, pT, passes jet ID, outside veto map, event not vetoed, leading
// jet, matched PF jet), same layout for every collection; PF jets are their
// own match
THnSparse *BookPopulation(const TString &name, const char *ptTitle,
                          const std::vector<Double_t> &ptEdges) {
  const int nPt = (int)ptEdges.size() - 1;
  const AxisBins flag = {2, 0.0, 2.0, ""};
  THnSparse *h = MakeTHnSparse<THnSparseD>(
      name, "",
      {{104, -5.2, 5.2, "#eta"},
       {72, -(Float_t)TMath::Pi(), (Float_t)TMath::Pi(), "#phi (rad)"},
       {nPt, (Float_t)ptEdges.front(), (Float_t)ptEdges.back(), ptTitle},
       flag,
       flag,
       flag,
       flag,
       flag});
  h->GetAxis(2)->Set(nPt, ptEdges.data());
  h->GetAxis(3)->SetTitle("passes jet ID");
  h->GetAxis(4)->SetTitle("outside veto map");
  h->GetAxis(5)->SetTitle("event not vetoed");
  h->GetAxis(6)->SetTitle("leading jet");
  h->GetAxis(7)->SetTitle("matched PF jet");
  h->Sumw2();
  return h;
}

// PF jet pT after each correction stage: (eta, raw, forest jtpt, L1, L1 +
// L2Relative, final), final = + L2Residual in data, + JER smearing in MC,
// whichever the config provides
THnSparse *BookStages(const TString &name, const std::vector<Double_t> &ptEdges,
                      bool isMC) {
  const int nPt = (int)ptEdges.size() - 1;
  const AxisBins pt = {nPt, (Float_t)ptEdges.front(), (Float_t)ptEdges.back(),
                       ""};
  THnSparse *h = MakeTHnSparse<THnSparseD>(
      name, "", {{104, -5.2, 5.2, "#eta"}, pt, pt, pt, pt, pt});
  const char *titles[5] = {"raw p_{T} (GeV)", "forest p_{T} (GeV)",
                           "L1 p_{T} (GeV)", "L1 + L2Relative p_{T} (GeV)",
                           isMC ? "L1 + L2Relative + smeared p_{T} (GeV)"
                                : "L1 + L2Relative + L2Residual p_{T} (GeV)"};
  for (int a = 0; a < 5; a++) {
    h->GetAxis(a + 1)->Set(nPt, ptEdges.data());
    h->GetAxis(a + 1)->SetTitle(titles[a]);
  }
  h->Sumw2();
  return h;
}

// gen or calo jets: position and pT of the jet itself, jet ID from the
// matched PF jet (fails when unmatched)
void FillMatched(THnSparse *h, const ForestEventLoop &loop, size_t c, int n,
                 const Float_t *pt, const Float_t *eta, const Float_t *phi,
                 double minPt, bool eventPass, float weight) {
  const auto &reco = loop.JetsIn(c).reco;
  const std::vector<int> match =
      JetSmearing::Match(n, eta, phi, reco.nref, reco.eta, reco.phi, kConeR,
                         JetSmearing::Mode::OneToOne, kMatchDRFraction);
  int lead = -1;
  for (int j = 0; j < n; j++) {
    if (lead < 0 || pt[j] > pt[lead]) {
      lead = j;
    }
  }
  for (int j = 0; j < n; j++) {
    if (pt[j] < minPt) {
      continue;
    }
    const bool matched = match[j] >= 0;
    const Double_t x[8] = {eta[j],
                           phi[j],
                           pt[j],
                           (Double_t)(matched && loop.PassesJetID(c, match[j])),
                           (Double_t)!loop.InVetoRegion(eta[j], phi[j]),
                           (Double_t)eventPass,
                           (Double_t)(j == lead),
                           (Double_t)matched};
    h->Fill(x, weight);
  }
}

} // namespace

void runPopulation(TString input, TString output, TString modeFlag,
                   Long64_t maxEvents) {

  const AnalysisConfig &cfg = Config();
  PrintConfigSummary(cfg);

  // forest, event selection, JEC, pT ordering -- no nref requirement
  // calo jets read from the PF jet tree
  // MC is smeared when the config has JER files (final stage of the stages
  // sparse, and the corrected pT everywhere else)
  const RunMode mode = ParseRunMode(modeFlag);
  const bool smear =
      mode == RunMode::MC && !cfg.jerScaleFactorFilesPerCone.empty();
  ForestEventLoop loop(input, mode, maxEvents, smear, 0, true);

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

  const TString &cone = cfg.coneLabels[c];
  TString caloName = cone;
  caloName.ReplaceAll("PF", "Calo");
  TString genName = cone;
  genName.ReplaceAll("PF", "Gen");
  THnSparse *h =
      BookPopulation(cone + "_population", "corrected p_{T} (GeV)", ptEdges);
  THnSparse *hCalo = BookPopulation(caloName + "_population",
                                    "corrected calo p_{T} (GeV)", ptEdges);
  THnSparse *hGen = loop.IsMC() ? BookPopulation(genName + "_population",
                                                 "gen p_{T} (GeV)", ptEdges)
                                : nullptr;

  THnSparse *hStages = BookStages(cone + "_stages", ptEdges, loop.IsMC());

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
      const Double_t x[8] = {reco.eta[j],
                             reco.phi[j],
                             corrPt[j],
                             (Double_t)loop.PassesJetID(c, j),
                             (Double_t)!loop.InVetoRegion(c, j),
                             (Double_t)eventPass,
                             (Double_t)(j == lead),
                             1.0};
      h->Fill(x, weight);

      const Double_t xs[6] = {
          reco.eta[j],           reco.rawpt[j],         reco.pt[j],
          loop.StagePt(c, j, 1), loop.StagePt(c, j, 2), corrPt[j]};
      hStages->Fill(xs, weight);
    }

    // calo jets: the PF chain of the same cone on the stored (raw) calo pT,
    // nominal area since the forest has none for them
    const auto &calo = loop.JetsIn(c).calo;
    std::vector<Float_t> caloPt(calo.n);
    for (int j = 0; j < calo.n; j++) {
      caloPt[j] =
          loop.CorrectPt(c, calo.pt[j], calo.eta[j], calo.phi[j], kConeArea);
    }
    FillMatched(hCalo, loop, c, calo.n, caloPt.data(), calo.eta, calo.phi,
                cfg.minJetPt, eventPass, weight);
    if (hGen) {
      const auto &gen = loop.JetsIn(c).gen;
      FillMatched(hGen, loop, c, gen.n, gen.pt, gen.eta, gen.phi, cfg.minJetPt,
                  eventPass, weight);
    }
  }

  // output
  TFile *fo = CreateOutputFile(output);
  fo->cd();
  h->Write();
  hStages->Write();
  hCalo->Write();
  if (hGen) {
    hGen->Write();
  }
  fo->Close();
}
