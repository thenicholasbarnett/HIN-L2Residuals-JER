#include "RunPopulation.h"

#include "TFile.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TString.h"
#include "Rtypes.h"

#include "Binning.h"
#include "ForestEventLoop.h"
#include "JetSelector.h"
#include "JetSmearer.h"
#include "Utilities.h"

#include "AnalysisConfig.h"

#include <algorithm>
#include <vector>

namespace {

// cone of the population collections (ak4PF, its gen jets, AK4 calo jets)
constexpr double kConeR = 0.4;
constexpr float kConeArea = (float)(TMath::Pi() * kConeR * kConeR);
// gen and calo jets take the jet ID decision of the PF jet within
// kMatchDRFraction * R, closest pairs first, no PF jet used twice; the gen
// match of PF and calo jets uses the same definition
constexpr double kMatchDRFraction = 0.5;
// MC pT-hat spike cut factor, as DiJetJERC's (runDijetProfiles)
constexpr double kSpikeFactor = 1.5;

// the two Run 3 veto map types, every population sparse carries both
const char *const kMapTypes[2] = {"jetvetomap", "jetvetomap_all"};
// best-guess selection (stages sparses): jetvetomap_all, event veto at
// p_T^veto JetSelector::kEventVetoMinPt (15 GeV)
constexpr int kBestMap = 1;

// one JetSelector per map type, same jet ID JSON as the event loop; the ID is
// taken from the loop, these only evaluate the maps
struct VetoMaps {
  std::vector<JetSelector> js;

  explicit VetoMaps(const AnalysisConfig &cfg) {
    for (const char *type : kMapTypes) {
      js.emplace_back(cfg.jetSystem == "ion" ? JetSelector::System::Ion
                                             : JetSelector::System::pp,
                      cfg.jetIdPath, cfg.vetoMapPath);
      js.back().SetVetoMap(type);
    }
  }

  bool Outside(int m, double eta, double phi) const {
    return !js[m].InVetoRegion(eta, phi);
  }
};

// p_T^veto per map: highest corrected pT among the event's jets that would
// veto it (jet ID, EMF < 0.9, inside the map), 0 when there is none. The
// event passes a veto threshold T when p_T^veto <= T (JetSelector::VetoEvent,
// T = 15 GeV)
std::vector<double> VetoPt(const ForestEventLoop &loop, size_t c,
                           const VetoMaps &maps) {
  const auto &reco = loop.JetsIn(c).reco;
  const float *corrPt = loop.CorrPt(c);
  std::vector<double> out(maps.js.size(), 0.0);
  for (int j = 0; j < reco.nref; j++) {
    if (reco.pf.CEF[j] + reco.pf.NEF[j] >= JetSelector::kEventVetoMaxEMF) {
      continue;
    }
    if (!loop.PassesJetID(c, j)) {
      continue;
    }
    for (size_t m = 0; m < maps.js.size(); m++) {
      if (corrPt[j] > out[m] && !maps.Outside(m, reco.eta[j], reco.phi[j])) {
        out[m] = corrPt[j];
      }
    }
  }
  return out;
}

bool PassesEventVeto(const std::vector<double> &vetoPt) {
  return vetoPt[kBestMap] <= JetSelector::kEventVetoMinPt;
}

// the best-guess selection, as a title for the plot legends
const char *kBestSelection =
    "jet ID + jetvetomap_all + p_{T}^{veto} < 15 GeV";

// p_T^veto axis: the population pT edges plus 15 GeV, so the default event
// veto is a bin edge; 0 (no vetoing jet) lands in the underflow
std::vector<Double_t> VetoPtEdges(const std::vector<Double_t> &ptEdges) {
  std::vector<Double_t> e = ptEdges;
  const double cut = JetSelector::kEventVetoMinPt;
  if (std::find(e.begin(), e.end(), cut) == e.end()) {
    e.insert(std::upper_bound(e.begin(), e.end(), cut), cut);
  }
  return e;
}

// MC pT-hat spikes: DiJetJERC's cut (leading gen jet < 1.5 pthat or leading
// jet < 1.5 leading gen jet), plus a leading jet without a forest gen match
// above 1.5 pthat -- pileup overlay in a low pthat event, which JME drops with
// Pileup_pthatmax > pthat (not in the forest)
bool IsSpike(const ForestEventLoop &loop, size_t c) {
  const int lead = loop.Sorted(c).lead;
  if (lead < 0) {
    return false;
  }
  const auto &jets = loop.JetsIn(c);
  const double pthat = loop.Event().pthat;
  const double pt1 = loop.CorrPt(c)[lead];
  double gen1 = 0.0;
  for (int j = 0; j < jets.gen.n; j++) {
    gen1 = TMath::Max(gen1, (double)jets.gen.pt[j]);
  }
  if (!(gen1 < kSpikeFactor * pthat || pt1 < kSpikeFactor * gen1)) {
    return true;
  }
  return jets.ref.pt[lead] < 0 && pt1 > kSpikeFactor * pthat;
}

// leading and subleading corrected jets (no jet ID) against the config's
// [population] cuts
bool PassesDijetCuts(const ForestEventLoop &loop, size_t c,
                     const AnalysisConfig &cfg) {
  const auto &s = loop.Sorted(c);
  const auto &reco = loop.JetsIn(c).reco;
  const float *pt = loop.CorrPt(c);
  const bool cutLead = cfg.popLeadPtMin > 0 || cfg.popLeadAbsEtaMax >= 0;
  const bool cutSublead =
      cfg.popSubleadPtMin > 0 || cfg.popSubleadAbsEtaMax >= 0;
  if (cutLead) {
    if (s.lead < 0 || pt[s.lead] <= cfg.popLeadPtMin) {
      return false;
    }
    if (cfg.popLeadAbsEtaMax >= 0 &&
        TMath::Abs(reco.eta[s.lead]) >= cfg.popLeadAbsEtaMax) {
      return false;
    }
  }
  if (cutSublead) {
    if (s.sublead < 0 || pt[s.sublead] <= cfg.popSubleadPtMin) {
      return false;
    }
    if (cfg.popSubleadAbsEtaMax >= 0 &&
        TMath::Abs(reco.eta[s.sublead]) >= cfg.popSubleadAbsEtaMax) {
      return false;
    }
  }
  return true;
}

// per jet: matched to a gen jet (MC), all false in data
std::vector<bool> GenMatched(const ForestEventLoop &loop, size_t c, int n,
                             const Float_t *eta, const Float_t *phi) {
  std::vector<bool> out(n, false);
  if (!loop.IsMC()) {
    return out;
  }
  const auto &gen = loop.JetsIn(c).gen;
  const std::vector<int> match =
      JetSmearing::Match(n, eta, phi, gen.n, gen.eta, gen.phi, kConeR,
                         JetSmearing::Mode::OneToOne, kMatchDRFraction);
  for (int j = 0; j < n; j++) {
    out[j] = match[j] >= 0;
  }
  return out;
}

// (eta, phi, pT, passes jet ID, outside jetvetomap, outside jetvetomap_all,
// p_T^veto jetvetomap, p_T^veto jetvetomap_all, leading jet, matched PF jet,
// matched gen jet), same layout for every collection; PF and gen jets are
// their own matches, nothing is gen matched in data. pT is the fully
// corrected one: the selection effect at the best-guess calibration
constexpr int kNPopAxes = 11;

THnSparse *BookPopulation(const TString &name, const char *ptTitle,
                          const std::vector<Double_t> &ptEdges) {
  const int nPt = (int)ptEdges.size() - 1;
  const std::vector<Double_t> vEdges = VetoPtEdges(ptEdges);
  const int nV = (int)vEdges.size() - 1;
  const AxisBins flag = {2, 0.0, 2.0, ""};
  const AxisBins veto = {nV, (Float_t)vEdges.front(), (Float_t)vEdges.back(),
                         ""};
  THnSparse *h = MakeTHnSparse<THnSparseD>(
      name, "",
      {{104, -5.2, 5.2, "#eta"},
       {72, -(Float_t)TMath::Pi(), (Float_t)TMath::Pi(), "#phi (rad)"},
       {nPt, (Float_t)ptEdges.front(), (Float_t)ptEdges.back(), ptTitle},
       flag,
       flag,
       flag,
       veto,
       veto,
       flag,
       flag,
       flag});
  h->GetAxis(2)->Set(nPt, ptEdges.data());
  h->GetAxis(3)->SetTitle("passes jet ID");
  h->GetAxis(4)->SetTitle("outside jetvetomap");
  h->GetAxis(5)->SetTitle("outside jetvetomap_all");
  h->GetAxis(6)->Set(nV, vEdges.data());
  h->GetAxis(6)->SetTitle("p_{T}^{veto} jetvetomap (GeV)");
  h->GetAxis(7)->Set(nV, vEdges.data());
  h->GetAxis(7)->SetTitle("p_{T}^{veto} jetvetomap_all (GeV)");
  h->GetAxis(8)->SetTitle("leading jet");
  h->GetAxis(9)->SetTitle("matched PF jet");
  h->GetAxis(10)->SetTitle("matched gen jet");
  h->Sumw2();
  return h;
}

// PF jet pT after each correction stage: (eta, raw, forest jtpt, L1, L1 +
// L2Relative, final), final = + L2Residual in data, + JER smearing in MC,
// whichever the config provides. Best-guess selection only (kBestSelection):
// the calibration effect at the best-guess selection
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

// MC response after each correction stage, for JES (mean) and JER (width):
// (|eta reco|, gen pT, raw / gen, forest jtpt / gen, L1 / gen, L1 +
// L2Relative / gen, final / gen), final = + JER smearing when configured.
// Gen-matched PF jets passing the best-guess selection; gen pT from the
// forest's own match (ref), the one the smearing uses. Binned in gen pT, not
// reco, against the falling-spectrum migration bias
const int kNRespPt = 40;

std::vector<Double_t> RespPtEdges() {
  const double ptLo = 5.0;
  const double ptHi = 1000.0;
  std::vector<Double_t> e(kNRespPt + 1);
  for (int i = 0; i <= kNRespPt; i++) {
    e[i] = ptLo * TMath::Power(ptHi / ptLo, (double)i / kNRespPt);
  }
  return e;
}

THnSparse *BookResponseStages(const TString &name, bool smeared) {
  const int nPt = kNRespPt;
  const std::vector<Double_t> ptEdges = RespPtEdges();
  const AxisBins resp = {100, 0.0, 2.0, ""};
  THnSparse *h = MakeTHnSparse<THnSparseD>(
      name, "",
      {{52, 0.0, 5.2, "|#eta^{reco}|"},
       {nPt, (Float_t)ptEdges.front(), (Float_t)ptEdges.back(),
        "p_{T}^{gen} (GeV)"},
       resp,
       resp,
       resp,
       resp,
       resp});
  h->GetAxis(1)->Set(nPt, ptEdges.data());
  const char *titles[5] = {
      "raw p_{T} / p_{T}^{gen}", "forest p_{T} / p_{T}^{gen}",
      "L1 p_{T} / p_{T}^{gen}", "L1 + L2Relative p_{T} / p_{T}^{gen}",
      smeared ? "L1 + L2Relative + smeared p_{T} / p_{T}^{gen}"
              : "L1 + L2Relative p_{T} / p_{T}^{gen}"};
  for (int a = 0; a < 5; a++) {
    h->GetAxis(a + 2)->SetTitle(titles[a]);
  }
  h->Sumw2();
  return h;
}

// MC final response against the selection, for JES and JER per selection:
// (|eta reco|, gen pT, final / gen, passes jet ID, outside jetvetomap,
// outside jetvetomap_all, p_T^veto jetvetomap, p_T^veto jetvetomap_all).
// Every gen-matched PF jet: the selection effect at the best-guess
// calibration
THnSparse *BookResponseSelection(const TString &name, bool smeared,
                                 const std::vector<Double_t> &vEdges) {
  const std::vector<Double_t> ptEdges = RespPtEdges();
  const int nV = (int)vEdges.size() - 1;
  const AxisBins flag = {2, 0.0, 2.0, ""};
  const AxisBins veto = {nV, (Float_t)vEdges.front(), (Float_t)vEdges.back(),
                         ""};
  THnSparse *h = MakeTHnSparse<THnSparseD>(
      name, "",
      {{52, 0.0, 5.2, "|#eta^{reco}|"},
       {kNRespPt, (Float_t)ptEdges.front(), (Float_t)ptEdges.back(),
        "p_{T}^{gen} (GeV)"},
       {100, 0.0, 2.0,
        smeared ? "L1 + L2Relative + smeared p_{T} / p_{T}^{gen}"
                : "L1 + L2Relative p_{T} / p_{T}^{gen}"},
       flag,
       flag,
       flag,
       veto,
       veto});
  h->GetAxis(1)->Set(kNRespPt, ptEdges.data());
  h->GetAxis(3)->SetTitle("passes jet ID");
  h->GetAxis(4)->SetTitle("outside jetvetomap");
  h->GetAxis(5)->SetTitle("outside jetvetomap_all");
  h->GetAxis(6)->Set(nV, vEdges.data());
  h->GetAxis(6)->SetTitle("p_{T}^{veto} jetvetomap (GeV)");
  h->GetAxis(7)->Set(nV, vEdges.data());
  h->GetAxis(7)->SetTitle("p_{T}^{veto} jetvetomap_all (GeV)");
  h->Sumw2();
  return h;
}

// gen or calo jets: position and pT of the jet itself, jet ID from the
// matched PF jet (fails when unmatched); isGen: the jets are gen jets
void FillMatched(THnSparse *h, const ForestEventLoop &loop, size_t c, int n,
                 const Float_t *pt, const Float_t *eta, const Float_t *phi,
                 double minPt, const VetoMaps &maps,
                 const std::vector<double> &vetoPt, float weight, bool isGen) {
  const auto &reco = loop.JetsIn(c).reco;
  const std::vector<int> match =
      JetSmearing::Match(n, eta, phi, reco.nref, reco.eta, reco.phi, kConeR,
                         JetSmearing::Mode::OneToOne, kMatchDRFraction);
  const std::vector<bool> genMatch =
      isGen ? std::vector<bool>(n, true) : GenMatched(loop, c, n, eta, phi);
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
    const Double_t x[kNPopAxes] = {
        eta[j],
        phi[j],
        pt[j],
        (Double_t)(matched && loop.PassesJetID(c, match[j])),
        (Double_t)maps.Outside(0, eta[j], phi[j]),
        (Double_t)maps.Outside(1, eta[j], phi[j]),
        vetoPt[0],
        vetoPt[1],
        (Double_t)(j == lead),
        (Double_t)matched,
        (Double_t)genMatch[j]};
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
  const VetoMaps maps(cfg);

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
  // the jet selection, for the plot legends
  hStages->SetTitle(kBestSelection);
  THnSparse *hResp = loop.IsMC()
                         ? BookResponseStages(cone + "_response_stages", smear)
                         : nullptr;
  if (hResp) {
    hResp->SetTitle(kBestSelection);
  }
  THnSparse *hRespSel =
      loop.IsMC() ? BookResponseSelection(cone + "_response_selection", smear,
                                          VetoPtEdges(ptEdges))
                  : nullptr;

  // event loop
  while (loop.Next()) {

    if (loop.IsMC() && IsSpike(loop, c)) {
      continue;
    }
    if (!PassesDijetCuts(loop, c, cfg)) {
      continue;
    }

    const float weight = loop.Weight();
    const auto &reco = loop.JetsIn(c).reco;
    const float *corrPt = loop.CorrPt(c);
    const int lead = loop.Sorted(c).lead;
    const std::vector<double> vetoPt = VetoPt(loop, c, maps);
    const bool eventPass = PassesEventVeto(vetoPt);
    const std::vector<bool> genMatch =
        GenMatched(loop, c, reco.nref, reco.eta, reco.phi);

    for (int j = 0; j < reco.nref; j++) {
      if (corrPt[j] < cfg.minJetPt) {
        continue;
      }
      const bool id = loop.PassesJetID(c, j);
      const bool outside[2] = {maps.Outside(0, reco.eta[j], reco.phi[j]),
                               maps.Outside(1, reco.eta[j], reco.phi[j])};
      const Double_t x[kNPopAxes] = {reco.eta[j],
                                     reco.phi[j],
                                     corrPt[j],
                                     (Double_t)id,
                                     (Double_t)outside[0],
                                     (Double_t)outside[1],
                                     vetoPt[0],
                                     vetoPt[1],
                                     (Double_t)(j == lead),
                                     1.0,
                                     (Double_t)genMatch[j]};
      h->Fill(x, weight);

      // gen-matched jets (forest ref, refpt < 0 = no match)
      const auto &ref = loop.JetsIn(c).ref;
      const bool best = id && outside[kBestMap] && eventPass;
      if (hRespSel && ref.pt[j] > 0) {
        const double g = ref.pt[j];
        const Double_t xs[8] = {TMath::Abs(reco.eta[j]),
                                g,
                                corrPt[j] / g,
                                (Double_t)id,
                                (Double_t)outside[0],
                                (Double_t)outside[1],
                                vetoPt[0],
                                vetoPt[1]};
        hRespSel->Fill(xs, weight);
      }

      if (!best) {
        continue;
      }
      if (hResp && ref.pt[j] > 0) {
        const double g = ref.pt[j];
        const Double_t xr[7] = {TMath::Abs(reco.eta[j]),
                                g,
                                reco.rawpt[j] / g,
                                reco.pt[j] / g,
                                loop.StagePt(c, j, 1) / g,
                                loop.StagePt(c, j, 2) / g,
                                corrPt[j] / g};
        hResp->Fill(xr, weight);
      }
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
                cfg.minJetPt, maps, vetoPt, weight, false);
    if (hGen) {
      const auto &gen = loop.JetsIn(c).gen;
      FillMatched(hGen, loop, c, gen.n, gen.pt, gen.eta, gen.phi, cfg.minJetPt,
                  maps, vetoPt, weight, true);
    }
  }

  // output
  TFile *fo = CreateOutputFile(output);
  fo->cd();
  h->Write();
  hStages->Write();
  if (hResp) {
    hResp->Write();
  }
  if (hRespSel) {
    hRespSel->Write();
  }
  hCalo->Write();
  if (hGen) {
    hGen->Write();
  }
  fo->Close();
}
