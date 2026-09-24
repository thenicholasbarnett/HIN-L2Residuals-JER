#include "RunDijetProfiles.h"

#include "TFile.h"
#include "TString.h"
#include "Rtypes.h"

#include "Binning.h"
#include "DiJetJERCHistograms.h"
#include "Dijet.h"
#include "DijetProfiles.h"
#include "ForestEventLoop.h"

#include "AnalysisConfig.h"

#include <algorithm>
#include <cmath>
#include <vector>

// Whether this dataset fills the bin holding pt: the triggered sample owns
// bins whose low edge is at or above the plateau threshold, the
// non-triggered one everything below (same rule as Step 3's merge), MC all.
// hadd of the triggered + non-triggered outputs is then the combined file.
static bool OwnsBin(RunMode mode, double pt, double thresh) {
  if (mode == RunMode::MC) {
    return true;
  }
  const auto &v = dijet2::PtEdges();
  int i = 0;
  while (i + 1 < (int)v.size() && pt >= v[i + 1]) {
    i++;
  }
  const double lo = (pt < v.front()) ? 0.0 : v[i];
  return (mode == RunMode::Triggered) == (lo >= thresh);
}

// DiJetJERC event selection + fills for one cone: jet cleaner (ID + veto,
// pT > 10), leading two cleaned jets, dphi > 2.7, MC pT-hat spike cut,
// barrel/probe by |eta| < 1.131 (coin flip if both or neither), alpha from
// the third cleaned jet if > 15 GeV (else 0 without one, 1 with one); in MC
// both jets gen-matched with back-to-back gen jets, else the event is
// dropped entirely, as their selector does. Coin: eventNumber % 2 instead of
// their run/lumi/event-seeded rand(), same role.
static void FillDiJetJERC(const ForestEventLoop &loop, size_t c,
                          DiJetJERCHistograms &h, float thresh) {
  namespace dj = dijetjerc;
  const auto &jets = loop.JetsIn(c);
  const auto &reco = jets.reco;
  const float *pt = loop.CorrPt(c);
  const bool isMC = loop.IsMC();
  const double w = loop.Weight();

  std::vector<int> good;
  for (int j = 0; j < reco.nref; j++) {
    if (pt[j] > dj::kJetPtMin && std::fabs(reco.eta[j]) < 5.2 &&
        loop.GoodJet(c, j)) {
      good.push_back(j);
    }
  }
  if (good.size() < 2) {
    return;
  }
  std::sort(good.begin(), good.end(),
            [&](int a, int b) { return pt[a] > pt[b]; });
  const int j1 = good[0], j2 = good[1];
  if (DijetDPhi(reco.phi[j1], reco.phi[j2]) < dj::kDphiMin) {
    return;
  }

  // gen jets above the gen cleaner threshold (pT ordered in the forest)
  int nGen = 0;
  if (isMC) {
    while (nGen < jets.gen.n && jets.gen.pt[nGen] > dj::kGenJetPtMin) {
      nGen++;
    }
    if (nGen < 2) {
      return;
    }
    const double pthat = loop.Event().pthat;
    if (!(jets.gen.pt[0] < 1.5 * pthat || pt[j1] < 1.5 * jets.gen.pt[0])) {
      return;
    }
  }

  const double ptave = 0.5 * (pt[j1] + pt[j2]);
  // trigger split on this pT_ave (Step 1 pT_avg slices align with it)
  const RunMode mode = loop.Mode();
  const bool ownsReco =
      mode == RunMode::MC || (mode == RunMode::Triggered) == (ptave >= thresh);
  double alpha = (good.size() < 3) ? 0.0 : 1.0;
  if (good.size() > 2 && pt[good[2]] > dj::kJet3PtMin) {
    alpha = pt[good[2]] / ptave;
  }

  const bool c1 = std::fabs(reco.eta[j1]) < dj::kBarrelEta;
  const bool c2 = std::fabs(reco.eta[j2]) < dj::kBarrelEta;
  int b = j1, p = j2;
  if (!c1 && c2) {
    b = j2;
    p = j1;
  } else if (c1 == c2 && loop.Event().event % 2 == 1) {
    b = j2;
    p = j1;
  }
  if (pt[b] < dj::kJet3PtMin && pt[p] < dj::kJet3PtMin) {
    return;
  }
  const int binB = dj::AbsEtaBin(reco.eta[b]);
  const int binP = dj::AbsEtaBin(reco.eta[p]);
  const bool isSM = binB >= 0 && binB == binP;
  const double A = (pt[b] - pt[p]) / (pt[b] + pt[p]);

  // gen side: matched gen jets of barrel and probe, back to back
  const auto &ref = jets.ref;
  if (isMC) {
    if (ref.pt[b] < dj::kGenThr && ref.pt[p] < dj::kGenThr) {
      return;
    }
    if (ref.pt[b] <= 0 || ref.pt[p] <= 0 ||
        DijetDPhi(ref.phi[b], ref.phi[p]) < dj::kDphiMin - 0.1) {
      return;
    }
  }

  // FE fills once per jet sitting in the barrel reference region: the
  // other jet's |eta| bin is filled
  auto fillCat = [&](THnSparse *const *hs, double etaB, double etaP, bool sm,
                     double ptX, double alphaX, double value) {
    if (sm) {
      DiJetJERCHistograms::Fill(hs[DiJetJERCHistograms::kSM], etaP, ptX, alphaX,
                                value, w);
      return;
    }
    if (std::fabs(etaB) < dj::kBarrelEta && dj::AbsEtaBin(etaP) >= 0) {
      DiJetJERCHistograms::Fill(hs[DiJetJERCHistograms::kFE], etaP, ptX, alphaX,
                                value, w);
    }
    if (std::fabs(etaP) < dj::kBarrelEta && dj::AbsEtaBin(etaB) >= 0) {
      DiJetJERCHistograms::Fill(hs[DiJetJERCHistograms::kFE], etaB, ptX, alphaX,
                                value, w);
    }
  };

  if (ownsReco && A != 0.0 && std::fabs(w / A) <= 5e6) {
    fillCat(h.reco, reco.eta[b], reco.eta[p], isSM, ptave, alpha, A);
    if (isMC) {
      auto dR = [&](int j) {
        return std::hypot(reco.eta[j] - ref.eta[j],
                          DijetDPhi(reco.phi[j], ref.phi[j]));
      };
      if (dR(p) < dj::kMatchDR) {
        fillCat(h.mctruth, reco.eta[b], reco.eta[p], isSM, ptave, alpha,
                pt[p] / ref.pt[p]);
      }
      if (dR(b) < dj::kMatchDR) {
        fillCat(h.mctruth, reco.eta[b], reco.eta[p], isSM, ptave, alpha,
                pt[b] / ref.pt[b]);
      }
    }
  }

  if (isMC) {
    const double genPtave = 0.5 * (jets.gen.pt[0] + jets.gen.pt[1]);
    double genAlpha = (nGen < 3) ? 0.0 : 1.0;
    if (nGen > 2 && jets.gen.pt[2] > dj::kGenThr) {
      genAlpha = jets.gen.pt[2] / ptave; // reco pT_ave, as upstream
    }
    const double genA = (ref.pt[b] - ref.pt[p]) / (ref.pt[b] + ref.pt[p]);
    // category from the reco jets, bins from the gen jets, as upstream
    const bool genSameBin =
        dj::AbsEtaBin(ref.eta[b]) >= 0 &&
        dj::AbsEtaBin(ref.eta[b]) == dj::AbsEtaBin(ref.eta[p]);
    if (!isSM || genSameBin) {
      fillCat(h.gen, ref.eta[b], ref.eta[p], isSM, genPtave, genAlpha, genA);
    }
  }
}

void runDijetProfiles(TString input, TString output, TString modeFlag,
                      Long64_t maxEvents, bool jerClosure) {

  const AnalysisConfig &cfg = Config();
  PrintConfigSummary(cfg);

  // forest, event selection, JEC (+ JER smearing), pT ordering
  ForestEventLoop loop(input, ParseRunMode(modeFlag), maxEvents, jerClosure);
  const RunMode mode = loop.Mode();
  const size_t nCones = loop.NCones();

  BinningConfig bins;
  std::vector<DijetProfiles> cones(nCones);
  std::vector<DiJetJERCHistograms> jerc(nCones);
  for (size_t c = 0; c < nCones; c++) {
    cones[c].Init();
    jerc[c].Init(cfg.coneLabels[c], bins, loop.IsMC());
  }

  while (loop.Next()) {
    const float weight = loop.Weight();

    for (size_t c = 0; c < nCones; c++) {
      const auto &reco = loop.JetsIn(c).reco;
      const float *pt = loop.CorrPt(c);
      const SortedJets &s = loop.Sorted(c);

      // DiJetJERC's own jet cleaning + selection, independent of Dijet2's
      FillDiJetJERC(loop, c, jerc[c], cfg.hltJ80Thresh);
      if (s.sublead == -1) {
        continue;
      }

      // every jet above 15 GeV within |eta| < 3 passes ID + veto map
      // ("allJetsGood", HF exempt, see kMaxIdEta); all jets above 15 GeV
      // but the leading two, HF included, build the jets-only FSR vector
      bool allGood = true;
      double mnX = 0, mnY = 0;
      for (int j = 0; j < reco.nref; j++) {
        if (pt[j] <= dijet2::kJetPtMin) {
          continue;
        }
        if (std::fabs(reco.eta[j]) < dijet2::kMaxIdEta && !loop.GoodJet(c, j)) {
          allGood = false;
          break;
        }
        if (j != s.lead && j != s.sublead) {
          mnX -= pt[j] * std::cos(reco.phi[j]);
          mnY -= pt[j] * std::sin(reco.phi[j]);
        }
      }
      if (!allGood) {
        continue;
      }

      // both leading jets act as tag and probe in turn
      for (int k = 0; k < 2; k++) {
        const int t = (k == 0) ? s.lead : s.sublead;
        const int p = (k == 0) ? s.sublead : s.lead;
        if (std::fabs(reco.eta[t]) >= dijet2::kTagEtaMax ||
            std::fabs(reco.eta[p]) >= dijet2::kMaxIdEta) {
          continue;
        }
        if (pt[t] <= dijet2::kJetPtMin || pt[p] <= dijet2::kJetPtMin) {
          continue;
        }
        if (DijetDPhi(reco.phi[t], reco.phi[p]) <= dijet2::kDphiMin) {
          continue;
        }
        if (std::fabs((pt[t] - pt[p]) / (pt[t] + pt[p])) >= dijet2::kAsymMax) {
          continue;
        }

        const dijet2::Balance b = dijet2::ComputeBalance(
            pt[t], reco.phi[t], pt[p], reco.phi[p], reco.eta[p], mnX, mnY);
        cones[c].Fill(b, weight, OwnsBin(mode, b.ptavp, cfg.hltJ80Thresh),
                      OwnsBin(mode, b.ptTag, cfg.hltJ80Thresh),
                      OwnsBin(mode, b.ptProbe, cfg.hltJ80Thresh));
      }
    }
  }

  // output: <cone>/Dijet2/ as DijetHistosFill writes it, <cone>/DiJetJERC/
  TFile *fo = CreateOutputFile(output);
  loop.WriteEventHists(fo);
  for (size_t c = 0; c < nCones; c++) {
    cones[c].Write(fo->mkdir(cfg.coneLabels[c] + "/Dijet2"));
    jerc[c].Write(fo->mkdir(cfg.coneLabels[c] + "/DiJetJERC"));
  }
  fo->Close();
}
