#ifndef DIJETJERCHISTOGRAMS_H
#define DIJETJERCHISTOGRAMS_H

#include "TDirectory.h"
#include "THnSparse.h"
#include "TString.h"

#include "Binning.h"
#include "Utilities.h"

#include <cmath>

// Inputs for JME's dijet JER SF code (UHH2 DiJetJERC, JERSF_Analysis/JER/
// wide_eta_binning/mainRun.cxx), filled from HiForest jets. Selection and
// fill coordinates follow DiJetJERC's preselection
// (src/AnalysisModule_DiJetTrg.cxx) and histogram step
// (JERSF_Analysis/hist_preparation/MC/wide_eta_bin/MySelector_full_DiJet.C),
// so an exporter only has to project, rebin and rename. Per category:
//   SM : both leading jets in the same |eta| bin
//   FE : one jet |eta| < kBarrelEta, the other in the bin (forward extension)
// sparse axes (|eta| of the jet in the bin, pT_ave, alpha, value), value =
// A = (pT_barrel - pT_probe) / sum (reco, gen) or pT_reco / pT_gen (mctruth,
// both jets). The SM/SM_control and FE_reference/FE_control/FE splits are
// |eta| ranges, left to the exporter.

namespace dijetjerc {

static constexpr double kBarrelEta = 1.131;  // s_eta_barr
static constexpr double kDphiMin = 2.7;      // s_delta_phi
static constexpr double kJetPtMin = 10.0;    // minJetPt, jet cleaner
static constexpr double kGenJetPtMin = 15.0; // minGenJetPt, gen jet cleaner
static constexpr double kJet3PtMin = 15.0;   // jet_thr, alpha
static constexpr double kGenThr = 10.0;      // gen_thr
static constexpr double kMatchDR = 0.3;      // s_delta_R, mctruth

// |eta| bin index in kAbsEtaEdges (their eta_L2R), -1 outside
inline int AbsEtaBin(double eta) {
  const double a = std::fabs(eta);
  for (size_t i = 0; i + 1 < kAbsEtaEdges.size(); i++) {
    if (a >= kAbsEtaEdges[i] && a < kAbsEtaEdges[i + 1]) {
      return (int)i;
    }
  }
  return -1;
}

} // namespace dijetjerc

struct DiJetJERCHistograms {
  enum Cat { kSM = 0, kFE = 1 };
  static constexpr int kNCat = 2;

  THnSparse *reco[kNCat] = {nullptr, nullptr};
  THnSparse *gen[kNCat] = {nullptr, nullptr};
  THnSparse *mctruth[kNCat] = {nullptr, nullptr};

  void Init(const TString &prefix, const BinningConfig &bins, bool isMC) {
    TDirectory::TContext detached(nullptr);
    const char *catName[kNCat] = {"sm", "fe"};
    const AxisBins asym = {200, -1.0, 1.0, "A"};
    for (int k = 0; k < kNCat; k++) {
      const TString base = prefix + "_jerc_" + catName[k];
      reco[k] = Make(base, bins, asym);
      if (isMC) {
        gen[k] = Make(base + "_gen", bins, asym);
        mctruth[k] = Make(base + "_mctruth", bins, bins.response);
      }
    }
  }

  static void Fill(THnSparse *h, double eta, double pt, double alpha,
                   double value, double w) {
    const double x[4] = {std::fabs(eta), pt, alpha, value};
    h->Fill(x, w);
  }

  void Write(TDirectory *dir) const {
    TDirectory::TContext ctx(dir);
    for (int k = 0; k < kNCat; k++) {
      for (THnSparse *h : {reco[k], gen[k], mctruth[k]}) {
        if (h) {
          h->Write();
        }
      }
    }
  }

private:
  static THnSparse *Make(const TString &name, const BinningConfig &bins,
                         const AxisBins &value) {
    THnSparse *h = MakeTHnSparse<THnSparseD>(
        name, "", {bins.abseta, bins.ptavg, bins.alpha, value});
    SetAbsEtaBins(h, 0);
    h->Sumw2();
    return h;
  }
};

#endif
