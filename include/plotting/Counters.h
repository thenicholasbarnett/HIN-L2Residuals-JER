#ifndef L2RESIDUALS_PLOTTING_COUNTERS_H
#define L2RESIDUALS_PLOTTING_COUNTERS_H

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TString.h"

#include "plotting/Style.h"
#include "plotting/Utilities.h"
#include "plotting/AsymmetryDistributions.h"
#include "plotting/Kinematics.h" // kNKinematicsCollections, kKinematicsCollections, kNKinematicsPtMins
#include "plotting/ResponsePlots.h" // kNResponseCollections, kResponseCollections, kNResponseVariants, kResponseVariants
#include "AnalysisConfig.h"
#include "Binning.h"
#include "Naming.h"

// progress accounting

inline int CountAsymDistPlots(TFile *fIn, const TString &cone,
                              const BinningConfig &bins, int minEntries) {
  TDirectory *dData = (TDirectory *)fIn->Get(cone + "/QA_data");
  TDirectory *dMC = (TDirectory *)fIn->Get(cone + "/QA_mc");
  if (!dData || !dMC)
    return 0;

  int n = 0;
  const int nEta = (int)kAbsEtaEdges.size() - 1;
  for (const auto &ptSl : bins.ptavgSlices) {
    const TString ptKey = L2Name::PtKey(ptSl);
    for (const auto &aSl : bins.alphaSlices) {
      const TString alphaKey = L2Name::AlphaKey(aSl);
      for (int ie = 0; ie < nEta; ie++) {
        const TString etaKey = L2Name::EtaKey(ie, false);
        const TString dname = L2Name::ObjectName(
            cone, "A_data",
            {L2Name::EtaModeKey(false), etaKey, ptKey, alphaKey});
        const TString mname = L2Name::ObjectName(
            cone, "A_mc", {L2Name::EtaModeKey(false), etaKey, ptKey, alphaKey});
        TH1D *hd = (TH1D *)dData->Get(dname);
        TH1D *hm = (TH1D *)dMC->Get(mname);
        if (hd && hm && hd->GetEntries() >= minEntries)
          n += 3;
      }
    }
  }
  return n;
}

inline int CountROverlayPlots(TFile *fIn, const TString &cone,
                              const BinningConfig &bins, bool useJer = false) {
  TDirectory *dRvals = (TDirectory *)fIn->Get(cone + "/Rvals");
  if (!dRvals)
    return 0;

  int n = 0;
  for (int m = 0; m < kNMethods; m++) {
    for (const auto &aSl : bins.alphaSlices) {
      const TString alphaKey = L2Name::AlphaKey(aSl);
      for (const auto &ptSl : bins.ptavgSlices) {
        const TString ptKey = L2Name::PtKey(ptSl);
        const TString rdName = L2Name::ObjectName(
            cone, CalibKind("R_data", useJer),
            {L2Name::EtaModeKey(false), ptKey, alphaKey}, {kMethodKeys[m]});
        const TString rmName = L2Name::ObjectName(
            cone, CalibKind("R_mc", useJer),
            {L2Name::EtaModeKey(false), ptKey, alphaKey}, {kMethodKeys[m]});
        if (HasH(dRvals, rdName) && HasH(dRvals, rmName))
          n++;
      }
    }
  }
  return n;
}

inline int CountAlphaFitPlots(TFile *fIn, const TString &cone,
                              const BinningConfig &bins, bool useJer = false) {
  TDirectory *dGraphs = (TDirectory *)fIn->Get(cone + "/graphs");
  if (!dGraphs)
    return 0;

  int n = 0;
  const int nEta = (int)kAbsEtaEdges.size() - 1;
  for (int m = 0; m < kNMethods; m++) {
    for (const auto &ptSl : bins.ptavgSlices) {
      const TString ptKey = L2Name::PtKey(ptSl);
      for (int ie = 0; ie < nEta; ie++) {
        const TString etaKey = L2Name::EtaKey(ie, false);
        const TString gname = L2Name::ObjectName(
            cone, CalibKind("R", useJer),
            {L2Name::EtaModeKey(false), etaKey, ptKey}, {kMethodKeys[m]});
        if (HasGraphWithN(dGraphs, {gname + "_norm", gname}, 2))
          n++;
      }
    }
  }

  // combined all-methods overlay: one more plot per (pT slice, eta bin)
  // wherever at least one method has a valid graph
  for (const auto &ptSl : bins.ptavgSlices) {
    const TString ptKey = L2Name::PtKey(ptSl);
    for (int ie = 0; ie < nEta; ie++) {
      const TString etaKey = L2Name::EtaKey(ie, false);
      bool any = false;
      for (int m = 0; m < kNMethods && !any; m++) {
        const TString gname = L2Name::ObjectName(
            cone, CalibKind("R", useJer),
            {L2Name::EtaModeKey(false), etaKey, ptKey}, {kMethodKeys[m]});
        if (HasGraphWithN(dGraphs, {gname + "_norm", gname}, 2))
          any = true;
      }
      if (any)
        n++;
    }
  }
  return n;
}

inline int CountPtFitPlots(TFile *fIn, const TString &cone) {
  TDirectory *dGraphs = (TDirectory *)fIn->Get(cone + "/graphs");
  if (!dGraphs)
    return 0;

  int n = 0;
  for (int m = 0; m < kNMethods; m++) {
    for (int mode = 0; mode < 2; mode++) {
      const bool fullEta = (mode == 1);
      const int nEta =
          fullEta ? (int)kEtaEdges.size() - 1 : (int)kAbsEtaEdges.size() - 1;
      const TString etaMode = L2Name::EtaModeKey(fullEta);
      for (int ie = 0; ie < nEta; ie++) {
        const TString etaKey = L2Name::EtaKey(ie, fullEta);
        const TString gname = L2Name::ObjectName(
            cone, "ptcorr", {etaMode, etaKey}, {kMethodKeys[m]});
        if (HasGraphWithN(dGraphs, {gname + "_norm", gname}, 2))
          n++;
      }
    }
  }
  return n;
}

inline int CountEtaSymPlots(TFile *fIn, const TString &cone,
                            const BinningConfig &bins, bool useJer = false) {
  int n = 0;
  for (int m = 0; m < kNMethods; m++) {
    for (const auto &sl : bins.ptavgSlices) {
      const TString ptKey = L2Name::PtKey(sl);
      const TString nameAbs = L2Name::ObjectName(
          cone, CalibKind("intercept", useJer),
          {L2Name::EtaModeKey(false), ptKey}, {kMethodKeys[m]});
      const TString nameFull = L2Name::ObjectName(
          cone, CalibKind("intercept", useJer),
          {L2Name::EtaModeKey(true), ptKey}, {kMethodKeys[m]});
      if (HasHAny(fIn, {cone + "/" + nameAbs}) &&
          HasHAny(fIn, {cone + "/" + nameFull}))
        n++;
    }
  }
  return n;
}

inline int CountMethodCompPlots(TFile *fIn, const TString &cone,
                                const BinningConfig &bins, bool fullEta,
                                bool useJer = false) {
  int n = 0;
  const TString etaMode = L2Name::EtaModeKey(fullEta);
  for (const auto &sl : bins.ptavgSlices) {
    const TString name =
        L2Name::ObjectName(cone, CalibKind("intercept", useJer),
                           {etaMode, L2Name::PtKey(sl)}, {kMethodKeys[0]});
    if (HasHAny(fIn, {cone + "/" + name}))
      n++;
  }
  return n;
}

inline int CountFinalsPlots(TFile *fIn, const TString &cone,
                            const BinningConfig &bins, bool useJer = false) {
  int n = 0;
  for (int m = 0; m < kNMethods; m++) {
    for (int mode = 0; mode < 2; mode++) {
      const bool fullEta = (mode == 1);
      const TString etaMode = L2Name::EtaModeKey(fullEta);

      bool anyValid = false;
      for (const auto &ptSl : bins.ptavgSlices) {
        const TString name = L2Name::ObjectName(
            cone, CalibKind("intercept", useJer),
            {etaMode, L2Name::PtKey(ptSl)}, {kMethodKeys[m]});
        if (HasHAny(fIn, {cone + "/" + name})) {
          anyValid = true;
          break;
        }
      }
      if (!anyValid) {
        const TString gridName = L2Name::ObjectName(
            cone, CalibKind("corrfinal", useJer), {etaMode}, {kMethodKeys[m]});
        anyValid =
            HasHAny(fIn, {cone + "/" + gridName + "_norm", gridName + "_norm",
                          cone + "/" + gridName, gridName});
      }
      if (anyValid)
        n++;
    }
  }
  return n;
}

// Companion to CountFinalsPlots: PlotFinalsGrid draws the same Step 3
// corrfinal/corrfinal_jer grid as a 2D color map, Step 3 output only.
inline int CountFinalsGridPlots(TFile *fIn, const TString &cone,
                                bool useJer = false) {
  int n = 0;
  for (int m = 0; m < kNMethods; m++) {
    for (int mode = 0; mode < 2; mode++) {
      const bool fullEta = (mode == 1);
      const TString etaMode = L2Name::EtaModeKey(fullEta);
      const TString gridName = L2Name::ObjectName(
          cone, CalibKind("corrfinal", useJer), {etaMode}, {kMethodKeys[m]});
      if (HasHAny(fIn, {cone + "/" + gridName + "_norm", gridName + "_norm",
                        cone + "/" + gridName, gridName}))
        n++;
    }
  }
  return n;
}

// Companion to CountFinalsPlots: PlotFinalsByCone overlays cone sizes
// instead of pT slices, called once (not once per cone).
inline int CountFinalsByConePlots(TFile *fIn, const BinningConfig &bins,
                                  bool useJer = false) {
  const std::vector<TString> &cones = Config().coneLabels;
  int n = 0;
  for (int m = 0; m < kNMethods; m++) {
    for (int mode = 0; mode < 2; mode++) {
      const bool fullEta = (mode == 1);
      const TString etaMode = L2Name::EtaModeKey(fullEta);
      for (const auto &ptSl : bins.ptavgSlices) {
        bool anyValid = false;
        for (const TString &cone : cones) {
          const TString name = L2Name::ObjectName(
              cone, CalibKind("intercept", useJer),
              {etaMode, L2Name::PtKey(ptSl)}, {kMethodKeys[m]});
          if (HasHAny(fIn, {cone + "/" + name})) {
            anyValid = true;
            break;
          }
        }
        if (!anyValid) {
          for (const TString &cone : cones) {
            const TString gridName =
                L2Name::ObjectName(cone, CalibKind("corrfinal", useJer),
                                   {etaMode}, {kMethodKeys[m]});
            if (HasHAny(fIn, {cone + "/" + gridName + "_norm",
                              gridName + "_norm", cone + "/" + gridName,
                              gridName})) {
              anyValid = true;
              break;
            }
          }
        }
        if (anyValid)
          n++;
      }
    }
  }
  return n;
}

// Companion to PlotClosureCompare: counts (method, etaMode, pT slice)
// combinations where both fIn (original) and fInClosure have a value,
// either as a direct Step 2 histogram or a Step 3 grid to project.
inline int CountClosureComparePlots(TFile *fIn, TFile *fInClosure,
                                    const TString &cone,
                                    const BinningConfig &bins,
                                    bool useJer = false) {
  if (!fInClosure) {
    return 0;
  }
  int n = 0;
  for (int m = 0; m < kNMethods; m++) {
    for (int mode = 0; mode < 2; mode++) {
      const bool fullEta = (mode == 1);
      const TString etaMode = L2Name::EtaModeKey(fullEta);
      const TString gridName = L2Name::ObjectName(
          cone, CalibKind("corrfinal", useJer), {etaMode}, {kMethodKeys[m]});
      const std::vector<TString> gridNames = {
          cone + "/" + gridName + "_norm", gridName + "_norm",
          cone + "/" + gridName, gridName};
      for (const auto &ptSl : bins.ptavgSlices) {
        const TString name = L2Name::ObjectName(
            cone, CalibKind("intercept", useJer),
            {etaMode, L2Name::PtKey(ptSl)}, {kMethodKeys[m]});
        const bool hasOrig =
            HasHAny(fIn, {cone + "/" + name}) || HasHAny(fIn, gridNames);
        const bool hasClos = HasHAny(fInClosure, {cone + "/" + name}) ||
                             HasHAny(fInClosure, gridNames);
        if (hasOrig && hasClos)
          n++;
      }
    }
  }
  return n;
}

inline int CountNormCompPlots(TFile *fIn, const TString &cone,
                              const BinningConfig &bins, bool fullEta,
                              bool useJer = false) {
  int n = 0;
  const TString etaMode = L2Name::EtaModeKey(fullEta);
  for (int m = 0; m < kNMethods; m++) {
    bool anyValid = false;
    for (const auto &sl : bins.ptavgSlices) {
      const TString ptKey = L2Name::PtKey(sl);
      const TString nameDirect =
          L2Name::ObjectName(cone, CalibKind("intercept", useJer),
                             {etaMode, ptKey}, {kMethodKeys[m]});
      if (HasHAny(fIn, {cone + "/" + nameDirect}) &&
          HasHAny(fIn, {cone + "/" + nameDirect + "_norm"})) {
        anyValid = true;
        break;
      }
    }
    if (anyValid)
      n++;
  }
  return n;
}

inline int CountKinematicsPlots(TFile *fIn, const TString &cone,
                                bool includeIncl = true) {
  if (!fIn->Get(cone + "/" + cone + "_incl"))
    return 0;

  const int plotsPerCollection = 3 + kNKinematicsPtMins;
  int n = 0;
  bool anyCollection = false;
  for (int ic = 0; ic < kNKinematicsCollections; ic++) {
    const TString coll = kKinematicsCollections[ic].inputKey;
    if (!includeIncl && coll == "incl")
      continue;
    if (fIn->Get(cone + "/" + cone + "_" + coll)) {
      n += plotsPerCollection;
      anyCollection = true;
    }
  }
  if (anyCollection)
    n += 3;
  return n;
}

inline int CountEventPlots(TFile *fIn) {
  int n = 0;
  if (fIn->Get("hvz_all") && fIn->Get("hvz"))
    n++;
  if (fIn->Get("hfilt"))
    n++;
  if (fIn->Get("h_hlt_j80"))
    n++;
  return n;
}

inline int CountResponsePlots(TFile *fIn, const TString &cone,
                              int minEntries) {
  bool any = false;
  for (int ic = 0; ic < kNResponseCollections && !any; ic++) {
    TString name = L2Name::ObjectName(cone, "JES", {"corr", "vs_ptgen"},
                                      {kResponseCollections[ic]});
    if (HasHAny(fIn, {cone + "/" + name}))
      any = true;
  }
  if (!any)
    return 0;

  int n = 0;
  for (int ic = 0; ic < kNResponseCollections; ic++) {
    const TString collection = kResponseCollections[ic];

    for (int iv = 0; iv < kNResponseVariants; iv++) {
      const TString variant = kResponseVariants[iv];
      TString jesName =
          L2Name::ObjectName(cone, "JES", {variant, "vs_ptgen"}, {collection});
      TH1D *hJes = (TH1D *)fIn->Get(cone + "/" + jesName);
      if (hJes) {
        TDirectory *dPt = (TDirectory *)fIn->Get(cone + "/QA_response_ptgen");
        for (int ip = 1; ip <= hJes->GetNbinsX(); ip++) {
          const double lo = hJes->GetXaxis()->GetBinLowEdge(ip);
          const double hi = hJes->GetXaxis()->GetBinUpEdge(ip);
          TString objName =
              L2Name::ObjectName(cone, "response_" + variant,
                                 {L2Name::PtGenKey(lo, hi)}, {collection});
          TH1D *h = dPt ? (TH1D *)dPt->Get(objName) : nullptr;
          if (h && h->GetEntries() >= minEntries)
            n++;
        }
      }
    }
  }

  // per-variant summary: up to 6 canvases per cone (JES/JER x corr/reco/raw),
  // each drawn whenever at least one collection has the corresponding object.
  auto countSummary = [&](const TString &quantity,
                          const std::vector<TString> &keys) -> int {
    for (int ic = 0; ic < kNResponseCollections; ic++) {
      TString name =
          L2Name::ObjectName(cone, quantity, keys, {kResponseCollections[ic]});
      if (HasHAny(fIn, {cone + "/" + name}))
        return 1;
    }
    return 0;
  };
  for (int iv = 0; iv < kNResponseVariants; iv++) {
    n += countSummary("JES", {kResponseVariants[iv], "vs_ptgen"});
    n += countSummary("JER", {kResponseVariants[iv], "vs_ptgen"});
  }

  // variant comparison: up to 6 canvases per cone (JES/JER x incl/tag/probe),
  // each drawn whenever at least one variant has the corresponding object.
  auto countVariantComparison = [&](const TString &collection,
                                    const TString &quantity) -> int {
    for (int iv = 0; iv < kNResponseVariants; iv++) {
      TString name = L2Name::ObjectName(
          cone, quantity, {kResponseVariants[iv], "vs_ptgen"}, {collection});
      if (HasHAny(fIn, {cone + "/" + name}))
        return 1;
    }
    return 0;
  };
  for (int ic = 0; ic < kNResponseCollections; ic++) {
    n += countVariantComparison(kResponseCollections[ic], "JES");
    n += countVariantComparison(kResponseCollections[ic], "JER");
  }

  // detector-region overlay: up to 2 (JES/JER) x 3 (incl/tag/probe) x 2
  // (disjoint/cumulative) canvases per cone, each drawn whenever any coarse
  // |eta| slice has the corresponding object.
  auto countEtaRange = [&](const TString &quantity, const TString &collection,
                           bool cumulative) -> int {
    const TString dirName =
        cumulative ? "JER_per_etarange_cumulative" : "JER_per_etarange";
    const std::vector<RangeBin> slices =
        cumulative ? BuildEtaRangeSlicesCumulative() : BuildEtaRangeSlices();
    for (const RangeBin &sl : slices) {
      TString name = L2Name::ObjectName(
          cone, quantity, {"corr", "vs_ptgen", L2Name::EtaKey(sl.lo, sl.hi)},
          {collection});
      if (HasHAny(fIn, {cone + "/" + dirName + "/" + name})) {
        return 1;
      }
    }
    return 0;
  };
  for (int ic = 0; ic < kNResponseCollections; ic++) {
    for (int icum = 0; icum < 2; icum++) {
      const bool cumulative = (icum == 1);
      n += countEtaRange("JES", kResponseCollections[ic], cumulative);
      n += countEtaRange("JER", kResponseCollections[ic], cumulative);
    }
  }

  return n;
}

#endif
