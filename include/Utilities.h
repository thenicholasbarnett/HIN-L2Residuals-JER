#ifndef UTILITIES_H
#define UTILITIES_H

#include "TH1F.h"
#include "TH1I.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TArrayD.h"
#include "TGraphErrors.h"

#include <initializer_list>
#include <utility>
#include <vector>

#include "Binning.h"

// Histograms //

template <typename T> void WriteAll(T *h) {
  if (h)
    h->Write();
}

template <typename T> void WriteAll(std::vector<T> &v) {
  for (auto &entry : v)
    WriteAll(entry);
}

// THnSparse

template <typename T>
T *MakeTHnSparse(const TString &name, const TString &title,
                 std::initializer_list<AxisBins> axes) {
  const std::vector<AxisBins> axVec(axes);
  const Int_t ndim = axVec.size();
  std::vector<Int_t> nbins;
  std::vector<Double_t> xmin, xmax;
  for (const auto &ax : axVec) {
    nbins.push_back(ax.nBins);
    xmin.push_back(ax.lo);
    xmax.push_back(ax.hi);
  }
  T *h = new T(name, title, ndim, nbins.data(), xmin.data(), xmax.data());
  for (Int_t i = 0; i < ndim; i++) {
    h->GetAxis(i)->SetTitle(axVec[i].title);
  }
  return h;
}

struct SparseRange {
  Int_t axis;
  Double_t lo;
  Double_t hi;
};

inline void SetAxisRange(THnSparse *h, Int_t axis, Double_t lo, Double_t hi) {
  h->GetAxis(axis)->SetRangeUser(lo, hi);
}

inline void ResetAxisRange(THnSparse *h, Int_t axis) {
  h->GetAxis(axis)->SetRange(0, 0);
}

inline TH1D *ProjectTHnSparse1D(THnSparse *h, Int_t projAxis,
                                const std::vector<SparseRange> &ranges = {},
                                const TString &suffix = "") {
  for (const auto &r : ranges)
    SetAxisRange(h, r.axis, r.lo, r.hi);
  TH1D *out = h->Projection(projAxis);
  out->SetName(out->GetName() + suffix);
  out->SetTitle("");
  for (const auto &r : ranges)
    ResetAxisRange(h, r.axis);
  return out;
}

inline TH2D *ProjectTHnSparse2D(THnSparse *h, Int_t xAxis, Int_t yAxis,
                                const std::vector<SparseRange> &ranges = {},
                                const TString &suffix = "") {
  for (const auto &r : ranges)
    SetAxisRange(h, r.axis, r.lo, r.hi);
  TH2D *out = h->Projection(yAxis, xAxis);
  out->SetName(out->GetName() + suffix);
  out->SetTitle("");
  for (const auto &r : ranges)
    ResetAxisRange(h, r.axis);
  return out;
}

// Fold a THnSparse along a symmetric eta axis (full eta → |eta|).
// The eta axis must be symmetric around 0 with an even number of bins.
// A new THnSparse is returned with the eta axis replaced by a |eta| axis
// of half the bin count, with the same positive-half bin edges.
// All other axes are copied unchanged.

inline THnSparse *FoldEtaAxis(THnSparse *h, Int_t etaAxis,
                              const TString &newName = "") {
  const Int_t ndim = h->GetNdimensions();
  const Int_t nEta = h->GetAxis(etaAxis)->GetNbins();
  const Int_t nAbsEta = nEta / 2;

  // build axis parameters for the new THnSparse
  std::vector<Int_t> nbins(ndim);
  std::vector<Double_t> xmin(ndim), xmax(ndim);
  for (Int_t i = 0; i < ndim; i++) {
    TAxis *ax = h->GetAxis(i);
    if (i == etaAxis) {
      nbins[i] = nAbsEta;
      xmin[i] = 0.0;
      xmax[i] = ax->GetXmax();
    } else {
      nbins[i] = ax->GetNbins();
      xmin[i] = ax->GetXmin();
      xmax[i] = ax->GetXmax();
    }
  }

  TString name = newName.IsNull() ? TString(h->GetName()) + "_abseta" : newName;
  THnSparse *hfold = new THnSparseD(name, h->GetTitle(), ndim, nbins.data(),
                                    xmin.data(), xmax.data());
  hfold->Sumw2();

  // copy axis titles and variable bin edges
  for (Int_t i = 0; i < ndim; i++) {
    TAxis *axOld = h->GetAxis(i);
    TAxis *axNew = hfold->GetAxis(i);
    if (i == etaAxis) {
      axNew->SetTitle("|" + TString(axOld->GetTitle()) + "|");
      // positive-half bin edges: the upper half of the original edge array
      const TArrayD *edges = axOld->GetXbins();
      if (edges && edges->GetSize() > 0) {
        std::vector<Double_t> absEdges(nAbsEta + 1);
        for (Int_t e = 0; e <= nAbsEta; e++)
          absEdges[e] = edges->At(nAbsEta + e);
        axNew->Set(nAbsEta, absEdges.data());
      }
    } else {
      axNew->SetTitle(axOld->GetTitle());
      const TArrayD *edges = axOld->GetXbins();
      if (edges && edges->GetSize() > 0)
        axNew->Set(axOld->GetNbins(), edges->GetArray());
    }
  }

  // iterate over all filled bins and fold
  std::vector<Int_t> coords(ndim);
  std::vector<Int_t> newCoords(ndim);

  for (Long64_t bin = 0; bin < h->GetNbins(); bin++) {
    Double_t val = h->GetBinContent(bin, coords.data());
    Double_t err = h->GetBinError(bin);

    // map eta bin → |eta| bin (both 1-indexed)
    // negative half: bin i ≤ nAbsEta maps to |eta| bin (nAbsEta - i + 1)
    // positive half: bin i > nAbsEta maps to |eta| bin (i - nAbsEta)
    Int_t etaBin = coords[etaAxis];
    Int_t absEtaBin =
        (etaBin <= nAbsEta) ? (nAbsEta - etaBin + 1) : (etaBin - nAbsEta);

    for (Int_t i = 0; i < ndim; i++)
      newCoords[i] = (i == etaAxis) ? absEtaBin : coords[i];

    Double_t oldVal = hfold->GetBinContent(newCoords.data());
    Double_t oldErr = hfold->GetBinError(newCoords.data());
    hfold->SetBinContent(newCoords.data(), oldVal + val);
    hfold->SetBinError(newCoords.data(),
                       TMath::Sqrt(oldErr * oldErr + err * err));
  }

  // hfold is built via SetBinContent/SetBinError, never Fill() -- so its own
  // entries counter would count unique filled-bin combos, not real entries.
  // Folding doesn't change the total, only remaps eta, so copy it from h
  // directly. Without this, Projection() off hfold underreports entries and
  // can silently fail downstream CanFit() guards.
  hfold->SetEntries(h->GetEntries());

  return hfold;
}

inline void DrawRefLine(float xmin, float xmax, float y = 1.0) {
  TLine *line = new TLine(xmin, y, xmax, y);
  line->SetLineWidth(1);
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);
  line->Draw("same");
}

// Fit guards: check whether a histogram or graph has enough data before fitting.
//   CanFit(obj)              : at least 1 entry/point anywhere
//   CanFit(obj, n)           : at least n entries/points anywhere
//   CanFit(obj, {lo, hi})    : at least 1 entry/point in [lo, hi]
//   CanFit(obj, {lo, hi}, n) : at least n entries/points in [lo, hi]

inline bool CanFit(const TH1 *h, int minEntries = 1) {
  return h && h->GetEntries() >= minEntries;
}

inline bool CanFit(TH1 *h, std::pair<double, double> range,
                   int minEntries = 1) {
  if (!h)
    return false;
  int blo = h->FindBin(range.first);
  int bhi = h->FindBin(range.second);
  return h->Integral(blo, bhi) >= minEntries;
}

inline bool CanFit(const TGraphErrors *gr, int minPts = 1) {
  return gr && gr->GetN() >= minPts;
}

inline bool CanFit(const TGraphErrors *gr, std::pair<double, double> range,
                   int minPts = 1) {
  if (!gr)
    return false;
  int n = 0;
  for (int i = 0; i < gr->GetN(); i++)
    if (gr->GetX()[i] >= range.first - 1e-9 &&
        gr->GetX()[i] <= range.second + 1e-9)
      n++;
  return n >= minPts;
}

// TH1

template <typename T> T *MakeTH1(const TString &hname, const AxisBins &bins) {
  T *h = new T(hname, hname, bins.nBins, bins.lo, bins.hi);
  h->SetTitle(bins.title);
  return h;
}

// TH3D(eta, phi, pT) with variable-width CMS JEC eta bins on X.
// Pass kEtaEdges for full eta or kAbsEtaEdges for |eta|.
// TH3D has no mixed uniform/variable constructor, so uniform phi and pT
// edge arrays are built from the AxisBins and the all-variable overload is used.
inline TH3D *MakeTH3DEtaPhiPt(const TString &name,
                              const std::vector<Double_t> &etaEdges,
                              const AxisBins &phi, const AxisBins &pt,
                              const TString &etaTitle = "#eta") {
  std::vector<Double_t> phiEdges(phi.nBins + 1), ptEdges(pt.nBins + 1);
  for (int i = 0; i <= phi.nBins; i++)
    phiEdges[i] = phi.lo + i * (phi.hi - phi.lo) / phi.nBins;
  for (int i = 0; i <= pt.nBins; i++)
    ptEdges[i] = pt.lo + i * (pt.hi - pt.lo) / pt.nBins;

  TH3D *h =
      new TH3D(name.Data(), "", (Int_t)etaEdges.size() - 1, etaEdges.data(),
               phi.nBins, phiEdges.data(), pt.nBins, ptEdges.data());
  h->GetXaxis()->SetTitle(etaTitle);
  h->GetYaxis()->SetTitle(phi.title);
  h->GetZaxis()->SetTitle(pt.title);
  h->Sumw2();
  return h;
}

inline void NormalizeTH1(TH1 *h) {
  const double integral = h->Integral();
  if (integral > 0.0)
    h->Scale(1.0 / integral);
}

#endif
