#include "TextFileWriter.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"

#include "Binning.h"
#include "Naming.h"
#include "AnalysisConfig.h"
#include "jetmet_jer/JetResolutionObject.h"

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cmath>

// min pT slices for 3-parameter fit
static constexpr int kMinSlices = 3;

// outputTag default and fixed output dir (see TextFileWriter.h)
static const char *const kDefaultTag = "L2Residual";
static const char *const kDefaultJerTag = "JER_SF";
static const char *const kTextOutputSubdir = "data/jec/preliminary";

// NEED TO CHANGE
// REFERENCE BINNING IN TOML
static constexpr double kPtLo = 40.0;
static constexpr double kPtHi = 1000.0;

static constexpr int kNPar = 3;

// Fit function: 1/(p0 + p1*log10(0.01*x) + p2/(x/10))
static Double_t FitFunc(Double_t *x, Double_t *par) {
  double denom =
      par[0] + par[1] * TMath::Log10(0.01 * x[0]) + par[2] / (x[0] / 10.0);
  return (denom != 0.0) ? 1.0 / denom : 1.0;
}

// JEC formula declaration
static const char *const kJECHeader =
    "{1 JetEta 1 JetPt 1./([p0]+[p1]*log10(0.01*x)+[p2]/(x/10.0)) Correction "
    "L2Residual}";

struct FitResult {
  double p[kNPar] = {1.0, 0.0, 0.0};
  bool valid = false;
};

// midpoint of slice, x-value in fit
static double SliceCenter(const RangeBin &sl) { return 0.5 * (sl.lo + sl.hi); }

// get intercept histogram
static TH1D *FetchIntercept(TFile *f, const TString &cone,
                            const TString &name) {
  TDirectory *coneDir = (TDirectory *)f->Get(cone);
  return coneDir ? (TH1D *)coneDir->Get(name) : nullptr;
}

// fit corr(pT_avg) per eta bin, write TGraphErrors to dGraphs, return the fit parameters
static FitResult FitPtSlices(const std::vector<double> &ptCenters,
                             const std::vector<double> &corr,
                             const std::vector<double> &corrErr,
                             const TString &graphName, TDirectory *dGraphs) {
  FitResult r;
  int n = (int)ptCenters.size();
  if (n < kMinSlices) {
    return r;
  }

  std::vector<double> ex(n, 0.0);
  TGraphErrors *gr = new TGraphErrors(n, ptCenters.data(), corr.data(),
                                      ex.data(), corrErr.data());
  gr->SetName(graphName);
  gr->SetTitle(";p_{T,avg} [GeV];Correction factor");

  TF1 *f = new TF1(graphName + "_fit", FitFunc, kPtLo, kPtHi, kNPar);
  // start at unity correction (1/(1+0+0)=1)
  // Starting at (1.5,1.5,1.5) gives a negative denominator at low pT
  f->SetParameter(0, 1.0);
  f->SetParameter(1, 0.0);
  f->SetParameter(2, 0.0);

  // fit function in graph so runPlotting can draw
  TFitResultPtr res = gr->Fit(f, "QSR");
  if (res.Get() && res->IsValid()) {
    for (int p = 0; p < kNPar; p++) {
      r.p[p] = res->Parameter(p);
    }
    r.valid = true;
  }

  if (dGraphs) {
    dGraphs->cd();
    gr->Write();
  }
  delete f; // clone in gr's function list is separately owned
  delete gr;
  return r;
}

// one CMS L2Residual JEC data line for the [etaLo, etaHi] range
static void WriteJECLine(std::ofstream &out, double etaLo, double etaHi,
                         const FitResult &fit) {
  out << etaLo << "\t" << etaHi << "\t" << (kNPar + 2) << "\t" << kPtLo << "\t"
      << kPtHi;
  if (fit.valid) {
    out << "\t" << fit.p[0] << "\t" << fit.p[1] << "\t" << fit.p[2];
  } else {
    out << "\t" << 1 << "\t" << 0 << "\t" << 0;
  }
  out << "\n";
}

// mirrored |eta| about eta = 0 for |eta| results text file
static bool WriteAbsEtaTextFile(const TString &path,
                                const std::vector<FitResult> &fits) {
  std::ofstream out(path.Data());
  if (!out.is_open()) {
    return false;
  }
  out << kJECHeader << "\n";
  const int nEta = (int)fits.size();
  for (int ieta = nEta - 1; ieta >= 0; ieta--) {
    WriteJECLine(out, -kAbsEtaEdges[ieta + 1], -kAbsEtaEdges[ieta], fits[ieta]);
  }
  for (int ieta = 0; ieta < nEta; ieta++) {
    WriteJECLine(out, kAbsEtaEdges[ieta], kAbsEtaEdges[ieta + 1], fits[ieta]);
  }
  out.close();
  return true;
}

// full-eta text file, each bin has its own fit
static bool WriteFullEtaTextFile(const TString &path,
                                 const std::vector<FitResult> &fits) {
  std::ofstream out(path.Data());
  if (!out.is_open()) {
    return false;
  }
  out << kJECHeader << "\n";
  const int nEta = (int)fits.size();
  for (int ieta = 0; ieta < nEta; ieta++) {
    WriteJECLine(out, kEtaEdges[ieta], kEtaEdges[ieta + 1], fits[ieta]);
  }
  out.close();
  return true;
}

// JER SF text output, via the vendored JME::JetResolutionObject
// JER SF file is a direct binned grid: one value per (eta bin, pT_avg slice) bin
// Definition line: "{1 JetEta 1 JetPt [0] Resolution}"
// 1 bin variable (JetEta), 1 structural variable (JetPt), [0] is JER SF value
// unc of SF  = error/value, retrievable via record.getParametersValues()[1]
// JER-smearing convention (s_up/down = sJER * (1 +/- unc))
struct JerRecord {
  double ptLo, ptHi, value, unc;
};

// JER SF values (eta_probe, pT_avg), skipping slices with no data
static std::vector<JerRecord>
CollectJerRecords(int ieta, const std::vector<RangeBin> &ptSlices,
                  const std::vector<TH1D *> &hSliceJer) {
  std::vector<JerRecord> out;
  for (size_t ip = 0; ip < hSliceJer.size(); ip++) {
    if (!hSliceJer[ip]) {
      continue;
    }
    double v = hSliceJer[ip]->GetBinContent(ieta + 1);
    double e = hSliceJer[ip]->GetBinError(ieta + 1);
    if (v == 0.0 && e == 0.0) {
      continue;
    }
    double unc = (v != 0.0) ? e / v : 0.0;
    out.push_back({ptSlices[ip].lo, ptSlices[ip].hi, v, unc});
  }
  return out;
}

// One record line per JerRecord
// "4" token is 2*nVariables + nParameters
// 2 for JetPt range, 2 for [value, unc] parameters
// see JetResolutionObject::Record
static void AppendJerLines(std::stringstream &buf, double etaLo, double etaHi,
                           const std::vector<JerRecord> &recs) {
  for (const auto &r : recs) {
    buf << etaLo << " " << etaHi << " 4 " << r.ptLo << " " << r.ptHi << " "
        << r.value << " " << r.unc << "\n";
  }
}

// JER SF temp text file first
// JME::JetResolutionObject: parses and re-emitts with saveToFile()
// mirrors abseta cells onto both eta halves
static bool WriteJerSfTextFile(const TString &path, bool fullEta,
                               const std::vector<RangeBin> &ptSlices,
                               const std::vector<TH1D *> &hSliceJer) {

  std::stringstream buf;
  buf << "{1 JetEta 1 JetPt [0] Resolution}\n";

  if (fullEta) {
    const int nEta = (int)kEtaEdges.size() - 1;
    for (int ieta = 0; ieta < nEta; ieta++) {
      AppendJerLines(buf, kEtaEdges[ieta], kEtaEdges[ieta + 1],
                     CollectJerRecords(ieta, ptSlices, hSliceJer));
    }
  } else {
    const int nEta = (int)kAbsEtaEdges.size() - 1;
    for (int ieta = nEta - 1; ieta >= 0; ieta--) {
      AppendJerLines(buf, -kAbsEtaEdges[ieta + 1], -kAbsEtaEdges[ieta],
                     CollectJerRecords(ieta, ptSlices, hSliceJer));
    }
    for (int ieta = 0; ieta < nEta; ieta++) {
      AppendJerLines(buf, kAbsEtaEdges[ieta], kAbsEtaEdges[ieta + 1],
                     CollectJerRecords(ieta, ptSlices, hSliceJer));
    }
  }

  TString tmpPath = path + ".tmp";
  {
    std::ofstream tmp(tmpPath.Data());
    if (!tmp.is_open()) {
      return false;
    }
    tmp << buf.str();
  }

  bool ok = true;
  try {
    JME::JetResolutionObject obj(tmpPath.Data());
    obj.saveToFile(path.Data());
  } catch (const std::exception &e) {
    std::cerr << "ERROR building JER SF text file " << path << ": " << e.what()
              << "\n";
    ok = false;
  }
  gSystem->Unlink(tmpPath);
  return ok;
}

// which residuals given pT_avg slice
// NonTriggeredOnly ignores threshold
// TriggeredOnly drops slices below threshold
enum class SourceMode { Merge, TriggeredOnly, NonTriggeredOnly };

static TFile *SelectSource(SourceMode mode, const RangeBin &ptSlice,
                           double threshold, TFile *fTrig, TFile *fNoTrig) {
  switch (mode) {
  case SourceMode::Merge:
    return (ptSlice.lo >= threshold) ? fTrig : fNoTrig;
  case SourceMode::TriggeredOnly:
    return (ptSlice.lo >= threshold) ? fTrig : nullptr;
  case SourceMode::NonTriggeredOnly:
    return fNoTrig;
  }
  return nullptr;
}

// fTrig/fNoTrig already-open files
// calibMode selects which calibration kind is read/written
// mutually exclusive within one call.
static void RunTextFileImpl(SourceMode srcMode, TFile *fTrig, TFile *fNoTrig,
                            CalibrationMode calibMode, TString outputRootFile,
                            TString outputTag, TString method, bool useNorm) {

  if (outputTag.IsNull()) {
    // JER SF isn't an L2Residual -- give it its own default, mirroring
    // runTextFilePtResolution's separately-scoped "JER_ptresolution" below
    outputTag =
        (calibMode == CalibrationMode::JER) ? kDefaultJerTag : kDefaultTag;
  }
  if (outputTag.Contains("/")) {
    std::cerr
        << "ERROR: TAG must be a plain name, not a path (it contains '/'): "
        << outputTag << "\n";
    return;
  }

  const TString suffix = useNorm ? "_norm" : "";

  const AnalysisConfig &cfg = Config();
  PrintConfigSummary(cfg);

  // calibration text files land here
  const TString textDir =
      TString(cfg.repoRoot.c_str()) + "/" + kTextOutputSubdir;
  if (gSystem->mkdir(textDir, kTRUE) < 0 && gSystem->AccessPathName(textDir)) {
    std::cerr << "ERROR: cannot create text output directory: " << textDir
              << "\n";
    return;
  }

  static const char *const kModeLabel[] = {"merge (triggered + non-triggered)",
                                           "single (triggered-only)",
                                           "single (non-triggered-only)"};
  const bool doJEC = (calibMode == CalibrationMode::JEC);
  std::cout << "Mode: " << kModeLabel[(int)srcMode] << "\n"
            << "Calibration:  " << (doJEC ? "JEC" : "JER") << "\n"
            << "Triggered residuals:     "
            << (fTrig ? fTrig->GetName() : "(none)") << "\n"
            << "Non-triggered residuals: "
            << (fNoTrig ? fNoTrig->GetName() : "(none)") << "\n"
            << "Output ROOT:  " << outputRootFile << "\n"
            << "Text output dir: " << textDir << "\n"
            << "Tag:          " << outputTag << "\n"
            << "Method:       " << method << "\n"
            << "Normalized:   " << (useNorm ? "yes (kFSR-norm)" : "no (direct)")
            << "\n";

  {
    Ssiz_t sl = outputRootFile.Last('/');
    if (sl != kNPOS) {
      gSystem->mkdir(TString(outputRootFile(0, sl)), kTRUE);
    }
  }

  TFile *fOut = new TFile(outputRootFile, "recreate");
  const bool wantAbsEta = (cfg.etaModeOutput != "eta");
  const bool wantFullEta = (cfg.etaModeOutput != "abseta");

  BinningConfig bins;
  const int nPt = (int)bins.ptavgSlices.size();

  // pT_avg bin edges for corrfinal grid
  std::vector<Double_t> ptEdges(nPt + 1);
  ptEdges[0] = bins.ptavgSlices[0].lo;
  for (int ip = 0; ip < nPt; ip++) {
    ptEdges[ip + 1] = bins.ptavgSlices[ip].hi;
  }

  const TString objKind = doJEC ? "intercept" : "intercept_jer";

  for (const TString &cone : cfg.coneLabels) {

    TDirectory *trigConeDir = fTrig ? (TDirectory *)fTrig->Get(cone) : nullptr;
    TDirectory *notrigConeDir =
        fNoTrig ? (TDirectory *)fNoTrig->Get(cone) : nullptr;
    if ((fTrig && !trigConeDir) || (fNoTrig && !notrigConeDir)) {
      std::cerr << "Missing " << cone
                << " directory in a residuals file, skipping\n";
      continue;
    }

    fOut->cd();
    TDirectory *coneDirOut = fOut->mkdir(cone.Data());
    TDirectory *dGraphs = doJEC ? coneDirOut->mkdir("graphs") : nullptr;

    // filled below for abseta then fulleta
    std::vector<FitResult> fitsAbsEta, fitsFullEta;
    std::vector<TH1D *> hSliceAbsEta, hSliceFullEta;

    for (int em = 0; em < 2; em++) {
      const bool fullEta = (em == 1);
      if (fullEta && !wantFullEta) {
        continue;
      }
      if (!fullEta && !wantAbsEta) {
        continue;
      }
      const std::vector<Double_t> &etaEdges =
          fullEta ? kEtaEdges : kAbsEtaEdges;
      const int nEta = (int)etaEdges.size() - 1;
      const TString etaMode = L2Name::EtaModeKey(fullEta);

      // one correction histogram per pT slice
      std::vector<TH1D *> &hSlice = fullEta ? hSliceFullEta : hSliceAbsEta;
      hSlice.assign(nPt, nullptr);
      int nMissingSlices = 0;
      for (int ip = 0; ip < nPt; ip++) {
        const auto &ptSlice = bins.ptavgSlices[ip];
        TFile *src =
            SelectSource(srcMode, ptSlice, cfg.hltJ80Thresh, fTrig, fNoTrig);
        if (src) {
          TString name =
              L2Name::ObjectName(cone, objKind,
                                 {etaMode, L2Name::PtKey(ptSlice)}, {method}) +
              suffix;
          hSlice[ip] = FetchIntercept(src, cone, name);
        }
        if (!hSlice[ip]) {
          nMissingSlices++;
        }
      }
      if (nMissingSlices > 0) {
        std::cerr << cone << " " << etaMode << ": " << nMissingSlices << "/"
                  << nPt
                  << " pT slices missing or excluded (below trigger "
                     "threshold in triggered-only mode, or "
                  << objKind << " histogram not found)\n";
      }

      // final corrections grid: x = eta/|eta|, y = pT_avg slices, z =
      // correction/JER SF -- direct binned values for both JEC and JER SF
      // (JER SF has no pT-dependence fit, see below)
      const TString gridKind = doJEC ? "corrfinal" : "corrfinal_jer";
      TString gridName =
          L2Name::ObjectName(cone, gridKind, {etaMode}, {method}) + suffix;
      TH2D *hGrid =
          new TH2D(gridName, "", nEta, etaEdges.data(), nPt, ptEdges.data());
      hGrid->GetXaxis()->SetTitle(fullEta ? "#eta_{probe}" : "|#eta_{probe}|");
      hGrid->GetYaxis()->SetTitle("p_{T,avg} [GeV]");
      hGrid->GetZaxis()->SetTitle(doJEC ? "Correction factor" : "JER SF");
      hGrid->Sumw2();
      for (int ip = 0; ip < nPt; ip++) {
        if (!hSlice[ip]) {
          continue;
        }
        for (int ieta = 0; ieta < nEta; ieta++) {
          hGrid->SetBinContent(ieta + 1, ip + 1,
                               hSlice[ip]->GetBinContent(ieta + 1));
          hGrid->SetBinError(ieta + 1, ip + 1,
                             hSlice[ip]->GetBinError(ieta + 1));
        }
      }
      coneDirOut->cd();
      hGrid->Write();
      delete hGrid;

      if (!doJEC) {
        continue; // JER SF has no pT-dependence fit -- flat grid only
      }

      // pT-dependence fit, one per eta bin
      std::vector<FitResult> &fits = fullEta ? fitsFullEta : fitsAbsEta;
      fits.assign(nEta, FitResult{});
      int nUnityFallback = 0;
      for (int ieta = 0; ieta < nEta; ieta++) {
        std::vector<double> ptX, corr, corrErr;
        for (int ip = 0; ip < nPt; ip++) {
          if (!hSlice[ip]) {
            continue;
          }
          double v = hSlice[ip]->GetBinContent(ieta + 1);
          double e = hSlice[ip]->GetBinError(ieta + 1);
          if (v == 0.0 && e == 0.0) {
            continue;
          }
          ptX.push_back(SliceCenter(bins.ptavgSlices[ip]));
          corr.push_back(v);
          corrErr.push_back(e > 0.0 ? e : 1e-4);
        }
        TString graphName =
            L2Name::ObjectName(cone, "ptcorr",
                               {etaMode, L2Name::EtaKey(ieta, fullEta)},
                               {method}) +
            suffix;
        fits[ieta] = FitPtSlices(ptX, corr, corrErr, graphName, dGraphs);
        if (!fits[ieta].valid) {
          nUnityFallback++;
        }
      }
      if (nUnityFallback > 0) {
        std::cerr << cone << " " << etaMode << ": " << nUnityFallback << "/"
                  << nEta << " eta bins fell back to unity (fewer than "
                  << kMinSlices << " pT slices had data)\n";
      }
    }

    if (doJEC) {
      TString absEtaTxt =
          textDir + "/" + outputTag + "_" + cone + "_abseta" + suffix + ".txt";
      TString etaTxt =
          textDir + "/" + outputTag + "_" + cone + "_eta" + suffix + ".txt";
      if (wantAbsEta && !WriteAbsEtaTextFile(absEtaTxt, fitsAbsEta)) {
        std::cerr << "Cannot open output file " << absEtaTxt << "\n";
      }
      if (wantFullEta && !WriteFullEtaTextFile(etaTxt, fitsFullEta)) {
        std::cerr << "Cannot open output file " << etaTxt << "\n";
      }

      std::cout << "Done. " << cone << ": ";
      if (wantAbsEta) {
        std::cout << 2 * (int)fitsAbsEta.size() << " eta bins -> " << absEtaTxt;
      }
      if (wantAbsEta && wantFullEta) {
        std::cout << ", ";
      }
      if (wantFullEta) {
        std::cout << (int)fitsFullEta.size() << " eta bins -> " << etaTxt;
      }
      std::cout << "\n";
    } else {
      TString absEtaJerTxt = textDir + "/" + outputTag + "_" + cone +
                             "_abseta_jer" + suffix + ".txt";
      TString etaJerTxt =
          textDir + "/" + outputTag + "_" + cone + "_eta_jer" + suffix + ".txt";
      if (wantAbsEta && !WriteJerSfTextFile(absEtaJerTxt, false,
                                            bins.ptavgSlices, hSliceAbsEta)) {
        std::cerr << "Cannot write JER SF output file " << absEtaJerTxt << "\n";
      }
      if (wantFullEta && !WriteJerSfTextFile(etaJerTxt, true, bins.ptavgSlices,
                                             hSliceFullEta)) {
        std::cerr << "Cannot write JER SF output file " << etaJerTxt << "\n";
      }

      std::cout << "Done (JER SF). " << cone << ": ";
      if (wantAbsEta) {
        std::cout << "-> " << absEtaJerTxt;
      }
      if (wantAbsEta && wantFullEta) {
        std::cout << ", ";
      }
      if (wantFullEta) {
        std::cout << "-> " << etaJerTxt;
      }
      std::cout << "\n";
    }
  }

  fOut->Close();
  if (fTrig) {
    fTrig->Close();
  }
  if (fNoTrig) {
    fNoTrig->Close();
  }
}

void runTextFile(TString triggeredResidualsFile,
                 TString nonTriggeredResidualsFile, TString outputRootFile,
                 CalibrationMode mode, TString outputTag, TString method,
                 bool useNorm) {
  TFile *fTrig = TFile::Open(triggeredResidualsFile, "read");
  TFile *fNoTrig = TFile::Open(nonTriggeredResidualsFile, "read");
  if (!fTrig || fTrig->IsZombie()) {
    std::cerr << "Cannot open " << triggeredResidualsFile << "\n";
    return;
  }
  if (!fNoTrig || fNoTrig->IsZombie()) {
    std::cerr << "Cannot open " << nonTriggeredResidualsFile << "\n";
    return;
  }
  RunTextFileImpl(SourceMode::Merge, fTrig, fNoTrig, mode, outputRootFile,
                  outputTag, method, useNorm);
}

void runTextFile(TString residualsFile, SingleDatasetKind kind,
                 TString outputRootFile, CalibrationMode mode,
                 TString outputTag, TString method, bool useNorm) {
  TFile *f = TFile::Open(residualsFile, "read");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << residualsFile << "\n";
    return;
  }
  if (kind == SingleDatasetKind::Triggered) {
    RunTextFileImpl(SourceMode::TriggeredOnly, f, nullptr, mode, outputRootFile,
                    outputTag, method, useNorm);
  } else {
    RunTextFileImpl(SourceMode::NonTriggeredOnly, nullptr, f, mode,
                    outputRootFile, outputTag, method, useNorm);
  }
}

// eta bin (|eta| if fullEta=false, signed eta_reco if fullEta=true) JER
// records read back from runResponse's JER_per_abseta/JER_per_eta directory
static std::vector<std::vector<JerRecord>>
CollectPtResEtaRecords(TDirectory *dPerEta, const TString &cone, int nEta,
                       bool fullEta, const TString &method) {
  std::vector<std::vector<JerRecord>> etaRecords(nEta);
  for (int ieta = 0; ieta < nEta; ieta++) {
    TString etaKey = L2Name::EtaKey(ieta, fullEta);
    TString hName = L2Name::ObjectName(cone, "JER", {"corr", "vs_ptgen", etaKey},
                                       {"incl", method});
    TH1D *h = (TH1D *)dPerEta->Get(hName);
    if (!h) {
      continue;
    }
    for (int ip = 1; ip <= h->GetNbinsX(); ip++) {
      const double v = h->GetBinContent(ip);
      const double e = h->GetBinError(ip);
      if (v == 0.0 && e == 0.0) {
        continue;
      }
      const double ptLo = h->GetXaxis()->GetBinLowEdge(ip);
      const double ptHi = h->GetXaxis()->GetBinUpEdge(ip);
      const double unc = (v != 0.0) ? e / v : 0.0;
      etaRecords[ieta].push_back({ptLo, ptHi, v, unc});
    }
  }
  return etaRecords;
}

// round-trips an assembled JER-grid buffer through JME::JetResolutionObject,
// same convention as WriteJerSfTextFile above
static bool WritePtResolutionFile(const TString &path,
                                  const std::stringstream &buf) {
  TString tmpPath = path + ".tmp";
  {
    std::ofstream tmp(tmpPath.Data());
    tmp << buf.str();
  }
  bool ok = true;
  try {
    JME::JetResolutionObject obj(tmpPath.Data());
    obj.saveToFile(path.Data());
  } catch (const std::exception &ex) {
    std::cerr << "ERROR writing pT resolution file " << path << ": "
              << ex.what() << "\n";
    ok = false;
  }
  gSystem->Unlink(tmpPath);
  return ok;
}

void runTextFilePtResolution(TString responseFile, TString outputTag,
                             TString method) {
  const AnalysisConfig &cfg = Config();

  TFile *fIn = TFile::Open(responseFile, "read");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "Cannot open response file: " << responseFile << "\n";
    return;
  }

  if (outputTag.IsNull()) {
    outputTag = "JER_ptresolution";
  }
  if (method.IsNull()) {
    method = cfg.defaultMethod;
  }
  if (outputTag.Contains("/")) {
    std::cerr << "ERROR: outputTag must not contain '/': " << outputTag << "\n";
    fIn->Close();
    return;
  }

  TString textDir = TString(cfg.repoRoot.c_str()) + "/" + kTextOutputSubdir;
  if (gSystem->mkdir(textDir, kTRUE) < 0 && gSystem->AccessPathName(textDir)) {
    std::cerr << "ERROR: cannot create text output directory: " << textDir
              << "\n";
    fIn->Close();
    return;
  }

  const bool wantAbsEta = (cfg.etaModeOutput != "eta");
  const bool wantFullEta = (cfg.etaModeOutput != "abseta");

  const int nAbsEta = (int)kAbsEtaEdges.size() - 1;
  const int nFullEta = (int)kEtaEdges.size() - 1;

  for (const TString &cone : cfg.coneLabels) {
    if (wantAbsEta) {
      TDirectory *dPerAbsEta = (TDirectory *)fIn->Get(cone + "/JER_per_abseta");
      if (!dPerAbsEta) {
        std::cerr
            << "No JER_per_abseta/ directory for " << cone
            << " -- was this file produced with abs-eta output enabled?\n";
      } else {
        auto etaRecords =
            CollectPtResEtaRecords(dPerAbsEta, cone, nAbsEta, false, method);

        TString path =
            textDir + "/" + outputTag + "_" + cone + "_abseta_ptresolution.txt";
        std::stringstream buf;
        buf << "{1 JetEta 1 JetPt [0] Resolution}\n";
        for (int ieta = nAbsEta - 1; ieta >= 0; ieta--) {
          AppendJerLines(buf, -kAbsEtaEdges[ieta + 1], -kAbsEtaEdges[ieta],
                         etaRecords[ieta]);
        }
        for (int ieta = 0; ieta < nAbsEta; ieta++) {
          AppendJerLines(buf, kAbsEtaEdges[ieta], kAbsEtaEdges[ieta + 1],
                         etaRecords[ieta]);
        }

        if (WritePtResolutionFile(path, buf)) {
          std::cout << "Done (pT resolution). " << cone << ": -> " << path
                    << "\n";
        }
      }
    }

    if (wantFullEta) {
      TDirectory *dPerEta = (TDirectory *)fIn->Get(cone + "/JER_per_eta");
      if (!dPerEta) {
        std::cerr
            << "No JER_per_eta/ directory for " << cone
            << " -- was this file produced with full-eta output enabled?\n";
      } else {
        auto etaRecords =
            CollectPtResEtaRecords(dPerEta, cone, nFullEta, true, method);

        TString path =
            textDir + "/" + outputTag + "_" + cone + "_eta_ptresolution.txt";
        std::stringstream buf;
        buf << "{1 JetEta 1 JetPt [0] Resolution}\n";
        for (int ieta = 0; ieta < nFullEta; ieta++) {
          AppendJerLines(buf, kEtaEdges[ieta], kEtaEdges[ieta + 1],
                         etaRecords[ieta]);
        }

        if (WritePtResolutionFile(path, buf)) {
          std::cout << "Done (pT resolution). " << cone << ": -> " << path
                    << "\n";
        }
      }
    }
  }

  fIn->Close();
}
