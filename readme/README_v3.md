# L2Residuals-2024ppRef

> Object-oriented framework for CMS L2 Residual jet energy corrections — 2024 pp reference run

---

## What This Does

Jet energies measured by CMS need several layers of calibration before they can be used in an analysis. This framework computes the **L2 Residual** correction: the final eta-dependent scale factor that brings data into agreement with simulation after MC-truth corrections have already been applied.

The measurement is based on dijet pT balance. In a back-to-back dijet event, a well-measured barrel jet constrains how well the detector sees a jet anywhere in eta. Running over millions of such events and comparing data to MC gives the correction as a function of jet eta and pT.

The output is a text file in the standard JEC format, ready to be chained into `JetCorrector`.

---

## Quick Start

**Requirements:** ROOT 6, CMake ≥ 3.10, C++17

```bash
git clone https://github.com/thenicholasbarnett/L2Residuals-2024ppRef
cd L2Residuals-2024ppRef
mkdir build && cd build
cmake ..
make
```

**Run the analysis:**

```bash
# Step 1 — fill asymmetry histograms
./runAsymmetry input.root output.root --data
./runAsymmetry input.root output.root --mc

# Step 2 — extract corrections
./runResiduals data_asymmetry.root mc_asymmetry.root residuals.root

# Step 3 — write JEC text file
./runTextFile residuals.root L2Residuals_2024ppRef.txt
```

---

## How It Works

The analysis runs in three steps:

**Step 1 — Asymmetry generation**

For each back-to-back dijet event, one jet is assigned as the tag (in the barrel, |eta| < 1.3) and one as the probe (anywhere). The asymmetry A = (pT_probe − pT_tag) / (pT_probe + pT_tag) is filled into a 4D sparse histogram binned in probe eta, pT_avg, and alpha (third-jet activity). A single executable handles both data and MC, configured by a flag.

**Step 2 — Residuals extraction**

The response R = (1 + <A>) / (1 − <A>) is computed from the mean asymmetry in each bin. The data/MC double ratio is extrapolated linearly to alpha = 0 to remove radiation bias. The result is the correction factor as a function of eta and pT_avg.

**Step 3 — Text file output**

The pT dependence in each eta bin is fit with a 3-parameter function and written to a JEC-standard text file. Negative eta is mirrored from positive. The file is immediately usable with `JetCorrector`.

---

## Repository Layout

```
include/    all headers — structs, utilities, jet corrections, binning
src/        implementations — compiled into a shared library
macros/     entry-point executables — thin wrappers around the library
configs/    run-specific settings (paths, triggers, correction files)
```

Key headers:

| Header | Purpose |
|--------|---------|
| `JetStruct.h` | reco/ref/gen jet arrays with automated branch binding |
| `EventStructs.h` | event-level variables (vz, weight, filters) |
| `BranchMapping.h` | one-line TTree branch setup from a struct |
| `Dijet.h` | tag/probe selection, asymmetry, alpha |
| `BinningConfig.h` | all histogram binning in one place |
| `JetCorrector.h` | Yi Chen v3.0 — chains text-file JEC levels |
| `JetSelection.h` | pp jet ID cuts and veto map |
| `Utilities.h` | THnSparse projections, plot style, canvas helpers |

---

## Adapting for a New Run Period

Add a new config header under `configs/` specifying the correction file paths, trigger pT thresholds, and jet algorithm. Pass it to `AsymmetryAnalysis` at construction. Everything else stays the same.

---

## Contact

Nicholas Barnett — thenicholasbarnett@gmail.com
