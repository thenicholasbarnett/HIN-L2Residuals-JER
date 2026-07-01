# L2Residuals-2024ppRef

Object-oriented C++17/ROOT framework for CMS L2 Residual jet energy corrections for the 2024 pp reference run. The core design is a rewrite of `DijetResiduals`: one DATA/MC executable, one 4D THnSparse per cone, and a separate extraction step that folds signed eta into abs(eta) only when needed.

L2 Residuals are eta-dependent data/MC scale factors applied after MC-truth JEC. This repo measures them with dijet pT balance: a barrel tag jet constrains a probe jet anywhere in eta. The asymmetry

```text
A = (pT_probe - pT_tag) / (pT_probe + pT_tag)
```

is measured in cumulative alpha slices, extrapolated to alpha -> 0, converted to R, and written as CMS JEC text files.

---

## Current Status

This repository is a work in progress.

- Step 1 (`runAsymmetry`) is implemented for hard-probe data, zero-bias data, and MC.
- Step 2 (`runResiduals`) is implemented for abs(eta) and full-eta extraction.
- Step 2 QA plotting (`plotResiduals`) is implemented.
- Step 3 (`runTextFile`) is implemented for JEC text output.
- Unit-style tests cover dijet construction, eta folding, text-file writing, and build/macro loadability.
- Full-statistics validation is still pending.

Known local state as of this v7 draft:

- `cfg/2024ppRef.h` points to all five cone JEC files: AK2, AK3, AK4, AK5, AK6.
- `kHLTJ80Branch` is currently `HLT_AK4PFJet80_v8`; confirm against the actual HiForest used for production.
- `condor/runtime_wrapper.sh` still requires `CMSSW_SRC` to be set before submission.
- A one-sided central-barrel `A` distribution has been observed in a small hard-probe test sample whose event numbers were all even. That may be a sample-parity artifact, but it needs a larger sample check.

---

## Build

Prerequisites:

- ROOT 6
- CMake >= 3.10
- C++17 compiler
- `nlohmann-json`

```bash
git clone https://github.com/thenicholasbarnett/L2Residuals-2024ppref
cd L2Residuals-2024ppref
mkdir build
cd build
cmake ..
make
```

Executables are written to `bin/`; the shared library is written to `lib/`. Run executables from the repository root so relative paths under `data/`, `cfg/`, and `lib/` resolve correctly.

---

## Pipeline

### Step 1: Fill Asymmetry Histograms

```bash
./bin/runAsymmetry input_data.root output_data.root
./bin/runAsymmetry input_data.root output_data.root --hard-probes
./bin/runAsymmetry input_zero_bias.root output_zero_bias.root --zero-bias
./bin/runAsymmetry input_mc.root output_mc.root --mc
```

`runAsymmetry` now fails loudly on invalid mode strings. Valid modes are:

- `--hard-probes`
- `--zero-bias`
- `--mc`

For each event and cone, Step 1 applies JEC, sorts jets by corrected pT, applies event and jet quality cuts, builds one selected dijet, and fills:

```text
(eta_probe, pT_avg, alpha, A)
```

The stored eta axis is signed eta with 36 CMS JEC bins. The alpha axis is a uniform 50-bin axis from 0 to 0.5; cumulative alpha thresholds are recovered during extraction with `SetRangeUser(0, threshold)`.

Primary selections:

- `abs(vz) < 15 cm`
- DATA: primary vertex filter
- DATA: golden JSON
- hard-probe DATA: `HLT_AK4PFJet80` branch from config plus leading AK4 corrected pT > 100 GeV
- MC: event weight from the `weight` branch
- Dijet: `DeltaPhi > 2.7`, subleading corrected pT > 10 GeV, at least one barrel jet, both leading jets pass tight jet ID and veto map

If a third jet exists but fails jet ID, it is treated as absent and alpha is set to 0.

### Step 2: Extract Residuals

```bash
./bin/runResiduals data_hadded.root mc_hadded.root residuals.root
```

For each cone, Step 2 runs two passes:

- abs(eta): folds the signed-eta sparse with `FoldEtaAxis()` and extracts the primary corrections.
- full eta: keeps signed eta for forward/backward symmetry checks.

For each pT_avg slice, alpha threshold, and eta bin, Step 2 projects the A distribution and estimates the mean with:

- Gaussian fit
- 90% truncated mean
- 95% truncated mean

Each method produces `R_data/R_MC` vs alpha threshold. A linear fit over thresholds up to 0.30 gives the alpha -> 0 intercept.

### Step 2 QA: Plot Residuals

```bash
./bin/plotResiduals residuals.root
./bin/plotResiduals residuals.root my_plots
```

Useful plot groups include:

- eta symmetry: full eta vs reflected abs(eta)
- method comparison: Gaussian vs truncated means
- A distributions by cone, pT_avg, alpha threshold, and eta bin

### Step 3: Write JEC Text File

```bash
./bin/runTextFile residuals.root L2Residuals_2024ppRef_AK4.txt
./bin/runTextFile residuals.root L2Residuals_2024ppRef_AK4.txt gauss ak4PF
./bin/runTextFile residuals.root L2Residuals_2024ppRef_AK4.txt trunc90 ak4PF
```

Methods:

- `gauss`
- `trunc90`
- `trunc95`

Cones:

- `ak2PF`
- `ak3PF`
- `ak4PF`
- `ak5PF`
- `ak6PF`

For now, failed fits or missing eta/pT information fall back to unity parameters:

```text
---
```

This behavior is intentionally marked for later discussion. A production version should probably fail fast or write an explicit QA summary when too many bins fall back to unity.

---

## Configuration

Main configuration lives in `cfg/2024ppRef.h`.

Important fields:

- `kJECFilesPerCone`: L2Relative text files for AK2-6 PF
- `kVetoMapPath`: jet veto map ROOT file
- `kJSONPath`: golden JSON
- `kHiTreePath`: event tree
- `kSkimTreePath`: skim/filter tree
- `kTrigTreePath`: trigger tree
- `kJetTreePaths`: one jet tree per cone
- `kHLTJ80Branch`: trigger branch name; currently `HLT_AK4PFJet80_v8`
- `kHLTJ80Thresh`: trigger efficiency pT threshold
- `kConeLabels`: cone labels used in output object names
- `kTrigCone`: cone used for event-level trigger efficiency kinematics

---

## Tests

```bash
cd build
ctest --output-on-failure
```

Individual binaries can also be run from the repo root:

```bash
./bin/TestDijet
./bin/TestFoldEtaAxis
./bin/TestTextFileWriter
```

Current focused coverage:

- `TestDijet`: leading-jet ordering, barrel assignment, delta-phi cut, alpha, pT_avg, A sign, phi wrapping
- `TestFoldEtaAxis`: signed eta to abs(eta) folding and uncertainty propagation
- `TestTextFileWriter`: text structure, eta ordering, unity fallback, fit round trip
- `TestBuild`: executable existence, usage behavior, shared-library loading, ROOT macro interpretability

---

## Notes For The Final README

Items to settle before promoting this draft to the top-level README:

- Confirm whether the hard-probe trigger branch should be `HLT_AK4PFJet80_v8` for the target HiForest.
- Decide how strict `runTextFile` should be when inputs are missing or fits fail.
- Decide whether barrel-barrel tag/probe assignment should remain event-parity based or move to a deterministic hash of run/lumi/event.
- Add a small end-to-end validation section once a larger hadded sample is available.
- Move stable project guidance into top-level `README.md` and keep scratch/debug notes in a local developer-notes file.
