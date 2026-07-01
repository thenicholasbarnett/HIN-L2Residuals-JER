# L2Residuals-2024ppRef

Object-oriented C++17/ROOT framework for CMS L2 Residual jet energy corrections for the 2024 pp reference run. The code measures eta-dependent data/MC scale factors with dijet pT balance: a barrel tag jet constrains a probe jet anywhere in eta.

The asymmetry is

```text
A = (pT_probe - pT_tag) / (pT_probe + pT_tag)
```

and is measured in cumulative alpha selections, extrapolated to alpha -> 0, converted to `R_data/R_MC`, and written as CMS JEC text files.

## Current Status

This repository is a work in progress.

- Step 1, `runAsymmetry`, is implemented for hard-probe data, zero-bias data, and MC. Short flags `-hp`, `-zb`, `-mc` are available alongside the long forms. Output is organized into per-cone TDirectories.
- Step 2, `runResiduals`, is implemented for abs(eta) and full-eta extraction. Output is organized into per-cone TDirectories with named sub-directories.
- Step 3, `runTextFile`, is implemented for JEC text output. Reads from cone TDirectories using the canonical naming policy.
- `plotResiduals` covers both Step 1 (`event`, `kinematics`) and Step 2 (`etasym`, `methods`, `finals`, `adist`, `roverlay`, `alpha`) diagnostics. The `all` flag runs everything. All canvases are 800×600 with `SetRealAspectRatio` active in interactive sessions.
- Tests cover dijet construction, eta folding, text-file writing, and build/macro loadability.
- Full-statistics validation is still pending.

Known local state:

- `cfg/2024ppRef.h` points to AK2, AK3, AK4, AK5, and AK6 L2Relative JEC files.
- `kHLTJ80Branch` is currently `HLT_AK4PFJet80_v8`; confirm this against the production HiForest.
- `condor/runtime_wrapper.sh` still requires `CMSSW_SRC` before Condor submission.
- A one-sided central-barrel A distribution was observed in a small hard-probe test sample whose event numbers were all even. This may be a sample-parity artifact and needs a larger-statistics check.
- `JSON_handler` intentionally exposes only `isGood(run, lumi)`.
- ROOT object and plot names use a shared naming policy in `include/Naming.h`.

## Repository Layout

```text
include/      public declarations, config-facing types, and header-only helpers
src/          compiled implementation for the main pipeline
macros/       ROOT/Cling entry points and compiled executable wrappers
cfg/          2024 pp reference configuration
data/jec/     tracked JEC text inputs
data/json/    tracked golden JSON input
data/veto/    tracked jet veto map
condor/       HTCondor submission helpers
tests/        focused unit-style and build tests
```

The main implementation files in `src/` are built into `lib/libl2residuals.dylib` on macOS or `lib/libl2residuals.so` on Linux. The command-line executables and ROOT macro entry points both call into that compiled library.

## Naming Convention

Naming helpers live in `include/Naming.h` under the `L2Name` namespace.

The canonical ordering is:

```text
cone / object-or-plot-kind / eta-mode-or-eta-bin / ptavg / alpha / method-or-detail
```

Examples of canonical ROOT object names:

```text
ak4PF_A_data_abseta_eta_0p0_0p261_ptavg_40_90_alpha_0p05
ak4PF_R_data_abseta_ptavg_40_90_alpha_0p05_gauss
ak4PF_intercept_abseta_ptavg_40_90_gauss
ak4PF_intercept_abseta_ptavg_40_90_gauss_norm
```

Examples of canonical plot paths:

```text
plots/ak4PF/adist/eta_0p0_0p261/ptavg_40_90/alpha_0p05/adist_ak4PF_eta_0p0_0p261_ptavg_40_90_alpha_0p05.png
plots/ak4PF/alpha/eta_0p0_0p261/ptavg_40_90/alphafit_ak4PF_gauss_eta_0p0_0p261_ptavg_40_90.png
plots/ak4PF/kinematics/incl/ptmin_40/kinematics_ak4PF_incl_eta_phi_ptmin_40.png
plots/event/event_vz.png
```

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
cmake --build .
```

Executables are written to `bin/`; the shared library is written to `lib/`. Run executables from the repository root so relative paths under `data/`, `cfg/`, and `lib/` resolve correctly. Output directories are created automatically if they do not exist.

## Step 1: Fill Asymmetry Histograms

```bash
./bin/runAsymmetry input_data.root output_data.root --hard-probes   # or -hp
./bin/runAsymmetry input_zb.root   output_zb.root   --zero-bias     # or -zb
./bin/runAsymmetry input_mc.root   output_mc.root   --monte-carlo   # or -mc
```

`runAsymmetry` fails loudly on invalid mode strings. For each event and cone, it applies JEC, sorts jets by corrected pT, applies quality cuts, builds one selected dijet, and fills one THnSparse per cone.

### Step 1 output file layout

Event-level QA histograms are written at the root of the file. Per-cone histograms are organized into TDirectories named after each cone:

```text
hvz_all           — vz before any cuts
hvz               — vz after vz + vertex filter
hfilt             — pprimaryVertexFilter (DATA only)
h_hlt_j80         — HLT_AK4PFJet80 bits (hard-probe DATA only)
ak2PF/
  ak2PF_asym      — THnSparse (eta_probe, pT_avg, alpha, A)
  ak2PF_incl      — TH3D (eta, phi, pT) for all corrected jets passing kMinPt
  ak2PF_tag       — TH3D (eta, phi, pT) for selected tag jets
  ak2PF_probe     — TH3D (eta, phi, pT) for selected probe jets
ak3PF/ ... ak6PF/ — same structure
```

Primary selections:

- `abs(vz) < 15 cm`
- DATA: primary vertex filter
- DATA: golden JSON via `isGood(run, lumi)`
- hard-probe DATA: configured `HLT_AK4PFJet80` branch plus leading AK4 corrected pT > 100 GeV
- MC: event weight from the `weight` branch
- Dijet: `DeltaPhi > 2.7`, subleading corrected pT > 10 GeV, at least one barrel jet, both leading jets pass tight jet ID and veto map
- If a third jet exists but fails jet ID, it is treated as absent and alpha is set to 0

## Step 2: Extract Residuals

```bash
./bin/runResiduals data_hadded.root mc_hadded.root residuals.root
```

For each cone, Step 2 runs two passes: abs(eta) and full eta. For each pT_avg slice, alpha threshold, and eta bin, Step 2 projects the A distribution and estimates the mean with a Gaussian fit, 90% truncated mean, and 95% truncated mean. Each method produces `R_data/R_MC` vs alpha threshold; a linear fit over thresholds up to 0.30 gives the alpha -> 0 intercept. A kFSR-normalized variant is also stored.

`runResiduals` reads the THnSparse from the cone TDirectory in the Step 1 file.

### Step 2 output file layout

All objects are organized into per-cone TDirectories. Within each cone directory, sub-directories hold intermediate objects and the intercept TH1Ds sit directly in the cone directory:

```text
ak4PF/
  ak4PF_asym_data_abseta    — folded THnSparse for abs(eta) extraction
  ak4PF_asym_mc_abseta
  QA_data/                  — A distribution TH1Ds for data, abs(eta)
  QA_mc/                    — A distribution TH1Ds for MC, abs(eta)
  QA_data_fulleta/          — A distribution TH1Ds for data, full eta
  QA_mc_fulleta/
  graphs/                   — TGraphErrors of R vs alpha with embedded fit (both eta modes)
  Rvals/                    — R_data and R_mc TH1Ds vs |eta|, abs(eta)
  Rvals_fulleta/            — R_data and R_mc TH1Ds vs eta, full eta
  ak4PF_intercept_abseta_ptavg_40_90_gauss         — intercept TH1D
  ak4PF_intercept_abseta_ptavg_40_90_gauss_norm    — kFSR-normalized variant
  ... (all methods, all pT slices, abs(eta) and fulleta)
ak2PF/ ak3PF/ ak5PF/ ak6PF/ — same structure
```

## Plotting

`plotResiduals` handles both Step 1 and Step 2 diagnostics through a single binary. The `all` flag runs everything; `event` and `kinematics` skip gracefully when passed a Step 2 file.

```bash
# Run everything on a Step 2 residuals file
./bin/plotResiduals residuals.root plots_out "all"

# Step 2 residual plots only
./bin/plotResiduals residuals.root plots_out "etasym methods finals adist roverlay alpha"

# Step 1 event-level QA (vz, filter, trigger)
./bin/plotResiduals step1_asym.root plots_step1 "event"

# Step 1 jet kinematics
./bin/plotResiduals step1_asym.root plots_step1 "kinematics"
```

Available flags:

```text
etasym     — full-eta vs |eta| reflected symmetry check
methods    — method comparison: gauss vs trunc90 vs trunc95
finals     — R_data/R_MC at alpha->0, all pT slices overlaid
adist      — A distributions per (eta, pT, alpha) bin, log-y, truncation lines
roverlay   — R_data and R_MC vs |eta| overlay with ratio panel
alpha      — R vs alpha threshold with linear fit, all 9 alpha points shown
kinematics — Step-1 inclusive/tag/probe pT, eta, phi projections and eta-phi maps
event      — Step-1 event-level QA: vz before/after cuts, ppvF, HLT trigger
all        — run all of the above (default)
```

Directory lookups use the canonical `cone/subdir` path.

When a flag's required histograms are absent — for example Step 2 flags on a Step 1 file, or `kinematics`/`event` on a Step 2 file — the function ticks the progress bar and returns with no output. The bar loads at the normal rate regardless of which file type is passed.

All canvases are 800×600 pixels. ROOT's `SetRealAspectRatio(4/3)` is called on each canvas and is active in interactive sessions; it is suppressed in batch mode where the canvas size is already fixed at creation.

Plot outputs follow the naming convention, organized cone-first:

```text
plots/
  event/
    event_vz.png
    event_ppvF.png        (DATA only)
    event_hlt_j80.png     (hard-probe DATA only)
  ak4PF/
    kinematics/incl/...
    adist/eta_0p0_0p261/ptavg_40_90/alpha_0p05/...
    alpha/eta_0p0_0p261/ptavg_40_90/...
    roverlay/ptavg_40_90/alpha_0p05/...
    etasym/ptavg_40_90/...
    methods/abseta/ptavg_40_90/...
    finals/abseta/...
```

## Step 3: Write JEC Text File

```bash
./bin/runTextFile residuals.root L2Residuals_2024ppRef_AK4.txt
./bin/runTextFile residuals.root L2Residuals_2024ppRef_AK4.txt gauss ak4PF
./bin/runTextFile residuals.root L2Residuals_2024ppRef_AK4.txt trunc90 ak4PF
```

Methods: `gauss`, `trunc90`, `trunc95`
Cones: `ak2PF`, `ak3PF`, `ak4PF`, `ak5PF`, `ak6PF`

Reads intercept TH1Ds from the cone TDirectory in the residuals file. For each of 18 |eta| bins, fits the pT dependence across 5 slice centers with `1/(p0 + p1*log10(0.01*pT) + p2/(pT/10))`. Writes a 36-line CMS JEC text file (negative eta outer->inner, then positive inner->outer). Failed fits or fewer than 3 valid pT slices fall back to unity.

## Condor (Step 1)

Step 1 runs one job per HiForest file on HTCondor at lxplus. Build the project on lxplus first (the Mac `.dylib` cannot be transferred — a Linux `.so` is required).

```bash
# Set CMSSW_SRC in condor/runtime_wrapper.sh, then:
bash condor/make_condor.sh JOBNAME filelist.txt /eos/cms/store/user/nbarnett/output --hard-probes
bash condor/make_condor.sh JOBNAME filelist.txt /eos/path --no-submit   # dry run
```

`filelist.txt` — one absolute EOS/AFS path per line. After all jobs finish, hadd before Step 2:

```bash
hadd data_hadded.root output_JOBNAME_*.root
hadd mc_hadded.root   output_mc_*.root
```

## Tests

```bash
cd build && ctest --output-on-failure
```

Individual binaries from repo root:

```bash
./bin/TestDijet
./bin/TestFoldEtaAxis
./bin/TestTextFileWriter
```

| Test | Count | Coverage |
|------|-------|----------|
| `TestDijet` | 37 | `FindLeadingJets`, `MakeDijet` — barrel assignment, Δφ cut, α, pT_avg, A sign, dphi wrap |
| `TestFoldEtaAxis` | 7 | Axis properties, symmetric fill, Sumw2 propagation, arbitrary axis index |
| `TestTextFileWriter` | 7 | File structure, eta ordering, mirror symmetry, unity fallback, fit round-trip, CMS eta extent |
| `TestBuild` | 9 | Executables exist and are runnable, usage message, library loadable, macro interpretable |

## Configuration

`cfg/2024ppRef.h` is the only file that needs editing to adapt to a new run period:

```text
kJECFilesPerCone   — L2Relative text file paths, one list per cone (AK2–6)
kVetoMapPath       — jet veto map ROOT file
kJSONPath          — golden JSON (currently runs 387474–387721)
kHiTreePath        — hiEvtAnalyzer/HiTree
kSkimTreePath      — skimanalysis/HltTree
kTrigTreePath      — hltanalysis/HltTree
kJetTreePaths      — ak{2,3,4,5,6}PFJetAnalyzer/t
kHLTJ80Branch      — HLT_AK4PFJet80_v8  (confirm version suffix from actual HiForest)
kHLTJ80Thresh      — 100.0 GeV
kConeLabels        — {ak2PF, ak3PF, ak4PF, ak5PF, ak6PF}
kTrigCone          — ak4PF
```

## Open Follow-Ups

- Confirm `HLT_AK4PFJet80_v8` against the production HiForest branch.
- Set `CMSSW_SRC` in `condor/runtime_wrapper.sh` before first Condor submission.
- Validate the full pipeline end-to-end on a large hadded sample.
- Decide whether `runTextFile` should fail fast when inputs are missing or fits fail.
- Confirm the one-sided central-barrel A distribution resolves with full-statistics data.
