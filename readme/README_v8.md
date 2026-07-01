# L2Residuals-2024ppRef

Object-oriented C++17/ROOT framework for CMS L2 Residual jet energy corrections for the 2024 pp reference run. The code measures eta-dependent data/MC scale factors with dijet pT balance: a barrel tag jet constrains a probe jet anywhere in eta.

The asymmetry is

```text
A = (pT_probe - pT_tag) / (pT_probe + pT_tag)
```

and is measured in cumulative alpha selections, extrapolated to alpha -> 0, converted to `R_data/R_MC`, and written as CMS JEC text files.

## Current Status

This repository is a work in progress.

- Step 1, `runAsymmetry`, is implemented for hard-probe data, zero-bias data, and MC.
- Step 2, `runResiduals`, is implemented for abs(eta) and full-eta extraction.
- Step 2 QA plotting is implemented in `plotResiduals`.
- Step 1 kinematics plotting is now available through `plotResiduals ... "kinematics"`.
- Step 3, `runTextFile`, is implemented for JEC text output.
- Tests cover dijet construction, eta folding, text-file writing, and build/macro loadability.
- Full-statistics validation is still pending.

Known local state as of this v8 draft:

- `cfg/2024ppRef.h` points to AK2, AK3, AK4, AK5, and AK6 L2Relative JEC files.
- `kHLTJ80Branch` is currently `HLT_AK4PFJet80_v8`; confirm this against the production HiForest.
- `condor/runtime_wrapper.sh` still requires `CMSSW_SRC` before Condor submission.
- A one-sided central-barrel A distribution was observed in a small hard-probe test sample whose event numbers were all even. This may be a sample-parity artifact, but it needs a larger-statistics check.
- `JSON_handler` intentionally exposes only `isGood(run, lumi)`. The previous run-only helper was removed so JSON selection always uses the run+lumi decision.
- ROOT object and plot names now use a shared naming policy in `include/Naming.h`.

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

The main implementation files in `src/` are built into `lib/libl2residuals.dylib` on macOS or `lib/libl2residuals.so` on Linux. The command-line executables and ROOT macro entry points both call into that same compiled library, so there is one physics implementation and two thin front doors.

## Naming Convention

Naming helpers live in `include/Naming.h` under the `L2Name` namespace. They centralize the key order for ROOT objects and plot paths so the convention is not reimplemented by hand in each macro or source file.

The canonical ordering is:

```text
cone / object-or-plot-kind / eta-mode-or-eta-bin / ptavg / alpha / method-or-detail
```

This mirrors the main 4D asymmetry histogram axes after the cone label:

```text
eta_probe, pT_avg, alpha, A
```

`A` itself is usually represented by the object kind, for example `A_data` or `A_mc`. For residual summaries, the method or final detail comes after the axis keys.

Examples of canonical ROOT object names:

```text
ak4PF_A_data_abseta_eta_0p0_0p261_ptavg_40_90_alpha_0p05
ak4PF_R_data_abseta_ptavg_40_90_alpha_0p05_gauss
ak4PF_R_abseta_eta_0p0_0p261_ptavg_40_90_gauss
ak4PF_intercept_abseta_ptavg_40_90_gauss
ak4PF_intercept_abseta_ptavg_40_90_gauss_norm
```

Examples of canonical plot paths:

```text
plots/ak4PF/adist/eta_0p0_0p261/ptavg_40_90/alpha_0p05/adist_ak4PF_eta_0p0_0p261_ptavg_40_90_alpha_0p05.png
plots/ak4PF/alpha/eta_0p0_0p261/ptavg_40_90/alphafit_ak4PF_gauss_eta_0p0_0p261_ptavg_40_90.png
plots/ak4PF/roverlay/ptavg_40_90/alpha_0p05/roverlay_ak4PF_gauss_ptavg_40_90_alpha_0p05.png
plots/ak4PF/kinematics/incl/ptmin_40/kinematics_ak4PF_incl_eta_phi_ptmin_40.png
```

The main helper functions are:

- `L2Name::EtaKey(ieta, fullEta)`: produces edge-stamped keys such as `eta_0p0_0p261`.
- `L2Name::PtKey(ptSlice)`: uses the configured pT slice short name, for example `ptavg_40_90`.
- `L2Name::AlphaKey(alphaSlice)`: uses the configured alpha short name, for example `alpha_0p05`.
- `L2Name::EtaModeKey(fullEta)`: returns `abseta` or `fulleta`.
- `L2Name::ObjectName(...)`: assembles canonical ROOT object names.
- `L2Name::Join(...)`: joins keys with `_` for object names or `/` for paths.

`plotResiduals` and `runTextFile` use the canonical names only; files outside this nomenclature should be regenerated.

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

Executables are written to `bin/`; the shared library is written to `lib/`. Run executables from the repository root so relative paths under `data/`, `cfg/`, and `lib/` resolve correctly.

## Step 1: Fill Asymmetry Histograms

```bash
./bin/runAsymmetry input_data.root output_data.root --hard-probes
./bin/runAsymmetry input_zero_bias.root output_zero_bias.root --zero-bias
./bin/runAsymmetry input_mc.root output_mc.root --mc
```

`runAsymmetry` applies JEC, sorts jets by corrected pT, applies event and jet quality cuts, builds one selected dijet per cone, and fills one THnSparse per cone:

```text
(eta_probe, pT_avg, alpha, A)
```

It also writes Step 1 kinematic control histograms for each cone:

```text
ak*PF_incl
ak*PF_tag
ak*PF_probe
```

These are `TH3D(eta, phi, pT)` histograms for inclusive jets, selected tag jets, and selected probe jets.

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

For each cone, Step 2 runs two passes:

- abs(eta): folds signed eta with `FoldEtaAxis()` and extracts the primary corrections.
- full eta: keeps signed eta for forward/backward symmetry checks.

For each pT_avg slice, alpha threshold, and eta bin, Step 2 projects the A distribution and estimates the mean with:

- Gaussian fit
- 90% truncated mean
- 95% truncated mean

Each method produces `R_data/R_MC` vs alpha threshold. A linear fit over thresholds up to 0.30 gives the alpha -> 0 intercept.

## Plotting

`plotResiduals` supports both Step 2 residual plotting and explicit Step 1 kinematics plotting.

```bash
./bin/plotResiduals residuals.root plots_residuals "all"
./bin/plotResiduals residuals.root plots_residuals "adist alpha roverlay"
./bin/plotResiduals step1_asymmetry.root plots_step1 "kinematics"
```

Available flags:

```text
all etasym methods finals adist roverlay alpha kinematics
```

`all` runs the Step 2 residual plot groups. `kinematics` is explicit because it expects a Step 1 `runAsymmetry` output file, not a Step 2 residual file.

Plot outputs use a cone-first directory structure:

```text
plots/
  ak4PF/
    adist/
      eta_0p0_0p261/
        ptavg_40_90/
          alpha_0p05/
            adist_ak4PF_eta_0p0_0p261_ptavg_40_90_alpha_0p05.png
    alpha/
      eta_0p0_0p261/
        ptavg_40_90/
          alphafit_ak4PF_gauss_eta_0p0_0p261_ptavg_40_90.png
    roverlay/
      ptavg_40_90/
        alpha_0p05/
          roverlay_ak4PF_gauss_ptavg_40_90_alpha_0p05.png
    kinematics/
      incl/
        kinematics_ak4PF_incl_pt.png
        kinematics_ak4PF_incl_eta.png
        kinematics_ak4PF_incl_phi.png
        ptmin_40/
          kinematics_ak4PF_incl_eta_phi_ptmin_40.png
      tag/
      probe/
```

The plot output convention follows the shared naming policy, omitting dimensions that do not apply. `plotResiduals` can still read older ROOT object names for compatibility with already-produced files.

## Step 3: Write JEC Text File

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

For now, failed fits or missing eta/pT information fall back to unity parameters. This is useful during development, but production use should include stricter validation or an explicit QA summary.

## Tests

```bash
cd build
ctest --output-on-failure
```

Individual binaries can also be run from the repository root:

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

## Open Follow-Ups

- Confirm whether `HLT_AK4PFJet80_v8` is the correct production branch.
- Decide whether `runTextFile` should fail fast when inputs are missing or fits fail.
- Decide whether barrel-barrel tag/probe assignment should remain event-parity based or move to a deterministic hash of run/lumi/event.
- Add an end-to-end validation section once larger hadded samples are available.
- Consider promoting this draft to a tracked top-level `README.md` or stop ignoring `readme/`.
