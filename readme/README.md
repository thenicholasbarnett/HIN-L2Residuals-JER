
<h1> L2Residual Jet Energy Corrections </h1>

<b> Dijet Residuals via p<sub>T</sub>-Balance </b>
Residual jet energy corrections are determined to account for differences in data and simulation by enforcing conservation of transverse momentum. This repository derives L2Residual corrections for the pp reference run collected for the heavy ion collisions in 2024 by CMS. Build this C++ project with cmake to generate standalone binaries and use provided bash sctipts to execute them with HTCondor. The versatility of this work flow is in the ability to rebin almost any dimension, overlay different methods of abstraction or calculation, and run any of the steps as compiled executables or interpret the code directly with cling using ROOT.

<b> Step 1 - Dijet Asymmetries from HiForest files </b>
<br>
<i> Plot various jet collection kinematics and event information. </i>

<b> Step 2 - Residual Corrections Determination </b>
<br>
<i> Plot asymmetries, k<sub>FSR</sub> extrapolations, data to simulation response ratios, method comparisons, and more! </i>

<b> Step 3 - Text file for Analysis Application </b>
<br>
<i> Plot correction factors vs p<sub>T</sub><sup>avg</sup> fitting and extrapolation. </i>
<br><br>

<h1> Quick Start </h1>

<strong> Install & Build </strong>

```
git clone git@github.com:thenicholasbarnett/L2Residuals-2024ppref.git
cd L2Residuals-2024ppref
mkdir build && cd build
cmake .. && make
cd ..
```

> Rebuild: `rm -rf build bin lib && mkdir build && cd build && cmake .. && make && cd ..`

<strong> Batch Process </strong>

```
bash ./condor/make_condor.sh <output_dir> <input_HiForest_filelist.txt>
```

```
bash ./condor/batch_hadd.sh <asymmetry_file.root> <input_glob> <batch_size> <N_parallel>
```

<strong> Derive Corrections </strong>

```
./bin/runResiduals <data_asymmetries_file.root> <mc_asymmetries_file.root> <output_residuals_file.root>
```

```
./bin/runTextFile.C <residuals_hp_file.root> <residuals_zb_file.root> <output_corrections_file.txt>
```

<strong> Plot </strong>

```
./macros/plotResiduals <output_plots_dir> <asymmetry_file.root>
```

```
./macros/plotResiduals <output_plots_dir> <residuals_file.root>
```

<h1> Usage </h1> 

<h3> Dependencies </h3>

- ROOT (6.x, with RIO, Tree, Hist, MathCore)
- CMake ≥ 3.10
- C++17
- [nlohmann/json](https://github.com/nlohmann/json) (for applying golden JSON)
<br>

<h3> Clone </h3>

```
git clone git@github.com:thenicholasbarnett/L2Residuals-2024ppref.git
```

<h3> Build </h3>

```
cmake -B build
cmake --build build
```

Binaries found in `bin/`, shared library in `lib/`, build files in `build/`

<h3> Configuration </h3>

Analysis-specific paths and labels live in:

```
cfg/2024ppRef.toml
```

This TOML file replaced the old pattern of hardcoding analysis values directly in
`cfg/2024ppRef.h`. The header still exists as a compatibility include, but the
actual values are loaded at runtime through:

```
include/AnalysisConfig.h
src/AnalysisConfig.cxx
```

The basic flow is:

```
cfg/2024ppRef.toml
  -> LoadAnalysisConfig(...)
  -> AnalysisConfig cfg
  -> Config()
  -> runAsymmetry / runResiduals / plotResiduals
```

For example, this TOML entry:

```toml
[trigger]
branch = "HLT_AK4PFJet80_v8"
threshold = 100.0
cone = "ak4PF"
```

is read by `LoadAnalysisConfig(...)` and stored as:

```cpp
cfg.hltJ80Branch
cfg.hltJ80Thresh
cfg.trigCone
```

Source files then use:

```cpp
const AnalysisConfig& cfg = Config();
```

and access values with fields like `cfg.coneLabels`, `cfg.jecFilesPerCone`,
`cfg.vetoMapPath`, and `cfg.jetTreePaths`.

The TOML file does not automatically "line up" with the C++ struct by name. The
mapping is explicit in `src/AnalysisConfig.cxx`, where keys such as
`trees.jets`, `paths.veto_map`, and `cones.labels` are assigned to the matching
`AnalysisConfig` fields. When adding a new config value, add it in three places:

1. `cfg/2024ppRef.toml`
2. `include/AnalysisConfig.h`
3. the assignment block inside `LoadAnalysisConfig(...)`

Some arrays are position-matched. In particular, the order of `cones.labels`,
`trees.jets`, and `jec.files` must stay aligned:

```toml
[cones]
labels = [ "ak2PF", "ak3PF", "ak4PF", "ak5PF", "ak6PF" ]

[trees]
jets = [
    "ak2PFJetAnalyzer/t",
    "ak3PFJetAnalyzer/t",
    "ak4PFJetAnalyzer/t",
    "ak5PFJetAnalyzer/t",
    "ak6PFJetAnalyzer/t",
]
```

The loader checks that these arrays have matching lengths before the analysis
runs.

By default, the code reads `cfg/2024ppRef.toml` relative to the directory where
the program is launched. For running from another directory, either set:

```bash
export L2RESIDUALS_HOME=/path/to/L2Residuals-2024ppref
```

or point directly to a config file:

```bash
export L2RESIDUALS_CONFIG=/path/to/2024ppRef.toml
```

<h3> <b> Step 1 </b> — Fill Asymmetry Histograms </h3>

Read HiForest ROOT files, apply L2Relative JEC, select dijets, and fill a 4D {η<sup>probe</sup>, p<sub>T</sub><sup>avg</sup>, α, A} THnSparse for each clustering algorithm.

```
./bin/runAsymmetry <input_HiForest.root> <output_asymmetry.root>
```

Flags for each type of dataset: `-hp`, `-zb`, `-mc`. 
An optional integer fourth argument limits the number of events processed.

<u> Output Structure: </u>

```
hvz_all, hvz, hfilt, h_hlt_j80          # TH1D Event Information
ak4PF/
  ak4PF_asym                            # THnSparse Asymmetries
  ak4PF_incl, ak4PF_tag, ak4PF_probe    # TH3D Kinematics
```

<h3> <b> Step 2 </b> — Extract Residuals </h3>

Read <b>Step 1</b> files after hadd (one data, one MC), project asymmetry distributions for each (η<sup>probe</sup>, p<sub>T</sub><sup>avg</sup>, α) bin, get asymmetry means with different methods (Gaussian, trunc90, trunc95), build response graphs, and extrapolate k<sub>FSR</sub> values as `α → 0`.

```
./bin/runResiduals <data_asymmetry.root> <mc_asymmetry.root> <output_residuals.root>
```

<u> Output Structure: </u>
```
ak4PF/
  QA_data/, QA_mc/                ← A distribution histograms (|η|)
  QA_data_fulleta/, QA_mc_fulleta/
  graphs/                         ← TGraphErrors of R vs α with fit
  Rvals/, Rvals_fulleta/          ← R_data and R_MC TH1Ds per alpha/pT slice
  ak4PF_intercept_abseta_ptavg_40_90_gauss       ← α→0 intercept TH1D
  ak4PF_intercept_abseta_ptavg_40_90_gauss_norm  ← kFSR-normalized variant
  ...  (all methods, all pT slices, abseta and fulleta)
```

<h3> <b> Step 3 </b> — Write Text file </h3>

In |η<sup>probe</sup>| or η<sup>probe</sup> ranges fit correction factors vs p<sub>T</sub><sup>avg</sup> with a 3-parameter function and write a plain text file that can be parsed with a header.

```
./bin/runTextFile <residuals.root> <output.txt> [method] [algorithm]

# methods: gauss (default) | trunc90 | trunc95
# algorithms: ak4PF (default) | ak2PF | ak3PF | ak5PF | ak6PF
```

<h2> Plotting </h2>

`plotResiduals` handles plotting for output files from each step. Flags that don't apply to the input file type skip silently, and all plots are made if no flags are given.

```
./bin/plotResiduals <input.root> <output_dir> [flags]
```

| Flag | Input | Description |
|------|-------|-------------|
| `event` | Step 1 | v<sub>z</sub>, primary vertex filter, HLT trigger |
| `kinematics` | Step 1 | η/φ/p<sub>T</sub> distributions and η-φ maps at different p<sub>T</sub> cuts for incl/tag/probe jets |
| `adist` | Step 2 | Asymmetry distributions for (η<sup>probe</sup>, p<sub>T</sub><sup>avg</sup>, α) bins, Data and MC overlay |
| `roverlay` | Step 2 | Responses for MC and Data vs \|η\| with ratio panel |
| `alpha` | Step 2 | k<sub>FSR</sub> extrapolations |
| `methods` | Step 2 | Gauss vs trunc90 vs trunc95 comparison |
| `etasym` | Step 2 | Full-η vs reflected \|η\| symmetry comparison |
| `normcomp` | Step 2 | normalized vs non-normalized extrapolated corrections comparison |
| `finals` | Step 2 | α→0 intercepts, all p<sub>T</sub><sup>avg</sup> slices overlaid |
| `all` | either | All applicable flags (default) |
<br>

Example of multiple flags being passed space-separated as a single quoted argument:
```bash
./bin/plotResiduals residuals.root plots/ "finals etasym methods"
```

<h2> Condor Submission </h2>

<b>Step 1</b> runs on HTCondor — one job per HiForest input file. <b>Step 2</b> and <b>Step 3</b> run locally after using hadd on the output files from <b>Step 1</b>.

<strong> Before First Submission: </strong>
1. Set `CMSSW_SRC` in `condor/runtime_wrapper.sh`
2. Build the project on LXPLUS: `cmake --build build`

<strong> Submitting HTCondor Jobs </strong>

```bash
# all filelists in data/txt/
bash condor/make_condor.sh /eos/cms/store/group/phys_heavyions/nbarnett/l2residuals [--all]

# specific filelist
bash condor/make_condor.sh /eos/.../l2residuals data/txt/filelist_HiForest_2024ppref_DATA_HP0.txt

# multiple filelists
bash condor/make_condor.sh /eos/.../l2residuals \
    data/txt/filelist_HiForest_2024ppref_DATA_HP0.txt \
    data/txt/filelist_HiForest_2024ppref_DATA_ZB0.txt

# dry run
bash condor/make_condor.sh /eos/.../l2residuals -n
```

Mode is auto-detected from the filelist filename (`*HP*` → `--hard-probes`, `*ZB*` → `--zero-bias`, `*MC*` → `--monte-carlo`). Output is found in `OUTPUT_DIR/condor/asymmetry/<timestamp>/<LABEL>/output_N.root`. Passing `OUTPUT_DIR/condor` or `OUTPUT_DIR/condor/asymmetry` is safe — the script normalizes them to the same path. Working directories and logs go to `condor/submissions/<timestamp>/`. A colored progress bar is displayed as each submission file is generated.

<h2> Hadd Many Files </h2>

```bash
bash condor/batch_hadd.sh \
  /eos/.../l2residuals/asym_hp0.root \
  "/eos/.../l2residuals/HP0/output_*.root" \
  10 2 [z|--zombie-check]
```

`batch_hadd.sh` used hadd on batches of files in parallel using tree reduction, and optionally scanning for zombie/corrupt files. Arguments: `OUT_FILE "IN_FILES" BATCH_SIZE NJOBS [-z]`.

<h2> Data Files </h2>

| Path | Contents |
|------|----------|
| `data/jec/` | L2Relative JEC text files (one for each clustering algorithm, all 5 cones) |
| `data/json/` | Golden JSON |
| `data/veto/` | Jet veto map |
| `data/txt/` | filelists of HiForest files from HardProbes, ZeroBias, and MonteCarlo datasets|
All data collected here is for the pp reference (5.36 TeV) collisions in 2024.

<h2> Tests </h2>

```
ctest --test-dir build
```
Tests exist for `FindLeadingJets` and `MakeDijet` logic, `FoldEtaAxis` correctness, `runTextFile` output format and η ordering, as well as build/library load checks.
