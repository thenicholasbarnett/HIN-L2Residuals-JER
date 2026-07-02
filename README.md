
<h1> L2Residual Jet Energy Corrections </h1>
<br>

<b> Dijet Residuals via p<sub>T</sub>-Balance </b>

Residual jet energy corrections are determined to account for differences in data and simulation by enforcing conservation of transverse momentum. This repository derives L2Residual corrections for the pp reference run collected for the heavy ion collisions in 2024 by CMS. Build this C++ project with CMake to generate standalone binaries and use provided bash scripts to execute them with HTCondor. The versatility of this workflow is in the ability to rebin almost any dimension, overlay different methods of abstraction or calculation, and run any of the steps as compiled executables or interpret the code directly with cling using ROOT.
<br><br>

<b> Step 1 - Find Dijet Asymmetries from HiForest files </b>
<br>
<i> Plot various jet collection kinematics and event information. </i>

<b> Step 2 - Determine Residual Corrections </b>
<br>
<i> Plot asymmetries, k<sub>FSR</sub> extrapolations, data to simulation response ratios, method comparisons, and more! </i>

<b> Step 3 - Write Corrections in Plain Text </b>
<br>
<i> Plot correction factors vs p<sub>T</sub><sup>avg</sup> fitting and extrapolation. </i>
<br><br>

<b>Workflow:</b> \
Fill Dijet Asymmetry Histograms → Hadd Asymmetry Histograms → Extract Residual Correction Values → Write Correction Text File
<br><br>

<h1> Quick Start </h1>

<strong> Clone and Build </strong>

```
git clone git@github.com:thenicholasbarnett/L2Residuals-2024ppref.git
cd L2Residuals-2024ppref
cmake -B build
cmake --build build
```

<strong> Rebuild </strong>

```
rm -rf build bin lib && cmake -B build && cmake --build build
```

<strong> Batch Process Asymmetries </strong>

```
bash ./condor/make_condor.sh <output_files-dir> <input_HiForest-filelist.txt> [CONFIG=cfg/2024ppRef.toml]
```

```
bash ./condor/batch_hadd.sh <output_asymmetry-file.root> <input_glob> <batch_size> <N_parallel>
```

<strong> Get Residual Corrections </strong>

```
./bin/runResiduals <input_data-asymmetries-file.root> <mc_asymmetries_file.root> <output_residuals-file.root> [CONFIG=cfg/2024ppRef.toml]
```

```
./bin/runTextFile <triggered_residuals-file.root> <nontriggered_residuals-file.root> <output_corrections-file.root> <output_text-prefix> [CONFIG=cfg/2024ppRef.toml]
```

<strong> Plot </strong>

```
./bin/runPlotting <input_asymmetry-file.root> <output_plots-dir> [CONFIG=cfg/2024ppRef.toml]
```

```
./bin/runPlotting <input_residuals-file.root> <output_plots-dir> [CONFIG=cfg/2024ppRef.toml]
```

```
./bin/runPlotting <input_corrections-file.root> <output_plots-dir> [CONFIG=cfg/2024ppRef.toml]
```

<h1> Usage </h1> 

<h3> Dependencies </h3>

- ROOT (6.x, with RIO, Tree, Hist, MathCore)
- CMake ≥ 3.10
- C++17
- [nlohmann/json](https://github.com/nlohmann/json) (for applying golden JSON)

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

Analysis parameters — file paths, tree names, cone labels, JEC files, trigger settings, the p<sub>T</sub><sup>avg</sup> slice binning, cut thresholds, and Condor runtime metadata live in `cfg/2024ppRef.toml`. Edit this file and rerun; no recompile needed.

The arrays `cones.labels`, `trees.jets`, and `jec.files` are position-matched and must stay in the same order. The loader validates lengths before running.

`[binning] ptavg_edges` sets the p<sub>T</sub><sup>avg</sup> slice boundaries used by Steps 2/3 and `runPlotting` (strictly ascending, at least 2 entries). Rebinning only needs a TOML edit and a Step 2/3 rerun — no recompile — but it's still a full rerun, not just a replot: Step 2 already collapses the fine-grained pT axis into per-slice fitted means, so the plotting macro only ever sees whatever slicing Step 2 committed to disk. Keep an edge exactly at (or only add edges strictly above/below) `trigger.threshold`, since Step 3 assigns each whole slice to HP or ZB based on its lower edge — a slice straddling the threshold goes entirely to ZB.

`[cuts]` holds `min_jet_pt` (global floor — inclusive jets and the dijet subleading-jet cut; the third jet used for α is exempt), `dphi` (back-to-back requirement), `max_abs_a` (Step 1 dijet acceptance), and `min_entries_per_bin` (Step 2 fit-attempt floor). `[trees] filter` and `[jet_id] veto_map_histogram` name the event-filter branch and which histogram to read from the veto map file.

To add a new C++ analysis config value, update three places: `cfg/2024ppRef.toml`, `include/AnalysisConfig.h`, and the assignment block in `src/AnalysisConfig.cxx`. Condor-only keys such as `[condor].cmssw_src` are consumed by `condor/make_condor.sh`.

Porting to a different collision system or run period: copy `cfg/default.toml` (a commented scaffold, not runnable as-is — every `REQUIRED_SET_ME` placeholder needs filling in) to a new file and pass it explicitly with `CONFIG=path`.

Each compiled entry point accepts an optional `CONFIG=path` argument and prints the resolved TOML path at startup. By default compiled binaries look for `cfg/2024ppRef.toml` under the compiled repo source directory, so they can be launched from outside the repo. Relative paths inside the TOML, such as `data/veto/...`, `data/json/...`, and `data/jec/...`, are resolved relative to that same repo root. Absolute paths are left unchanged. To override globally:

```bash
export L2RESIDUALS_CONFIG=/path/to/2024ppRef.toml  # point to a specific file
export L2RESIDUALS_HOME=/path/to/repo               # looks in $L2RESIDUALS_HOME/cfg/
```

<h3> Naming Convention </h3>

ROOT object names and plot paths are generated through `include/Naming.h`. Treat this header as the naming policy: when a histogram, graph, canvas, or directory needs a name, the code asks `L2Name` what the canonical name should be and writes that name. Later plotting and text-file steps ask the same question again to retrieve the object. If an object does not fit this nomenclature, it is considered outside the current workflow and should be regenerated.

Names are built from ordered keys:

```
cone_object_etaMode_etaBin_ptavg_alpha_methodOrDetail
```

Only keys that apply are included. Examples:

```
ak4PF_A_data_abseta_eta_0p0_0p261_ptavg_30_70_alpha_0p05
ak4PF_R_abseta_eta_0p261_0p522_ptavg_100_175_gauss
ak4PF_intercept_fulleta_ptavg_250_500_trunc90_norm
ak4PF_corrfinal_abseta_gauss_norm
```

The main helper functions are:

| Helper | Purpose |
| :- | :- |
| `L2Name::ObjectName(cone, kind, orderedKeys, detailKeys)` | Assembles the full ROOT object name |
| `L2Name::EtaModeKey(fullEta)` | Returns `abseta` or `fulleta` |
| `L2Name::EtaKey(ieta, fullEta)` | Converts the configured eta-bin edges into keys such as `eta_0p0_0p261` or `eta_m0p261_0p0` |
| `L2Name::PtKey(ptSlice)` | Converts configured p<sub>T</sub><sup>avg</sup> slices into keys such as `ptavg_30_70` |
| `L2Name::AlphaKey(alphaSlice)` | Converts alpha thresholds into keys such as `alpha_0p05` |

Plot directories use the same keys, but joined as path components instead of one long object name. For example:

```
plots/ak4PF/adist/eta_0p0_0p261/ptavg_30_70/alpha_0p05/
```

<h3> <b> Step 1 </b> — Fill Asymmetry Histograms </h3>

Read HiForest ROOT files, apply L2Relative JEC (plus `jec.residual_files`, if set, for data only — never for `-mc`), select dijets, and fill a 4D {η<sup>probe</sup>, p<sub>T</sub><sup>avg</sup>, α, A} THnSparse for each clustering algorithm.

```
./bin/runAsymmetry <input_HiForest-file.root> <output_asymmetry-file.root> [mode] [maxEvents] [CONFIG=cfg/2024ppRef.toml]
```

Flags for each type of dataset: `-t`/`--triggered` (e.g. HardProbes), `-nt`/`--non-triggered` (e.g. ZeroBias or MinBias), `-mc`/`--monte-carlo`. Triggered and non-triggered are both data — the only difference is whether an HLT decision and efficiency-plateau cut apply; which physical dataset plays which role is a per-run-period choice, not something the code hardcodes.
An optional integer fourth argument limits the number of events processed.

<u> Output Structure: </u>

```
hvz_all, hvz, hfilt, h_hlt_j80          # TH1D Event Information
ak4PF/
  ak4PF_asym                            # THnSparse Asymmetries
  ak4PF_incl, ak4PF_tag, ak4PF_probe    # TH3D Kinematics
```

<h3> <b> Step 2 </b> — Extract Residuals </h3>

Read <b>Step 1</b> files after hadd (one data, one MC), project asymmetry distributions for each (η<sup>probe</sup>, p<sub>T</sub><sup>avg</sup>, α) bin, get asymmetry means with different methods (Gaussian, double-Gaussian, trunc90, trunc95), build response graphs, and extrapolate k<sub>FSR</sub> values as `α → 0`.

```
./bin/runResiduals <input_data-asymmetry-file.root> <input_mc-asymmetry-file.root> <output_residuals-file.root> [CONFIG=cfg/2024ppRef.toml]
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

<h3> <b> Step 3 </b> — Write Text Files </h3>

Merges the triggered and non-triggered <b>Step 2</b> outputs (e.g. HardProbes and ZeroBias/MinBias): for each p<sub>T</sub><sup>avg</sup> slice, uses the triggered file if the slice starts at or above the trigger threshold (`trigger.threshold` in `cfg/2024ppRef.toml`, i.e. `cfg.hltJ80Thresh`) and the non-triggered file otherwise — the triggered dataset is trigger-biased below its efficiency plateau, the non-triggered dataset fills in the rest. Processes every cone in `cones.labels` in one call.

```
./bin/runTextFile <triggered_residuals-file.root> <nontriggered_residuals-file.root> <output_corrections-file.root> <output_text-prefix> [method] [direct] [CONFIG=cfg/2024ppRef.toml]

# methods:  gauss (default) | doubleGauss | trunc90 | trunc95
# [direct]: kFSR-normalized intercepts are used by default (the standard method);
#           pass the literal word "direct" for the non-normalized variant instead
```

For a dataset with no triggered/non-triggered split (e.g. one min-bias or single-trigger sample — every pT slice reads from the same file regardless of the trigger threshold):

```
./bin/runTextFile --single <residuals-file.root> <output_corrections-file.root> <output_text-prefix> [method] [direct] [CONFIG=cfg/2024ppRef.toml]
```

For each cone, in |η<sup>probe</sup>| or η<sup>probe</sup> ranges, fits correction factors vs p<sub>T</sub><sup>avg</sup> with a 3-parameter function and writes two plain text files that can be parsed with a header. Since the normalized variant is the default, both filenames get a `_norm` suffix unless `direct` is passed:

```
<output_text-prefix>_<cone>_abseta[_norm].txt   ← fit on |η|, mirrored onto both eta halves
<output_text-prefix>_<cone>_eta[_norm].txt      ← independent fit per full-η bin, no mirroring
```

<u> Output ROOT Structure: </u>
```
ak4PF/
  ak4PF_corrfinal_abseta_gauss[_norm]   ← TH2D: |η| vs p_{T,avg}, z = final merged correction (method/variant used)
  ak4PF_corrfinal_fulleta_gauss[_norm]  ← same, vs full η
  graphs/                                ← TGraphErrors of correction vs p_{T,avg} per eta bin, with the 3-parameter fit embedded
```

<h2> Plotting </h2>

`runPlotting` handles plotting for output files from each step. Flags that don't apply to the input file type skip silently. `flags` has three modes: omitted (empty) runs a curated smart default — NOT every applicable plot, see the table below for which flags that includes per step; `"all"` runs every applicable flag unconditionally; a space-separated list runs exactly those flags.

```
./bin/runPlotting <input_asymmetries-file.root> <output_plots-dir> [flags] [CONFIG=cfg/2024ppRef.toml]
./bin/runPlotting <input_residuals-file.root> <output_plots-dir> [flags] [CONFIG=cfg/2024ppRef.toml]
./bin/runPlotting <input_corrections-file.root> <output_plots-dir> [flags] [CONFIG=cfg/2024ppRef.toml]
```

| Flag | Input | Description | In smart default? |
| :-:| :-: | - | :-: |
| `event` | Step 1 | v<sub>z</sub>, primary vertex filter, HLT trigger | yes |
| `kinematics` | Step 1 | η/φ/p<sub>T</sub> distributions and η-φ maps at different p<sub>T</sub> cuts for incl/tag/probe jets | yes, but tag/probe only — named explicitly (or via `all`) also includes inclusive jets |
| `adist` | Step 2 | Asymmetry distributions for (η<sup>probe</sup>, p<sub>T</sub><sup>avg</sup>, α) bins, Data and MC overlay | yes |
| `roverlay` | Step 2 | Responses for MC and Data vs \|η\| with ratio panel | yes |
| `alpha` | Step 2 | k<sub>FSR</sub> extrapolations | yes |
| `methods` | Step 2 | Gauss vs trunc90 vs trunc95 comparison | no — explicit only |
| `etasym` | Step 2 | Full-η vs reflected \|η\| symmetry comparison | no — explicit only |
| `normcomp` | Step 2 | normalized vs non-normalized extrapolated corrections comparison | no — explicit only |
| `finals` | Step 2 or 3 | α→0 intercepts (Step 2) or final merged corrections (Step 3), all p<sub>T</sub><sup>avg</sup> slices overlaid | Step 3 only (via the corrfinal-grid path); suppressed by default on a Step 2 file even though the flag itself would find data there — pass it explicitly to get it from Step 2 |
| `ptfit` | Step 3 | Correction factor vs p<sub>T</sub><sup>avg</sup> per eta bin, with the 3-parameter fit drawn | yes |
| `all` | any | Every applicable flag, unconditionally (including inclusive-jet kinematics and `finals` from either source) | n/a |

<br>

Example of multiple flags being passed space-separated as a single quoted argument:
```bash
./bin/runPlotting residuals.root plots/ "finals etasym methods" CONFIG=cfg/2024ppRef.toml
```

<h2> Condor Submission </h2>

<b>Step 1</b> runs on HTCondor — one job per HiForest input file. <b>Step 2</b> and <b>Step 3</b> run locally after using hadd on the output files from <b>Step 1</b>.

<strong> Before First Submission: </strong>
1. Set `[condor].cmssw_src` in the TOML you will pass with `CONFIG=...`
2. On LXPLUS, activate CMSSW before configuring/building so the binary links against CMSSW ROOT:

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd <CMSSW_RELEASE>/src
cmsenv
cd /path/to/L2Residuals-2024ppref
cmake -S . -B build
cmake --build build
```

CMake refuses to configure on lxplus without `cmsenv`, and `make_condor.sh` checks both `[condor].cmssw_src` and the built binary's cvmfs ROOT RPATH before submitting. The runtime wrapper copied into each submission is stamped with the TOML's `cmssw_src` value.

<strong> Submitting HTCondor Jobs </strong>

```bash
# all filelists in data/txt/
bash condor/make_condor.sh /eos/cms/store/group/phys_heavyions/nbarnett/l2residuals --all-txt

# specific filelist
bash condor/make_condor.sh /eos/.../l2residuals data/txt/filelist_HiForest_2024ppref_DATA_HP0.txt

# multiple filelists
bash condor/make_condor.sh /eos/.../l2residuals \
    data/txt/filelist_HiForest_2024ppref_DATA_HP0.txt \
    data/txt/filelist_HiForest_2024ppref_DATA_ZB0.txt

# dry run
bash condor/make_condor.sh /eos/.../l2residuals -n

# tagged pass — e.g. a closure rerun using an abs-eta-derived residual file
bash condor/make_condor.sh /eos/.../l2residuals data/txt/filelist_HiForest_2024ppref_DATA_HP0.txt TAG=clos_abseta

# submit with a specific TOML (default: cfg/2024ppRef.toml) — independent of TAG,
# mix and match freely; useful for closure passes or a different run period/system
bash condor/make_condor.sh /eos/.../l2residuals data/txt/filelist_HiForest_2024ppref_DATA_HP0.txt \
    TAG=clos_abseta CONFIG=cfg/2024ppRef_clos_abseta.toml
```

Mode for each filelist is looked up from `[condor.filelist_modes]` in the selected TOML, keyed by the filelist's basename (no `.txt`) — e.g. `filelist_HiForest_2024ppref_DATA_HP0 = "triggered"`. The physical dataset name (HP0, ZB0, MinBias, Jet80, ...) is free-form; only the mapped value (`triggered` | `non-triggered` | `mc`) matters, so this generalizes across run periods/collision systems without touching the script. Comment out (or omit) a filelist's entry in the TOML to skip submitting jobs for it without moving or renaming the file — `make_condor.sh` reports it as `SKIP` and moves on. Output is found in `OUTPUT_DIR/condor/asymmetry/<timestamp>/<LABEL>/output_N.root` — or `OUTPUT_DIR/condor/asymmetry_<TAG>/<timestamp>/...` if `TAG=value` is given, so separate reprocessing/closure passes don't land in the same output tree. Passing `OUTPUT_DIR/condor` or `OUTPUT_DIR/condor/asymmetry[_<TAG>]` is safe — the script normalizes them to the same path. Working directories and logs go to `condor/submissions/<timestamp>/`. A colored progress bar is displayed as each submission file is generated.

`CONFIG=path` picks which TOML gets submitted with the jobs (default `cfg/2024ppRef.toml`), independent of `TAG` — mix and match freely. Whatever's selected is transferred to the sandbox under a fixed name (`analysis_config.toml`), and its `[condor].cmssw_src` value is stamped into the worker `runtime_wrapper.sh`; for a different run period or collision system, point `CONFIG` at that system's TOML (e.g. copied from `cfg/default.toml`) with no other changes needed. The config actually used is echoed at the end of the run and archived alongside the generated `.condor` submission files in `condor/submissions/<timestamp>/`.

<h2> Hadd Many Files </h2>

```bash
bash condor/batch_hadd.sh \
  /eos/.../l2residuals/asym_hp0.root \
  "/eos/.../l2residuals/HP0/output_*.root" \
  10 2 [z|--zombie-check]
```

`batch_hadd.sh` uses hadd on batches of files in parallel using tree reduction, and optionally scanning for zombie/corrupt files. Arguments: `OUT_FILE "IN_FILES" BATCH_SIZE NJOBS [-z]`.

<h2> Data Files </h2>

| Path | Contents |
| :-: | - |
| `data/jec/` | L2Relative JEC text files (one for each clustering algorithm, all 5 cones) |
| `data/json/` | Golden JSON |
| `data/veto/` | Jet veto map |
| `data/txt/` | filelists of HiForest files — triggered (e.g. HardProbes), non-triggered (e.g. ZeroBias/MinBias), and MC datasets |

The data files here are for the pp reference (5.36 TeV) collisions in 2024.

<h2> Tests </h2>

```
ctest --test-dir build
```
Tests exist for `FindLeadingJets` and `MakeDijet` logic, `FoldEtaAxis` correctness, `runTextFile` triggered/non-triggered merge selection and output format/η ordering for both text files, as well as build/library load checks.
