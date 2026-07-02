
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

<strong> Local CMake Build </strong>

Use this on a laptop or anywhere you want standalone binaries outside CMSSW:

```bash
git clone git@github.com:thenicholasbarnett/L2Residuals-2024ppref.git
cd L2Residuals-2024ppref
cmake -B build
cmake --build build
```

CMake writes executables to `build/bin/` and the shared library to `build/lib/`. Rebuild from scratch with only:

```bash
rm -rf build
cmake -B build
cmake --build build
```

<strong> LXPLUS SCRAM Build </strong>

Use this when you want `scram b -j4` and CMSSW/Condor compatibility. The package path is currently fixed:

```bash
cd <CMSSW_RELEASE>/src
cmsenv
mkdir -p Analysis
git clone git@github.com:thenicholasbarnett/L2Residuals-2024ppref.git Analysis/L2Residuals
scram b -j4
which runAsymmetry
```

SCRAM writes executables to `$CMSSW_BASE/bin/$SCRAM_ARCH/`, which `cmsenv` puts on `PATH`, so commands are called as `runAsymmetry`, `runResiduals`, etc. If `cmsenv` is not available in a fresh shell, first run `source /cvmfs/cms.cern.ch/cmsset_default.sh`.

<strong> Command Tokens </strong>

All compiled commands use `KEY=value` tokens. There are no positional arguments. `CONFIG=path` is required by every compiled entry point, and typos such as `OUTPT=...` are rejected immediately.

```bash
# Step 1, usually on lxplus/Condor or from a SCRAM shell
runAsymmetry INPUT=<input_HiForest.root> OUTPUT=<output_asymmetry.root> MODE=triggered|non-triggered|mc CONFIG=cfg/2024ppRef.toml

# Step 2, often local CMake after hadd
./build/bin/runResiduals DATA=<data_asymmetry.root> MC=<mc_asymmetry.root> OUTPUT=<residuals.root> CONFIG=cfg/2024ppRef.toml

# Step 3, split triggered/non-triggered residuals
./build/bin/runTextFile TRIGGERED=<triggered_residuals.root> NONTRIGGERED=<nontriggered_residuals.root> OUTPUT=<corrections.root> [PREFIX=<name>] CONFIG=cfg/2024ppRef.toml

# Step 3, one residual file for every pT slice
./build/bin/runTextFile SINGLE=<residuals.root> OUTPUT=<corrections.root> [PREFIX=<name>] CONFIG=cfg/2024ppRef.toml

# JES/JER extraction, on a hadded Step 1 MC file (independent of Steps 2/3)
./build/bin/runResponse INPUT=<mc_asymmetry.root> OUTPUT=<response.root> CONFIG=cfg/2024ppRef.toml

# Plot any workflow output
./build/bin/runPlotting INPUT=<input.root> OUTDIR=<output_plots_dir> CONFIG=cfg/2024ppRef.toml
```

<strong> Batch Process Asymmetries </strong>

Build with SCRAM or with CMake under `cmsenv`, set `[condor].cmssw_src` in the selected TOML, then submit Step 1 jobs:

```bash
bash ./condor/make_condor.sh <output_files_dir> <input_HiForest_filelist.txt> CONFIG=cfg/2024ppRef.toml
```

Merge Step 1 outputs before Step 2:

```bash
bash ./condor/batch_hadd.sh <output_asymmetry.root> <input_glob> <batch_size> <N_parallel>
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

For CMSSW/SCRAM builds on lxplus, place the checkout at `$CMSSW_BASE/src/Analysis/L2Residuals` and build from `src`:

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd <CMSSW_RELEASE>/src
cmsenv
mkdir -p Analysis
git clone git@github.com:thenicholasbarnett/L2Residuals-2024ppref.git Analysis/L2Residuals
scram b -j4
```

SCRAM builds the same library code from `src/*.cxx` and the executable wrappers in `bin/*.cc`. The SCRAM package name is currently hard-coded as `Analysis/L2Residuals`, so keep the checkout at `$CMSSW_BASE/src/Analysis/L2Residuals` for now. The CMake workflow remains the preferred standalone/laptop build; SCRAM is the CMSSW/lxplus build front door.

TODO: figure out whether this can be made package-path independent without losing the clean SCRAM library/executable split.

CMake binaries are written to `build/bin/`, the shared library to `build/lib/`, and all generated build files stay under `build/`. The source `bin/` directory is reserved for SCRAM executable wrappers.

<h3> Configuration </h3>

Analysis parameters — file paths, tree names, cone labels, JEC files, trigger settings, the p<sub>T</sub><sup>avg</sup> slice binning, cut thresholds, and Condor runtime metadata live in `cfg/2024ppRef.toml`. Edit this file and rerun; no recompile needed.

The arrays `cones.labels`, `trees.jets`, and `jec.files` are position-matched and must stay in the same order. The loader validates lengths before running.

`[binning] ptavg_edges` sets the p<sub>T</sub><sup>avg</sup> slice boundaries used by Steps 2/3 and `runPlotting` (strictly ascending, at least 2 entries). Rebinning only needs a TOML edit and a Step 2/3 rerun — no recompile — but it's still a full rerun, not just a replot: Step 2 already collapses the fine-grained pT axis into per-slice fitted means, so the plotting macro only ever sees whatever slicing Step 2 committed to disk. Keep an edge exactly at (or only add edges strictly above/below) `trigger.threshold`, since Step 3 assigns each whole slice to the triggered or non-triggered dataset based on its lower edge — a slice straddling the threshold goes entirely to the non-triggered dataset.

`[cuts]` holds `min_jet_pt` (global floor — inclusive jets and the dijet subleading-jet cut; the third jet used for α is exempt), `dphi` (back-to-back requirement), `max_abs_a` (Step 1 dijet acceptance), and `min_entries_per_bin` (Step 2 fit-attempt floor). `[trees] filter` and `[jet_id] veto_map_histogram` name the event-filter branch and which histogram to read from the veto map file.

To add a new C++ analysis config value, update three places: `cfg/2024ppRef.toml`, `include/AnalysisConfig.h`, and the assignment block in `src/AnalysisConfig.cxx`. Condor-only keys such as `[condor].cmssw_src` are consumed by `condor/make_condor.sh`.

Porting to a different collision system or run period: copy `cfg/default.toml` (a commented scaffold, not runnable as-is — every `REQUIRED_SET_ME` placeholder needs filling in) to a new file and pass it explicitly with `CONFIG=path`.

Every compiled entry point takes `CONFIG=path` as a **required** argument and prints the resolved TOML path at startup — there is no implicit default TOML. Which config gets loaded is a physics-affecting choice (e.g. the main run-period config vs. a closure config with different `residual_files`), so a missing `CONFIG=` is always a hard CLI error, never a silent fallback to `cfg/2024ppRef.toml`. Relative paths inside the TOML, such as `data/veto/...`, `data/json/...`, and `data/jec/...`, are resolved relative to the repo root — the compiled-in source directory, or `L2RESIDUALS_HOME` if set. `L2RESIDUALS_HOME` is a separate, narrower mechanism than `CONFIG`: it only affects *where relative paths inside an already-chosen TOML resolve from*, and does not itself select or default a TOML.

`L2RESIDUALS_CONFIG` is the environment-variable equivalent of `CONFIG=path` — set it once in a shell so interpreted-ROOT usage (which bypasses the compiled CLI entirely, see below) and repeated invocations don't need `CONFIG=...` retyped every time. A compiled binary's own `CONFIG=` token still takes precedence when both are present:

```bash
export L2RESIDUALS_CONFIG=/path/to/2024ppRef.toml   # required by every entry point, one way or another
export L2RESIDUALS_HOME=/path/to/repo                # relative-path resolution only; does not select a TOML
```

<h3> Command-Line Convention </h3>

Every compiled binary's arguments are `KEY=value` tokens — there are no positional arguments anywhere in this CLI, and argument order never matters. This removes an entire class of silent argument-order bugs (e.g. swapping `DATA=`/`MC=`, or `TRIGGERED=`/`NONTRIGGERED=`), since every value is self-identifying regardless of where it appears on the command line. Required tokens — `CONFIG=` on every binary, plus whichever input/output tokens that binary needs — are validated immediately: a missing required token, an unrecognized token (e.g. a typo like `OUTPT=`), or any bare argument that isn't a valid `KEY=value` token at all is an immediate usage error, never something silently ignored or defaulted. See `include/CliTokens.h` for the parser shared by every entry point.

| Binary | Required tokens | Optional tokens |
| :- | :- | :- |
| `runAsymmetry` | `INPUT=`, `OUTPUT=`, `MODE=`, `CONFIG=` | `MAXEVENTS=` |
| `runResiduals` | `DATA=`, `MC=`, `OUTPUT=`, `CONFIG=` | none |
| `runTextFile` split mode | `TRIGGERED=`, `NONTRIGGERED=`, `OUTPUT=`, `CONFIG=` | `PREFIX=`, `METHOD=`, `NORM=` |
| `runTextFile` single-file mode | `SINGLE=`, `OUTPUT=`, `CONFIG=` | `PREFIX=`, `METHOD=`, `NORM=` |
| `runResponse` | `INPUT=`, `OUTPUT=`, `CONFIG=` | none |
| `runPlotting` | `INPUT=`, `CONFIG=` | `OUTDIR=`, `FLAGS=`, `CLOSURE=` |

`MODE=` must be `triggered`, `non-triggered`, or `mc`. `PREFIX=` for `runTextFile` is a plain filename prefix, not a path — it must not contain `/`, and defaults to `L2Residual` when omitted; see the Step 3 section below for where the resulting text files go. `METHOD=` may be `gauss`, `doubleGauss`, `trunc90`, or `trunc95`. `NORM=false` selects the direct non-normalized Step 3 variant; omitted or any value other than `false` uses the normalized standard variant. `FLAGS=` for plotting must be quoted if it contains spaces, e.g. `FLAGS="finals etasym methods"`. `CLOSURE=true` for `runPlotting` fixes the `finals` plot's y-axis to 0.95–1.05 with red dotted guide lines at 0.99/1.01, for checking a closure pass's R<sub>MC</sub>/R<sub>data</sub> ≈ 1 at a glance; omitted (or any other value) keeps the normal auto-scaled range used for the correction derivation itself. No effect on any flag other than `finals`.

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

Read HiForest ROOT files, apply L2Relative JEC (plus `jec.residual_files`, if set, for data only — never for `MODE=mc`), select dijets, and fill a 4D {η<sup>probe</sup>, p<sub>T</sub><sup>avg</sup>, α, A} THnSparse for each clustering algorithm.

```
./build/bin/runAsymmetry INPUT=<input_HiForest-file.root> OUTPUT=<output_asymmetry-file.root> MODE=triggered|non-triggered|mc [MAXEVENTS=n] CONFIG=cfg/2024ppRef.toml
```

`MODE=triggered` (e.g. HardProbes), `MODE=non-triggered` (e.g. ZeroBias or MinBias), or `MODE=mc`. Triggered and non-triggered are both data — the only difference is whether an HLT decision and efficiency-plateau cut apply; which physical dataset plays which role is a per-run-period choice, not something the code hardcodes.
`MAXEVENTS=n` optionally limits the number of events processed.

<u> Output Structure: </u>

```
hvz_all, hvz, hfilt, h_hlt_j80          # TH1D Event Information
ak4PF/
  ak4PF_asym                            # THnSparse Asymmetries
  ak4PF_incl, ak4PF_tag, ak4PF_probe    # TH3D Kinematics
  ak4PF_incl_resp, ak4PF_tag_resp,
    ak4PF_probe_resp                    # THnSparse (η_reco, p_T^reco, η_gen, p_T^gen, response), MC mode only
```

<h3> <b> JES/JER Extraction </b> — <code>runResponse</code> </h3>

Reads a hadded Step 1 MC file (no data-mode equivalent — the response THnSparses above are MC-only) and, for each cone and each matched jet collection (incl/tag/probe), extracts JES (mean of a Gaussian fit to the p<sub>T</sub><sup>reco</sup>/p<sub>T</sub><sup>gen</sup> response) and JER (σ/mean of the same fit, the standard fractional-resolution convention) as a function of η<sub>gen</sub> and, separately, p<sub>T</sub><sup>gen</sup>. Binning is always on the gen quantity, never reco — conditioning on truth directly avoids the falling-spectrum migration bias that binning by a resolution-smeared reco quantity would introduce.

```
./build/bin/runResponse INPUT=<input_mc-asymmetry-file.root> OUTPUT=<output_response-file.root> CONFIG=cfg/2024ppRef.toml
```

<u> Output Structure: </u>
```
ak4PF/
  QA_response_abseta/, QA_response_fulleta/, QA_response_ptgen/  ← raw per-bin response distributions (no embedded fit)
  ak4PF_JES_abseta_vs_etagen_incl, _tag, _probe    ← TH1D vs |η_gen|
  ak4PF_JER_abseta_vs_etagen_incl, _tag, _probe
  ak4PF_JES_fulleta_vs_etagen_incl, _tag, _probe   ← TH1D vs η_gen
  ak4PF_JER_fulleta_vs_etagen_incl, _tag, _probe
  ak4PF_JES_vs_ptgen_incl, _tag, _probe            ← TH1D vs p_T^gen
  ak4PF_JER_vs_ptgen_incl, _tag, _probe
```

<h3> <b> Step 2 </b> — Extract Residuals </h3>

Read <b>Step 1</b> files after hadd (one data, one MC), project asymmetry distributions for each (η<sup>probe</sup>, p<sub>T</sub><sup>avg</sup>, α) bin, get asymmetry means with different methods (Gaussian, double-Gaussian, trunc90, trunc95), build response graphs, and extrapolate k<sub>FSR</sub> values as `α → 0`.

```
./build/bin/runResiduals DATA=<input_data-asymmetry-file.root> MC=<input_mc-asymmetry-file.root> OUTPUT=<output_residuals-file.root> CONFIG=cfg/2024ppRef.toml
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
./build/bin/runTextFile TRIGGERED=<triggered_residuals-file.root> NONTRIGGERED=<nontriggered_residuals-file.root> OUTPUT=<output_corrections-file.root> [PREFIX=<name>] [METHOD=gauss] [NORM=true] CONFIG=cfg/2024ppRef.toml

# METHOD: gauss (default) | doubleGauss | trunc90 | trunc95
# NORM:   true (default) uses the kFSR-normalized intercepts (the standard method);
#         NORM=false uses the direct, non-normalized variant instead
```

For a dataset with no triggered/non-triggered split (e.g. one min-bias or single-trigger sample — every pT slice reads from the same file regardless of the trigger threshold), pass `SINGLE=` instead of `TRIGGERED=`/`NONTRIGGERED=` (mutually exclusive with them):

```
./build/bin/runTextFile SINGLE=<residuals-file.root> OUTPUT=<output_corrections-file.root> [PREFIX=<name>] [METHOD=gauss] [NORM=true] CONFIG=cfg/2024ppRef.toml
```

For each cone, in |η<sup>probe</sup>| or η<sup>probe</sup> ranges, fits correction factors vs p<sub>T</sub><sup>avg</sup> with a 3-parameter function and writes two plain text files per cone that can be parsed with a header. These always go to `data/jec/preliminary/` (relative to the repo root, created automatically if missing) — a gitignored directory reserved for locally-generated/preliminary corrections, not a caller-chosen path. `PREFIX=` is a plain filename prefix, not a path (it must not contain `/`); it defaults to `L2Residual` when omitted. Since the normalized variant is the default, both filenames get a `_norm` suffix unless `NORM=false` is passed:

```
data/jec/preliminary/<prefix>_<cone>_abseta[_norm].txt   ← fit on |η|, mirrored onto both eta halves
data/jec/preliminary/<prefix>_<cone>_eta[_norm].txt      ← independent fit per full-η bin, no mirroring
```

<u> Output ROOT Structure: </u>
```
ak4PF/
  ak4PF_corrfinal_abseta_gauss[_norm]   ← TH2D: |η| vs p_{T,avg}, z = final merged correction (method/variant used)
  ak4PF_corrfinal_fulleta_gauss[_norm]  ← same, vs full η
  graphs/                                ← TGraphErrors of correction vs p_{T,avg} per eta bin, with the 3-parameter fit embedded
```

<h2> Plotting </h2>

`runPlotting` handles plotting for output files from each step. Flags that don't apply to the input file type skip silently. `FLAGS=` has three modes: omitted (empty) runs a curated smart default — NOT every applicable plot, see the table below for which flags that includes per step; `FLAGS=all` runs every applicable flag unconditionally; a space-separated value runs exactly those flags.

```
./build/bin/runPlotting INPUT=<input_asymmetries-file.root> [OUTDIR=<output_plots-dir>] [FLAGS="..."] CONFIG=cfg/2024ppRef.toml
./build/bin/runPlotting INPUT=<input_residuals-file.root> [OUTDIR=<output_plots-dir>] [FLAGS="..."] [CLOSURE=true] CONFIG=cfg/2024ppRef.toml
./build/bin/runPlotting INPUT=<input_corrections-file.root> [OUTDIR=<output_plots-dir>] [FLAGS="..."] [CLOSURE=true] CONFIG=cfg/2024ppRef.toml
./build/bin/runPlotting INPUT=<input_response-file.root> [OUTDIR=<output_plots-dir>] [FLAGS="..."] CONFIG=cfg/2024ppRef.toml
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
| `finals` | Step 2 or 3 | α→0 intercepts (Step 2) or final merged corrections (Step 3), all p<sub>T</sub><sup>avg</sup> slices overlaid. `CLOSURE=true` fixes the y-range to 0.95–1.05 with 0.99/1.01 guide lines for closure checks | Step 3 only (via the corrfinal-grid path); suppressed by default on a Step 2 file even though the flag itself would find data there — pass it explicitly to get it from Step 2 |
| `ptfit` | Step 3 | Correction factor vs p<sub>T</sub><sup>avg</sup> per eta bin, with the 3-parameter fit drawn | yes |
| `response` | `runResponse` | Per-bin p<sub>T</sub><sup>reco</sup>/p<sub>T</sub><sup>gen</sup> distributions with a Gaussian guide fit, plus JES/JER vs η<sub>gen</sub> and vs p<sub>T</sub><sup>gen</sup> summary overlays (incl/tag/probe) | yes |
| `all` | any | Every applicable flag, unconditionally (including inclusive-jet kinematics and `finals` from either source) | n/a |

<br>

Example of multiple flags being passed space-separated as a single quoted `FLAGS=` value:
```bash
./build/bin/runPlotting INPUT=residuals.root OUTDIR=plots/ FLAGS="finals etasym methods" CONFIG=cfg/2024ppRef.toml
```

<h2> Condor Submission </h2>

<b>Step 1</b> runs on HTCondor — one job per HiForest input file. <b>Step 2</b> and <b>Step 3</b> run locally after using hadd on the output files from <b>Step 1</b>.

<strong> Before First Submission: </strong>
1. Set `[condor].cmssw_src` in the TOML you will pass with `CONFIG=...`; this should point to the CMSSW `src` directory used by worker jobs.
2. Build `runAsymmetry` in a CMSSW environment. SCRAM is preferred on lxplus:

```bash
cd <CMSSW_RELEASE>/src
cmsenv
cd Analysis/L2Residuals
scram b -j4
which runAsymmetry
```

`make_condor.sh` first looks for the SCRAM-built executable at `$CMSSW_BASE/bin/$SCRAM_ARCH/runAsymmetry`. If it exists, that binary is copied into the Condor sandbox. If no SCRAM binary exists, the script falls back to the CMake binary at `build/bin/runAsymmetry` and copies `build/lib/libl2residuals.so` with it. For that fallback, configure/build under `cmsenv` so the binary links against CMSSW ROOT:

```bash
cd <CMSSW_RELEASE>/src
cmsenv
cd Analysis/L2Residuals
cmake -S . -B build
cmake --build build
```

The runtime wrapper copied into each submission is stamped with the TOML's `cmssw_src` value and runs `scramv1 runtime -sh` on the worker before executing `./runAsymmetry`.

<strong> Submitting HTCondor Jobs </strong>

```bash
# all filelists in data/txt/
bash condor/make_condor.sh /eos/cms/store/group/phys_heavyions/nbarnett/l2residuals --all-txt CONFIG=cfg/2024ppRef.toml

# specific filelist
bash condor/make_condor.sh /eos/.../l2residuals data/txt/filelist_HiForest_2024ppref_DATA_HP0.txt CONFIG=cfg/2024ppRef.toml

# multiple filelists
bash condor/make_condor.sh /eos/.../l2residuals \
    data/txt/filelist_HiForest_2024ppref_DATA_HP0.txt \
    data/txt/filelist_HiForest_2024ppref_DATA_ZB0.txt \
    CONFIG=cfg/2024ppRef.toml

# dry run
bash condor/make_condor.sh /eos/.../l2residuals --all-txt CONFIG=cfg/2024ppRef.toml -n

# tagged pass — e.g. a closure rerun using an abs-eta-derived residual file
bash condor/make_condor.sh /eos/.../l2residuals data/txt/filelist_HiForest_2024ppref_DATA_HP0.txt \
    CONFIG=cfg/2024ppRef.toml TAG=clos_abseta

# a different TOML entirely — independent of TAG, mix and match freely;
# useful for closure passes or a different run period/system
bash condor/make_condor.sh /eos/.../l2residuals data/txt/filelist_HiForest_2024ppref_DATA_HP0.txt \
    TAG=clos_abseta CONFIG=cfg/2024ppRef_clos_abseta.toml
```

Mode for each filelist is looked up from `[condor.filelist_modes]` in the selected TOML, keyed by the filelist's basename (no `.txt`) — e.g. `filelist_HiForest_2024ppref_DATA_HP0 = "triggered"`. The physical dataset name (HP0, ZB0, MinBias, Jet80, ...) is free-form; only the mapped value (`triggered` | `non-triggered` | `mc`) matters, so this generalizes across run periods/collision systems without touching the script. Comment out (or omit) a filelist's entry in the TOML to skip submitting jobs for it without moving or renaming the file — `make_condor.sh` reports it as `SKIP` and moves on. Output is found in `OUTPUT_DIR/condor/asymmetry/<timestamp>/<LABEL>/output_N.root` — or `OUTPUT_DIR/condor/asymmetry_<TAG>/<timestamp>/...` if `TAG=value` is given, so separate reprocessing/closure passes don't land in the same output tree. Passing `OUTPUT_DIR/condor` or `OUTPUT_DIR/condor/asymmetry[_<TAG>]` is safe — the script normalizes them to the same path. Working directories and logs go to `condor/submissions/<timestamp>/`. A colored progress bar is displayed as each submission file is generated.

`CONFIG=path` picks which TOML gets submitted with the jobs — **required, no default**, since which TOML gets used is a physics-affecting choice (e.g. main run-period config vs. a closure config) — independent of `TAG`, mix and match freely. Whatever's selected is transferred to the sandbox under a fixed name (`analysis_config.toml`), and its `[condor].cmssw_src` value is stamped into the worker `runtime_wrapper.sh`; for a different run period or collision system, point `CONFIG` at that system's TOML (e.g. copied from `cfg/default.toml`) with no other changes needed. The config actually used is echoed at the end of the run and archived alongside the generated `.condor` submission files in `condor/submissions/<timestamp>/`.

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
| `data/jec/preliminary/` | Step 3 output — locally-generated L2Residual correction text files (gitignored; created automatically by `runTextFile`) |
| `data/json/` | Golden JSON |
| `data/veto/` | Jet veto map |
| `data/txt/` | filelists of HiForest files — triggered (e.g. HardProbes), non-triggered (e.g. ZeroBias/MinBias), and MC datasets |

The data files here are for the pp reference (5.36 TeV) collisions in 2024.

<h2> Tests </h2>

```
ctest --test-dir build
```
Tests exist for `FindLeadingJets` and `MakeDijet` logic, `FoldEtaAxis` correctness, `ConeHistograms`' MC response histogram plumbing/matching, `runTextFile` triggered/non-triggered merge selection and output format/η ordering for both text files, `runResponse`'s JES/JER extraction (Gaussian fit recovery, per-bin entry-count guard, vs-η<sub>gen</sub>/vs-p<sub>T</sub><sup>gen</sup> binning), as well as build/library load checks.
