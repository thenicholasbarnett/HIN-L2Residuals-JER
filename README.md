# L2Residuals-2024ppRef

Object-oriented C++17/ROOT framework for CMS L2 Residual jet energy corrections — 2024 pp reference run at CERN/CMS. Rewrite of [DijetResiduals](https://github.com/thenicholasbarnett/DijetResiduals): unified DATA/MC executable, single 4D THnSparse instead of 45 TH2Ds.

L2 Residuals are the final η-dependent data/MC scale factors applied after MC-truth JEC. Measured via dijet pT balance: a barrel tag jet (|η| < 1.3) constrains a probe jet anywhere in η. Asymmetry A = (pT\_probe − pT\_tag) / (pT\_probe + pT\_tag) is extrapolated to zero radiation (α→0), then R = data/MC gives the correction.

---

## Prerequisites

- ROOT 6, CMake ≥ 3.10, C++17
- **nlohmann-json** (for golden JSON in DATA mode):
  - Arch: `sudo pacman -S nlohmann-json`
  - Ubuntu/Debian: `sudo apt install nlohmann-json3-dev`

---

## Build

```bash
git clone https://github.com/thenicholasbarnett/L2Residuals-2024ppref
cd L2Residuals-2024ppref
mkdir build && cd build && cmake .. && make
```

Outputs to `bin/` (executables) and `lib/` (shared library). Run all executables from the **repo root** so `data/`, `cfg/`, and `lib/` resolve correctly.

---

## Pipeline

Status: Steps 1–3 code complete; none yet validated on real data.

### Step 1 — Fill asymmetry histograms

```bash
# DATA (hard probes — default)
./bin/runAsymmetry input_data.root output_data.root

# MC
./bin/runAsymmetry input_mc.root output_mc.root --mc

# zero-bias
./bin/runAsymmetry input.root output.root --zero-bias
```

One job per HiForest input file on Condor (see [Condor](#condor) below). Hadd the outputs offline before Step 2:

```bash
hadd data_hadded.root output_data_*.root
hadd mc_hadded.root   output_mc_*.root
```

**What it does:** For each event across 5 cone sizes (AK2–6 PF), applies L2Relative JEC, finds the leading dijet, enforces quality cuts (see table), and fills one entry in a 4D THnSparse `(η_probe, pT_avg, α, A)`.

| Cut | DATA | MC |
|-----|------|----|
| \|vz\| < 15 cm | ✓ | ✓ |
| Primary vertex filter (`ppvF`) | ✓ | — |
| Golden JSON (runs 387474–387721) | ✓ | — |
| HLT\_AK4PFJet80\_v1, lead pT > 100 GeV | ✓ | — |
| pthat weight on all fills | — | ✓ |

Dijet requirements: Δφ > 2.7, sublead pT > 10 GeV, both jets pass tight jet ID + veto map, at least one in barrel. If third jet fails jet ID → α = 0, treated as clean 2-jet event.

### Step 2 — Extract corrections

```bash
./bin/runResiduals data_hadded.root mc_hadded.root residuals.root
```

Reads the THnSparses from Step 1. For each cone, runs two extraction passes:

- **|η| pass:** folds the full-η sparse to |η| via `FoldEtaAxis`, then for each (α threshold, pT\_avg slice, |η| bin) → project A → extract ⟨A⟩ by 3 methods → R\_data/R\_MC → linear fit vs α → intercept at α→0.
- **Full-η pass:** same loop on the unfolded sparse, produces 36-bin signed-η corrections for the symmetry check.

Three mean estimators: Gaussian fit (A ∈ [−0.5, 0.5]), truncated mean 90% (trim 5% each tail), truncated mean 95% (trim 2.5% each tail). Alpha fit range: thresholds 0.05–0.30 only (above 0.30, trigger bias distorts).

### Step 2 QA — Plots

```bash
./bin/plotResiduals residuals.root            # writes to timestamped plots_residuals_YYYYMMDD/ dir
./bin/plotResiduals residuals.root my_plots   # explicit output dir
```

Writes PNG files per cone:

- **`etasym_{cone}_{method}_{ptavg}.png`** — top: full-η vs |η| reflected; bottom: ratio (flat at 1.0 = symmetric)
- **`methods_{cone}_{ptavg}_abseta.png`** — top: all 3 methods overlaid; bottom: trunc/Gauss ratios (flat = stable)
- **`methods_{cone}_{ptavg}_fulleta.png`** — same for signed η

### Step 3 — Write JEC text file

```bash
./bin/runTextFile residuals.root L2Residuals_2024ppRef_AK4.txt          # default: gauss, ak4PF
./bin/runTextFile residuals.root L2Residuals_2024ppRef_AK4.txt gauss ak4PF
./bin/runTextFile residuals.root L2Residuals_2024ppRef_AK4.txt trunc90 ak4PF
```

Reads intercept TH1Ds from Step 2. For each of 18 |η| bins, fits the pT dependence across 5 slice centers (65, 105, 155, 225, 630 GeV) with `1/(p0 + p1·log10(0.01·pT) + p2/(pT/10))`. Writes a 36-line CMS JEC text file (negative η outer→inner, then positive inner→outer). Failed fits or fewer than 3 valid slices fall back to unity (p0=1, p1=0, p2=0).

Args: `method` = `gauss` | `trunc90` | `trunc95`; `cone` = `ak4PF` | `ak2PF` | `ak3PF` | `ak5PF` | `ak6PF`.

---

## Condor

Step 1 is designed to run one job per input HiForest file on HTCondor at CERN/lxplus.

**Before submitting:**
1. Build the project (see Build above)
2. Stage all 5 cone JEC files into `data/jec/` (currently only AK4 is present — **blocker**)
3. Set `CMSSW_SRC` in `condor/runtime_wrapper.sh`:
   ```bash
   CMSSW_SRC="/afs/cern.ch/user/n/nbarnett/public/condor/workArea/CMSSW_13_2_4/src"
   ```

**Submit:**
```bash
bash condor/make_condor.sh JOBNAME filelist.txt /eos/cms/store/user/nbarnett/output_dir
bash condor/make_condor.sh JOBNAME filelist.txt /eos/path --mc               # MC mode
bash condor/make_condor.sh JOBNAME filelist.txt /eos/path --no-submit        # generate only
```

`filelist.txt` — one absolute path per line to HiForest ROOT files.

**What the script does:**
1. Creates `condor_JOBNAME_TIMESTAMP/` in the current directory
2. Copies `bin/runAsymmetry`, `lib/libl2residuals.so`, `data/`, `runtime_wrapper.sh` into it
3. Generates `submit_JOBNAME.condor` with one queue entry per input file
4. Submits via `condor_submit` (or skips with `--no-submit`)

**Per-job output:** `OUTPUT_DIR/output_JOBNAME_N.root`. Logs in `condor_*/logs/{out,err,log}/job_N.*`.

JobFlavour is `longlunch` (≤2 hours). Each job uses 1 CPU. The worker's `LD_LIBRARY_PATH` is set to `.` so it finds the copied `libl2residuals.so`.

After all jobs finish, hadd the outputs before running Step 2.

---

## Tests

```bash
# Run all tests via CMake (from build/)
cd build && ctest --output-on-failure

# Run individual test binaries directly (from repo root)
./bin/TestDijet
./bin/TestFoldEtaAxis
./bin/TestTextFileWriter
```

| Test binary | Count | What it covers |
|-------------|-------|----------------|
| `TestDijet` | 37 | `FindLeadingJets` (0/1/2/3/4-jet cases, ordering) and `MakeDijet` (barrel assignment, Δφ cut, α computation, pT\_avg, A sign, dphi wrap) |
| `TestFoldEtaAxis` | 7 | Axis properties before/after folding, symmetric fill (content + Sumw2 error), positive-only, negative-only, different-|η| fills stay separate, fold on arbitrary axis index |
| `TestTextFileWriter` | 7 | File structure (header, 36 lines, Npar, pT range), η ordering (negative outer→inner then positive inner→outer), mirror symmetry of fit params, unity fallback for <3 slices, unity fallback for empty bins, fit round-trip (unit input → f≈1.0 at all pT centers), CMS η extent (±5.191) |
| `TestBuild` (bash) | 9 | Executables exist and are runnable, usage message, shared library loadable by ROOT, macro interpretable via `.L` |

All tests exit 0 on success, non-zero on any failure.

---

## Configuration

`cfg/2024ppRef.h` — the only file that needs editing to adapt the analysis:

```
kJECFilesPerCone   — L2Relative text file paths, one list per cone (AK2–6)
kVetoMapPath       — jet veto map ROOT file
kJSONPath          — golden JSON (currently runs 387474–387721)
kHiTreePath        — hiEvtAnalyzer/HiTree
kSkimTreePath      — skimanalysis/HltTree
kTrigTreePath      — hltanalysis/HltTree
kJetTreePaths      — ak{2,3,4,5,6}PFJetAnalyzer/t
kHLTJ80Branch      — HLT_AK4PFJet80_v1  ← confirm version suffix from actual HiForest
kHLTJ80Thresh      — 100.0 GeV (lead pT threshold for trigger efficiency region)
kConeLabels        — {ak2PF, ak3PF, ak4PF, ak5PF, ak6PF}
kTrigCone          — ak4PF  (trigger efficiency evaluated on this cone's kinematics)
```

---

## Remaining Work

| Task | Status |
|------|--------|
| Stage AK2/3/5/6 L2Relative JEC files into `data/jec/` | **Blocker for Condor** — only AK4 present |
| Confirm `kTrigTreePath` / `kHLTJ80Branch` version suffix | Check from actual HiForest file |
| Validate Step 1 on test MC sample | Not done |
| Validate Step 2 on Step 1 output | Not done |
| Validate Step 3 on Step 2 output | Not done |

---

<details>
<summary>Design choices</summary>

**One 4D THnSparse, one fill per event.** The original filled 9 alpha-threshold TH2Ds per event (and a 10th for 2-jet events) — 45 histograms total, non-independent entries, inflated counts. Here one fill at the event's actual alpha goes into `(η_probe, pT_avg, α, A)`. Cumulative alpha slices are recovered at extraction time with `SetRangeUser(alphaAxis, 0, threshold)` — a clean 2-jet event (α=0) falls below every positive threshold automatically.

**Full η stored; |η| folded at extraction.** The sparse keeps signed η across all 36 CMS JEC bins. `FoldEtaAxis()` (Utilities.h) produces a new |η| sparse at analysis time. Both versions are used: |η| for primary corrections, signed η for the forward-backward symmetry check.

**`ExtractAndFit` called twice.** Rather than duplicating the extraction loop, `ResidualsExtractor.cxx` has one static `ExtractAndFit(hData, hMC, ..., etaEdges, nameSuffix)`. `runResiduals` calls it once with |η| sparses + `kAbsEtaEdges`, and once with full-η sparses + `kEtaEdges`. The `nameSuffix` (`""` or `"_fulleta"`) distinguishes all output names.

**Index-based dijet logic.** `FindLeadingJets(corrPt, nref)` returns indices `{lead, sublead, third}` into the corrected-pT array — no copies. JEC changes the pT ordering, so all selections run on corrected pT. Jet ID (needs PF fractions + multiplicities) runs between sorting and dijet construction.

**`eventNumber % 2` for barrel-barrel tag assignment.** CMS event numbers are sequential integers — exact reproducible 50/50 split, independent of job restart or sample order. The original used `TRandom2` seeded at 1, which is sequence-dependent.

**Third-jet ID veto → α = 0.** If the third jet fails jet ID, `hasThird = false` and `alpha = 0`. A fake jet passing pT sorting but failing quality criteria would inflate α, suppressing the event from tight alpha slices. Treating it as a clean 2-jet event is the conservative choice.

**Alpha fit range capped at 0.30.** Only the first 6 of 9 alpha thresholds (0.05–0.30) enter the linear fit. Above 0.30, trigger turn-on bias distorts the A distribution and the linear extrapolation breaks down.

</details>

<details>
<summary>Binning</summary>

| Quantity | Configuration |
|----------|---------------|
| η\_probe axis | 36 CMS JEC bins, −5.191 to 5.191, variable width |
| \|η\_probe\| axis | 18 bins, 0 to 5.191 (folded at extraction) |
| pT\_avg axis | 990 bins × 1 GeV, 10–1000 GeV (coarse slices via SetRangeUser) |
| pT\_avg slices | 40–90, 90–120, 120–190, 190–260, 260–1000 GeV |
| Alpha axis | 50 bins, 0–0.5 |
| Alpha thresholds (fit) | 0.05, 0.10, 0.15, 0.20, 0.25, 0.30 (6 of 9 stored) |
| Asymmetry A | 100 bins, −1 to 1 |
| Gaussian fit window | A ∈ [−0.5, 0.5] |
| Min entries for fit | 100 |

`kAbsEtaEdges`: 0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191

Step 3 pT slice centers used for the pT-dependence fit: 65, 105, 155, 225, 630 GeV.

</details>

<details>
<summary>Step 2 output structure</summary>

For each cone (e.g. `ak4PF`) in `residuals.root`:

| Object | Description |
|--------|-------------|
| `ak4PF_asym_data_abseta` | Folded \|η\| THnSparse from DATA |
| `ak4PF_asym_mc_abseta` | Folded \|η\| THnSparse from MC |
| `ak4PF_QA_data/` | 1D A distributions, \|η\| pass |
| `ak4PF_QA_mc/` | same for MC |
| `ak4PF_QA_data_fulleta/` | 1D A distributions, full-η pass |
| `ak4PF_QA_mc_fulleta/` | same for MC |
| `ak4PF_graphs/` | TGraphErrors of R\_data/R\_MC vs α threshold + linear fit |
| `ak4PF_intercept_{method}{ptavg}` | TH1D of α→0 corrections vs \|η\|, 18 bins |
| `ak4PF_intercept_{method}{ptavg}_fulleta` | same vs signed η, 36 bins |

`{method}` ∈ `gauss`, `trunc90`, `trunc95`  
`{ptavg}` ∈ `_ptavg_40_90`, `_ptavg_90_120`, `_ptavg_120_190`, `_ptavg_190_260`, `_ptavg_260_1000`

Total per cone: 15 intercept TH1Ds for |η| + 15 for full η + QA directories + graphs.

</details>

<details>
<summary>File map</summary>

```
cfg/
  2024ppRef.h           JEC paths, trigger, TTree paths, cone labels

include/
  RunAsymmetry.h        Step 1 entry point declaration
  ResidualsExtractor.h  Step 2 entry point declaration
  TextFileWriter.h      Step 3 entry point declaration
  Binning.h             BinningConfig — all axis defs, pT/alpha slices, eta edges
  Dijet.h               FindLeadingJets (O(n) single-pass), MakeDijet — header-only
  DijetHistograms.h     ConeHistograms struct — owns the 4D THnSparse + control TH3Ds
  Utilities.h           MakeTHnSparse, FoldEtaAxis, canvas/style helpers
  JetCorrector.h        Yi Chen v3.0 — chains text-file JEC levels
  JetSelection.h        JetSelect — tight pp jet ID (η-dependent) + veto map
  JSON_handler.h        golden JSON checking (nlohmann-json)
  BranchMapping.h       SetBranches() — one-call TTree branch binding
  EventStructs.h        EventStruct (vz, weight, run/lumi), FiltersStruct
  JetStruct.h           JetStruct<MAXNREF> — reco/ref/gen arrays, PF info
  ProgressBar.h         terminal progress bar
  Colors.h              Hiroshige + Klimt color palettes

src/
  RunAsymmetry.cxx       Step 1 implementation
  ResidualsExtractor.cxx Step 2 — ExtractAndFit called twice per cone
  TextFileWriter.cxx     Step 3 — pT-dependence fit + JEC text file writer

macros/
  runAsymmetry.C    Step 1 entry point (compiled or ROOT interpreted)
  runResiduals.C    Step 2 entry point
  plotResiduals.C   Step 2 QA — eta symmetry + method comparison plots
  runTextFile.C     Step 3 entry point

tests/
  TestDijet.cxx          37 tests — FindLeadingJets + MakeDijet
  TestFoldEtaAxis.cxx     7 tests — FoldEtaAxis
  TestTextFileWriter.cxx  7 tests — JEC text file structure, ordering, fallbacks
  TestBuild.sh            9 checks — build artifacts + ROOT interpretability

data/
  jec/    L2Relative text files (only AK4 staged — AK2/3/5/6 needed)
  json/   Cert_Collisions2024_ppref_387474_387721_golden.json
  veto/   Winter25Prompt25_RunCDE.root (jet veto map)

condor/
  make_condor.sh      generates + submits Condor jobs for Step 1
  runtime_wrapper.sh  executes on worker nodes — set CMSSW_SRC before submitting
```

Every macro can also be interpreted by ROOT (without recompiling):
```bash
root -l -b -q 'macros/runAsymmetry.C("in.root","out.root","--mc")'
```
Build the library first; run from repo root.

</details>

<details>
<summary>Adapting for a new run period</summary>

Edit `cfg/2024ppRef.h` (or create a new config header under `cfg/` and swap the `#include`):

1. Update `kJECFilesPerCone` with new L2Relative text file paths for each cone
2. Update `kVetoMapPath` with the new veto map ROOT file
3. Update `kJSONPath` with the new golden JSON
4. Update `kTrigTreePath`, `kHLTJ80Branch`, `kHLTJ80Thresh` for the new trigger
5. Confirm `kHiTreePath` and `kSkimTreePath` match the new HiForest schema

All analysis logic, histogram structure, and extraction code stays unchanged.

</details>

---

Nicholas Barnett — thenicholasbarnett@gmail.com
