# L2Residuals-2024ppRef

> Object-oriented C++/ROOT framework for CMS L2 Residual jet energy corrections — 2024 pp reference run

---

## What This Does

Jet energies measured by CMS need several layers of calibration. This framework computes the **L2 Residual** correction: the final eta-dependent scale factor that brings data into agreement with simulation after MC-truth corrections have already been applied.

The measurement uses **dijet pT balance**. In a back-to-back dijet event, a well-measured barrel jet (|η| < 1.3) constrains the response of a probe jet anywhere in η. Comparing the asymmetry A = (pT_probe − pT_tag) / (pT_probe + pT_tag) between data and MC across millions of events gives the correction as a function of probe η and pT. The output is a text file in the standard CMS JEC format, ready to be chained into `JetCorrector`.

This is a rewrite of the procedural [DijetResiduals](https://github.com/thenicholasbarnett/DijetResiduals) codebase. The two separate DATA and MC scripts have been unified into a single executable, and the 9×5 array of 2D histograms has been replaced with a single 4D sparse histogram.

---

## Quick Start

**Requirements:** ROOT 6, CMake ≥ 3.10, C++17

For the golden JSON (DATA): `sudo pacman -S nlohmann-json` (Arch) or install `nlohmann-json3-dev` (Debian/Ubuntu).

```bash
git clone https://github.com/thenicholasbarnett/L2Residuals-2024ppref
cd L2Residuals-2024ppref
mkdir build && cd build
cmake ..
make
```

**Run the analysis:**

```bash
# Step 1 — fill asymmetry histograms (one job per input file on Condor; hadd outputs offline)
./bin/runAsymmetry input_data.root output.root
./bin/runAsymmetry input_mc.root   output.root --mc

# Step 2 — extract L2 Residual corrections (both |eta| and full eta)
./bin/runResiduals data_hadded.root mc_hadded.root residuals.root

# Step 2 QA — plot eta symmetry and method comparisons
./bin/plotResiduals residuals.root           # timestamped output dir
./bin/plotResiduals residuals.root my_plots  # explicit output dir

# Step 3 (not started) — write JEC-format text file
./bin/runTextFile residuals.root L2Residuals_2024ppRef.txt
```

**Run unit tests:**

```bash
# from the build directory
ctest --output-on-failure

# or directly, without CMake (workaround for vdt packaging issue on Arch)
g++ -std=c++17 $(root-config --cflags) tests/TestDijet.cxx -o /tmp/TestDijet $(root-config --libs) -I include -I cfg
/tmp/TestDijet
```

---

## How It Works

### Step 1 — Asymmetry Generation

For each event, raw jet pT values are corrected offline using `JetCorrector` (Yi Chen v3.0, chaining L2Relative + L2Residual text files). Jets are ranked by corrected pT in a single pass.

Events must satisfy:

| Cut | DATA | MC |
|-----|------|----|
| \|vz\| < 15 cm | ✓ | ✓ |
| Primary vertex filter (`ppvF`) | ✓ | — |
| Golden JSON (certified lumi sections) | ✓ | — |
| Trigger efficiency region | ✓ | — |
| pthat weight on all fills | — | ✓ |

The dijet is selected by requiring δφ > 2.7, at least one jet in the barrel, and sublead pT > 10 GeV. Both leading and subleading jets must pass tight jet ID and the jet veto map. The tag jet is the barrel jet; if both are in the barrel, `eventNumber % 2` gives a deterministic 50/50 assignment.

The result — probe η, pT_avg, alpha = pT_third/pT_avg, and asymmetry A — is filled once into a 4D THnSparse. At extraction time, cumulative alpha slices are recovered with `SetRangeUser`.

### Step 2 — Residuals Extraction *(code complete, not yet validated on data)*

`runResiduals` takes hadded DATA and MC asymmetry files from Step 1. For each cone size it runs two complete extraction passes:

**|η| pass** — folds the full-eta THnSparse to |η| once via `FoldEtaAxis`, then for each (alpha_threshold, pT_avg_slice, |η| bin):
1. `SetRangeUser` on alpha and pT_avg axes
2. Project to 1D asymmetry distribution
3. Extract ⟨A⟩ by three methods: Gaussian fit (±0.5), truncated mean 90%, truncated mean 95%
4. Compute R = (1 + ⟨A⟩) / (1 − ⟨A⟩) for DATA and MC; form R_data/R_MC
5. Fit R_data/R_MC vs alpha threshold (up to 0.30) with a linear function; intercept = alpha→0 correction

**Full-η pass** — identical loop on the unfolded (signed η) THnSparse, producing a second set of corrections with 36 bins covering −5.191 to 5.191.

Output per cone: 15 intercept TH1Ds for |η| (`{cone}_intercept_{method}{ptavg}`) and 15 for full η (`{cone}_intercept_{method}{ptavg}_fulleta`), plus TGraphErrors + fits and QA A-distributions.

**Three mean estimators** provide a cross-check on systematic sensitivity:
- *Gauss fit*: fits a Gaussian in A ∈ [−0.5, 0.5]; robust against non-Gaussian tails
- *Trunc. 90%*: trims 5% from each tail before computing the mean; model-independent
- *Trunc. 95%*: trims 2.5% from each tail; intermediate

### Step 2 QA — Comparison Plots

`plotResiduals` reads the residuals ROOT file and writes PNG files to a timestamped directory. Two plot types, each with a main panel and a ratio/difference panel below:

**η symmetry check** (`etasym_{cone}_{method}_{ptavg}.png`):
- Top: full-η corrections vs. |η| corrections reflected about η = 0
- Bottom: (full η) / (reflected |η|) — a flat line at 1.0 means perfect forward-backward symmetry

**Method comparison** (`methods_{cone}_{ptavg}_{abseta|fulleta}.png`):
- Top: Gauss, trunc90, trunc95 overlaid
- Bottom: trunc90/Gauss and trunc95/Gauss — flat at 1.0 means the correction is method-stable

### Step 3 — Text File Output *(not started)*

The pT dependence in each |η| bin is fit with `1 / (p0 + p1·log10(0.01·pT) + p2/(pT/10))`. The text file has one line per eta bin (negative eta mirrored from positive), in the standard JEC format understood by `JetCorrector`.

---

## Repository Layout

```
include/
  Dijet.h            FindLeadingJets (O(n) single-pass), MakeDijet — all index-based
  DijetHistograms.h  ConeHistograms struct — one per cone size, owns the 4D THnSparse
  Binning.h          BinningConfig — all axis definitions, pT/alpha extraction slices
  EventStructs.h     EventStruct (vz, weight, run/lumi), FiltersStruct (ppvF)
  JetStruct.h        JetStruct<MAXNREF> — reco/ref/gen arrays, PF info, BranchMap()
  BranchMapping.h    SetBranches() — one-call TTree branch binding from a struct map
  JetCorrector.h     Yi Chen v3.0 — chains text-file JEC levels sequentially
  JetSelection.h     JetSelect — tight pp jet ID (eta-dependent) + veto map from ROOT file
  JSON_handler.h     golden JSON checking (requires nlohmann-json)
  ResidualsExtractor.h  Step 2 entry point declaration
  RunAsymmetry.h        Step 1 entry point declaration
  Utilities.h        MakeTHnSparse, FoldEtaAxis, canvas/style helpers
  Colors.h           Hiroshige + Klimt color palettes

cfg/
  2024ppRef.h        JEC file paths, trigger branch names + thresholds, TTree paths, cone labels

macros/
  runAsymmetry.C     Step 1 entry point — unified DATA/MC, configured by --mc/--zero-bias flags
  runResiduals.C     Step 2 entry point — |eta| and full-eta extraction
  plotResiduals.C    Step 2 QA — eta symmetry check and method comparison plots

src/
  RunAsymmetry.cxx       Step 1 implementation
  ResidualsExtractor.cxx Step 2 implementation — ExtractAndFit helper called twice per cone

tests/
  TestDijet.cxx          37 tests — FindLeadingJets and MakeDijet
  TestFoldEtaAxis.cxx     7 tests — FoldEtaAxis

readme/              versioned READMEs (this is v5)
```

---

## Step 2 Output Structure

For each cone (e.g. `ak4PF`):

| Object | Description |
|--------|-------------|
| `ak4PF_asym_data_abseta` | Folded \|η\| THnSparse from DATA |
| `ak4PF_asym_mc_abseta` | Folded \|η\| THnSparse from MC |
| `ak4PF_QA_data/` | 1D A distributions, \|η\| pass (one per ptavg × alpha × eta cell) |
| `ak4PF_QA_mc/` | same for MC |
| `ak4PF_QA_data_fulleta/` | 1D A distributions, full-η pass |
| `ak4PF_QA_mc_fulleta/` | same for MC |
| `ak4PF_graphs/` | TGraphErrors of R_data/R_MC vs alpha threshold, with linear fit |
| `ak4PF_intercept_{method}{ptavg}` | TH1D of alpha→0 corrections vs \|η\|, 18 bins |
| `ak4PF_intercept_{method}{ptavg}_fulleta` | same vs signed η, 36 bins |

`{method}` ∈ {`gauss`, `trunc90`, `trunc95`}  
`{ptavg}` ∈ {`_ptavg_40_90`, `_ptavg_90_120`, `_ptavg_120_190`, `_ptavg_190_260`, `_ptavg_260_1000`}

---

## Key Design Choices

### One 4D THnSparse, one fill per event

The original code used `TH2D[9][5]` — 9 alpha thresholds × 5 pT_avg slices = 45 histograms (plus 45 more for full η). A 2-jet event was filled into all 9 alpha histograms; a 3-jet event was filled into as many alpha bins as it passed. This inflated entry counts and made histograms non-independent.

Here, one 4D THnSparse with axes `(η_probe, pT_avg, alpha, A)` is filled exactly once per event at the event's actual alpha value (0 for 2-jet events). Cumulative alpha slices — "all events with alpha < 0.10", etc. — are recovered at analysis time with:

```cpp
h->GetAxis(kAlphaAxis)->SetRangeUser(0, threshold);
TH1D* proj = (TH1D*)h->Projection(kAAxis);
```

A 2-jet event with alpha = 0 falls below every positive threshold automatically, reproducing the original cumulative behavior without any special-casing.

### Full η stored; |η| folded at extraction

The THnSparse stores full η_probe across all 36 CMS JEC bins. At extraction, `FoldEtaAxis()` in `Utilities.h` produces a new THnSparse with |η| as the axis. Both the folded and unfolded versions are used in Step 2: the |η| intercepts are the primary result; the full-η intercepts provide the forward-backward symmetry check via `plotResiduals`.

### `ExtractAndFit` — one helper for both eta types

Rather than duplicating the extraction loop, `ResidualsExtractor.cxx` defines a static `ExtractAndFit(hData, hMC, ..., etaEdges, nameSuffix)` function. `runResiduals` calls it twice: once with the folded |η| sparses and `kAbsEtaEdges`, and once with the raw full-η sparses and `kEtaEdges`. The `nameSuffix` (`""` or `"_fulleta"`) distinguishes all output object names.

### Index-based dijet logic

`FindLeadingJets(corrPt, nref)` returns a `SortedJets{lead, sublead, third}` — indices into the corrected-pT array, not copies. `MakeDijet` takes those indices plus the raw kinematic arrays. Jet corrections change the pT ordering, so indices computed on corrected pT are the only sensible reference. Jet ID (which needs PF fractions and multiplicities) is applied in the event loop between sorting and dijet construction, without touching `Dijet.h`.

### `event % 2` for barrel-barrel tag assignment

CMS event numbers are sequential integers, so `eventNumber % 2` gives an exact, reproducible 50/50 tag/probe split per event. The original used `TRandom2::Integer(100) % 2` with seed 1 — this is sequence-dependent across job restarts and not guaranteed to be 50/50 for any finite sample.

### Jet ID on the third jet

The original used the third jet's corrected pT unconditionally. A fake jet that passes pT sorting but fails quality criteria would artificially inflate alpha, making the event appear radiation-rich and suppressing it from tight alpha slices. Here, if the third jet fails jet ID, `hasThird = false` and `alpha = 0` — the event is treated as a clean 2-jet event for the purpose of alpha computation.

---

## Binning

| Quantity | Configuration |
|----------|--------------|
| pT_avg axis | 990 bins × 1 GeV, 10–1000 GeV (coarse slices applied at extraction) |
| pT_avg slices | 40–90, 90–120, 120–190, 190–260, 260–1000 GeV |
| η_probe | 36 CMS JEC standard bins, −5.191 to 5.191 (variable width) |
| \|η_probe\| | 18 bins, 0 to 5.191 (folded at extraction) |
| Alpha axis | 50 bins, 0–0.5 |
| Alpha thresholds | 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45 |
| Asymmetry A | 100 bins, −1 to 1 |

Alpha fits for the extrapolation use thresholds up to 0.30 (first 6 of 9); above this, trigger bias distorts the distribution.

---

## Remaining Work

| Task | Status |
|------|--------|
| Stage AK2/3/5/6 L2Relative JEC files into `data/jec/` | Config references them, files not staged |
| Confirm `kTrigTreePath` and `kHLTJ80Branch` version suffix | TODO: check from actual HiForest file |
| Validate Step 1 on test MC sample | Not done |
| Validate Step 2 on Step 1 output | Not done |
| Step 3 — `src/TextFileWriter.cxx` + `macros/runTextFile.C` | Not started |

---

## Adapting for a New Run Period

Add a new config header under `cfg/` specifying:
- JEC text file paths (L2Relative + L2Residual for the new period)
- Veto map ROOT file path
- Golden JSON path
- Trigger TTree path and HLT branch names with efficiency thresholds
- Cone size labels

Pass the new config at build time (or swap the `#include` in the macro). All analysis logic, histogram structure, and extraction code stays unchanged.

---

## Contact

Nicholas Barnett — thenicholasbarnett@gmail.com
