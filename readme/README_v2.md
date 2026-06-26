# L2Residuals-2024ppRef

An object-oriented C++/ROOT framework for L2 Residual jet energy corrections (JEC) using dijet pT balance, developed for the CMS 2024 pp reference run (√s = 13.6 TeV).

---

## Physics Background

L2 Residual corrections remove the residual eta-dependence of the jet energy scale in data after MC-truth (L2Relative) corrections have been applied. The measurement uses dijet pT balance: in a back-to-back dijet event, a well-calibrated **tag** jet in the barrel (|eta| < 1.3) constrains the response of a **probe** jet anywhere in the detector.

The dijet asymmetry is defined as:

```
A = (pT_probe - pT_tag) / (pT_probe + pT_tag)
```

The response ratio R = (1 + <A>) / (1 - <A>) is computed for data and MC separately. Their ratio — the double ratio R_MC / R_data — gives the multiplicative correction to apply to data jets at each eta and pT. A linear extrapolation of the double ratio as a function of alpha (pT_third / pT_avg) to alpha = 0 removes the bias from soft radiation.

---

## Analysis Steps

### Step 0 — Trigger Turn-On

Determine the fully efficient pT threshold for each jet trigger before measuring asymmetries. For the 2024 pp reference run:

| Trigger | Efficient above |
|---------|----------------|
| HLT_AK4PFJet40 | ~45 GeV (leading jet pT) |
| HLT_AK4PFJet80 | ~95 GeV (leading jet pT) |

### Step 1 — Asymmetry Generation

A single unified generator runs over both data and MC HiForest files, configured by a `bool isData` flag. Output is a `THnSparseD` with axes:

- axis 0: |eta_probe|
- axis 1: asymmetry A
- axis 2: pT_avg
- axis 3: alpha

### Step 2 — Residuals Extraction

Projects the sparse histogram into 1D asymmetry distributions per (pT_avg, |eta|, alpha) bin, extracts <A>, forms the double ratio, and extrapolates to alpha = 0 with a linear fit.

### Step 3 — Text File Generation

Fits the pT dependence of the correction in each |eta| bin with a 3-parameter functional form and writes a JEC-standard text file compatible with `JetCorrector`/`SingleJetCorrector`.

---

## Binning (2024 pp Reference)

| Quantity | Values |
|----------|--------|
| pT_avg | 40–90, 90–120, 120–190, 190–260, 260–1000 GeV |
| \|eta\| | 18 bins, 0 to 5.191 (CMS JEC standard) |
| Alpha | 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45 |
| Asymmetry | 100 bins, −1 to 1 |

---

## Requirements

- ROOT 6
- CMake ≥ 3.10
- C++17

## Build

```bash
mkdir build && cd build
cmake ..
make
```

## Running

```bash
# Step 1 — asymmetry generation
./build/runAsymmetry input.root output.root --data    # data
./build/runAsymmetry input.root output.root --mc      # MC

# Step 2 — residuals extraction
./build/runResiduals data_asymmetry.root mc_asymmetry.root residuals.root

# Step 3 — text file
./build/runTextFile residuals.root L2Residuals_2024ppRef.txt
```

For Condor batch jobs, ship the compiled executables and configure `condor/runtime_wrapper.sh` with the appropriate CMSSW/ROOT environment setup.

---

## Repository Structure

```
L2Residuals-2024ppRef/
├── include/
│   ├── BinningConfig.h       # pT_avg, eta, alpha, asymmetry bin definitions
│   ├── BranchMapping.h       # automated TTree branch binding
│   ├── Colors.h              # plot color palettes (Hiroshige, Klimt)
│   ├── Dijet.h               # tag/probe assignment, asymmetry, alpha
│   ├── EventStructs.h        # EventStruct, FiltersStruct with BranchMap()
│   ├── JetCorrector.h        # Yi Chen v3.0 — chains text-file JEC levels
│   ├── JetSelection.h        # pp jet ID cuts + veto map
│   ├── JetStruct.h           # reco/ref/gen jet arrays with BranchMap()
│   ├── ProgressBar.h         # terminal progress bar
│   └── Utilities.h           # THnSparse helpers, canvas/legend/style
├── src/
│   ├── AsymmetryAnalysis.cxx # unified DATA/MC event loop
│   ├── ResidualsExtractor.cxx
│   └── TextFileWriter.cxx
├── macros/
│   ├── runAsymmetry.C
│   ├── runResiduals.C
│   └── runTextFile.C
├── configs/
│   └── 2024ppRef.h           # paths, trigger thresholds, correction files
└── CMakeLists.txt
```

---

## Contact

Nicholas Barnett — thenicholasbarnett@gmail.com
