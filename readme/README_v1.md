# L2Residuals-2024ppRef

Object-oriented C++/ROOT framework for computing L2 Residual jet energy corrections from dijet pT balance, written for the 2024 pp reference run at CMS.

## Overview

L2 Residual corrections account for residual eta-dependent jet energy scale differences between data and MC after L2Relative (MC-truth) corrections have been applied. The method uses back-to-back dijet events: a tag jet in the barrel constrains the response of a probe jet at arbitrary eta.

The analysis runs in three steps:

1. **Asymmetry generation** — runs over data and MC forests, fills high-dimensional asymmetry histograms binned in pT_avg, eta, and alpha
2. **Residuals extraction** — slices histograms, extracts mean asymmetry, forms the data/MC double ratio, extrapolates alpha → 0
3. **Text file generation** — fits the pT dependence and writes a JEC-format correction text file

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

Executables are placed in `build/`.

## Running

```bash
./build/runAsymmetry input.root output.root --data
./build/runAsymmetry input.root output.root --mc

./build/runResiduals data_asymmetry.root mc_asymmetry.root residuals.root
./build/runTextFile residuals.root L2Residuals_2024ppRef.txt
```

For batch processing, ship the compiled executables to Condor via the runtime wrapper in `condor/`.

## Repository Structure

```
include/    headers — structs, branch mapping, jet corrections, utilities
src/        implementations — compiled into libl2residuals.so
macros/     entry points — thin executables linked against the library
configs/    run-specific configuration (paths, trigger thresholds, correction files)
```

## Configuration

Run-specific parameters (correction text file paths, trigger pT thresholds, jet algorithm) are set in `configs/2024ppRef.h`. To adapt the framework for a different run period, add a new config header and pass it to the analysis class.

## Contact

Nicholas Barnett — thenicholasbarnett@gmail.com
