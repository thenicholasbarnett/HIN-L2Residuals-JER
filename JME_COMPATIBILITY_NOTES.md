# JME Compatibility Notes

This repo should feel familiar to the JME POG wherever that helps review,
maintenance, or shared calibration work. Prefer JetMETAnalysis conventions and
well-tested utilities when they fit this analysis cleanly, especially for
command-line syntax, correction naming, text payload formats, and plotting
style. Compatibility is the point; blind copying is not. When JetMET code is
tightly coupled to JRA ntuples, CMSSW packages, or a different histogram layout,
use it as a reference unless importing it makes the local workflow simpler and
more robust.

## Current Decision: Command-Line Interface

**Superseded (2026-07-02): fully vendored, `include/CliTokens.h` deleted.**
The compiled binaries now use JetMET's real `CommandLine` class directly —
`JetUtilities/interface/CommandLine.h`/`.cc` copied verbatim into
`external/jetmet/CommandLine.h`/`.cc` (attribution banner in both files), no
compatibility shim on top. This was Nicky's explicit call: "if we can delete
CliTokens at the end of this I'll be happy" — reliance on JetMET's own
architecture was the goal, not a lookalike reimplementation.

```bash
./build/bin/runResponse \
  -input mc.root \
  -output response.root \
  -config cfg/2024ppRef.toml
```

Also accepts an initial argument file with JetMET's own `key = value` syntax:

```text
input = mc.root
output = response.root
config = cfg/2024ppRef.toml
flags = finals etasym methods
```

Two real, deliberate behavior changes from the old `CliTokens.h` came with
this: keys are now **case-sensitive** (`CommandLine`'s own convention — no
more `-INPUT`/`-Input` tolerance), and the legacy `KEY=value` shell-token form
is gone entirely, including from `condor/runtime_wrapper.sh`'s
script-generated invocation line, which now emits `-input`/`-output`/
`-mode`/`-config`. `condor/make_condor.sh` is unaffected — it's a bash script
and never called into `CliTokens.h` or `CommandLine.h`, both C++.

## JetMET Code Survey

Strong candidates to borrow or adapt:

| JetMET code | Why it helps | Suggested use here |
| --- | --- | --- |
| `JetUtilities/interface/CommandLine.h` and `src/CommandLine.cc` | Standard JME `params.config -key value` parser. Supports config includes, vectors, defaults, required values, and unused-option checks. | Consider vendoring if exact JME behavior matters more than keeping a tiny local parser. |
| `JetUtilities/interface/JetInfo.hh` and `src/JetInfo.cc` | Encodes JME jet algorithm labels, correction level names, legend labels, detector regions, and string helpers. | Adapt a self-contained subset for official correction level naming and plot labels. Avoid importing package-level globals unless needed. |
| `JetAnalyzers/bin/jet_apply_jec_x.cc` | Reference implementation for applying JEC with CMSSW `FactorizedJetCorrector` and `JetCorrectorParameters`. | Use as a model for an optional CMSSW-backed corrector or validation binary. It is tied to JRA inputs, so do not copy directly into the HiForest path. |
| `JetUtilities/interface/HistogramUtilities.hh` | Shared response metric helpers: histogram mean/RMS, fit mean/RMS/sigma, median, and errors. | Useful reference for unifying asymmetry/response fit metric handling and naming. |
| `JetAnalyzers/bin/JERWriter.h` | Writes official-style JER text payloads from fitted functions and bin definitions. | Study for future JER payload export. Current L2Residual output is JEC-focused, so this is not immediate. |
| `JetUtilities/interface/RootStyle.h` and `Style.h` | CMS/JME plotting style helpers: TDR style, legends, CMS labels, marker/color conventions. | Good candidate for plot style alignment, either by copying a small subset or using it as a style reference. |
| `JetUtilities/interface/Variables.hh` and `src/Variables.cc` | Common variable names and axis titles for `refpt`, `jtpt`, `refeta`, `rho`, etc. | Adapt small mappings for plot labeling consistency. |
| `JetUtilities/interface/ProgressBar.hh` | Tiny terminal progress helpers. | Optional convenience only; easy to reimplement locally if wanted. |

Useful primarily as reference:

| JetMET code | Why it is useful | Caveat |
| --- | --- | --- |
| `JetAnalyzers/bin/jet_response_fitter_x.cc` | Shows JME response fitting flow, Gaussian/DSCB options, histogram naming, and fit failure handling. | Assumes JRA response histogram layout rather than this repo's sparse/asymmetry workflow. |
| `JetAnalyzers/bin/jet_response_and_resolution_x.cc` | Documents JME response/resolution production conventions. | Mostly coupled to JetMET input/output naming. |
| `JetUtilities/interface/L2Creator.hh` and `src/L2Creator.cc` | Full JME L2Relative correction creation machinery, including L3 response inputs, graph fitting, splines, and text output. | Heavy CMSSW/JetMET/GSL/ObjectLoader dependency stack. Best used as validation/reference before importing pieces. |
| `JetUtilities/interface/ClosureMaker.hh`, `RatioMaker.hh`, and related sources | Closure and ratio plot patterns that may match future QA needs. | Depend on JetMET histogram organization and plotting assumptions. |

Low-priority for this repo:

| JetMET code | Reason |
| --- | --- |
| JRA event/object loaders and pileup ntuple helpers | They target JRA/PU trees, while this repo reads HiForest and analysis-specific ROOT products. |
| Synch, mass, and weighted-spectrum binaries | Useful within JetMETAnalysis, but not directly part of the L2Residual dijet balance path. |

## Working Rule

When the repo needs functionality that JetMETAnalysis already implements, first
ask whether the JetMET version improves compatibility with JME review and shared
operations. If yes, prefer adopting it or matching its behavior. If the JetMET
implementation brings in broad dependencies or incompatible data assumptions,
keep the local implementation but document the JetMET behavior it is matching.
