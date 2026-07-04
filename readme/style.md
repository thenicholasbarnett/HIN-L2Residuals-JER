# Style Guide

Reference for the cleanup pass in `cleanup_plan.md`, and for any new code written in this repo going forward. One target instead of drifting file to file. Style/voice reference repo: `CMS_2026PbPb`.

## Hard rule: never touch vendored code

**Nothing under `external/` gets touched, ever — not formatting, not comments, not casts, nothing.** This covers all vendored JetMET/JetMETAnalysis and CMSSW `CondFormats/JetMETObjects` source (`external/jetmet/`, `external/jetmet_jer/`, `external/nlohmann/`, `external/toml.hpp`). Not a style preference — a correctness boundary. Vendored code stays byte-identical to upstream so it can be diffed/updated against its real source later. A deliberate local patch (e.g. `JetResolutionObject.cc`'s `STANDALONE` fix) is a rare, clearly-flagged exception, never part of a sweep.

## No AI attribution, anywhere

No AI/assistant-tool name, no "generated with," no watermark of any kind — not in code, comments, commit messages, or docs. This repo's code should read as entirely human-authored. (Same rule as every repo Nicky owns — see the global CLAUDE.md / `~/.dotfiles/CODE_STYLE.md`.)

## Control flow and spacing

```cpp
if ( condition ) { continue; }
if ( !ptr ) { return; }
for ( size_t i = 0; i < n; i++ ) { ... }
while ( cond ) { ... }
functionCall( arg1, arg2 );
```

- Space between the control keyword and `(`: `if (`, `for (`, `while (` — not `if(`.
- Space just inside both the condition parens and the body braces: `if ( x ) { y; }`.
- No space between a function/method name and its `(`: `functionCall( arg )`, not `functionCall ( arg )`.
- **Every `if`/`for`/`while`, including one-line bodies and guard clauses, uses `{ }`.** No bare `if ( x ) return;`.
- This is the target going forward — most of the existing codebase currently uses `if( x ){` with no space before the paren. Bringing existing code in line with the new spacing is Phase 3 (Mechanical Spacing Cleanup) of `cleanup_plan.md`, not something to chase file-by-file outside that phase.

No column alignment on declarations, assignments, or struct members — one line each, no padding to a shared column. Alignment is fine only for genuine tables (bin-edge arrays, lookup tables) where the columns really are data.

```cpp
// prefer
std::string input = cl.getValue<std::string>( "input" );
std::string output = cl.getValue<std::string>( "output" );

// avoid
std::string input   = cl.getValue<std::string>( "input" );
std::string output  = cl.getValue<std::string>( "output" );
```

Casts stay C-style — `( Type )value`, not `static_cast<Type>( value )`. Matches ROOT/CMS convention and the rest of this codebase.

Variable names like `fi`, `fo`, `dcs`, `js` stay as-is. That's the author's real working style, not a mistake to fix.

## Structure

Include order, seen consistently across this repo — keep it:
1. The file's own header (`.cxx` files only).
2. ROOT headers (`TFile.h`, `TString.h`, ...).
3. This project's own headers, roughly dependency order.
4. Standard library headers last.

Short labeled dividers for grouping related constants/axes/structs, e.g.:
```cpp
// axes
static constexpr int kEtaAxis = 0;
static constexpr int kPtAvgAxis = 1;

// methods
static constexpr int kNMethods = 4;
```
or the heavier `// ---- result structs ----` form for a bigger section break. Either is fine; don't invent a third style.

Guard clauses go first in a function, one condition per line, each returning/continuing immediately — don't nest error handling.

## Comments

Default is no comment. Casual, first-person voice is fine ("we skip this because...") but the target is brevity — most functions need zero comments, not a friendlier version of the same paragraph. Terse fragments over complete sentences. One line is normal; past three lines it probably belongs in docs instead.

Keep a comment only for:
- A real ROOT ownership/lifetime trap.
- A CMS/JEC/JER text format detail that isn't obvious from the code.
- A build or environment quirk that can't be named away.
- Short section labels in long event loops or long files: `// event loop`, `// trigger`, `// output`.

Avoid "must," "never," "always," "not stale," "physics-affecting choice" — legalistic hedging, not information.

**No physics explanations in source comments.** Don't write "we take the limit as alpha goes to zero because..." in code — that reasoning lives in docs. A short label naming *what* something is (`// alpha extrapolation`, `// third-jet pT`) is fine; a paragraph justifying *why* is not.

Debugging/incident history (annobin, RPATH, SCRAM path quirks, CMSSW compatibility fights) moves to one docs file, not scattered inline comments or repeated across files.

## Macro/CLI headers

Keep the compiled usage string and one compact interpreted-ROOT example. Don't repeat the full `CommandLine` parser philosophy in every macro file — that lives once in the README.

## Before / after (real examples from this repo)

**Control-flow spacing + braces** (`src/RunAsymmetry.cxx`):
```cpp
// before
if( modeFlag == "mc" ) mode = RunMode::MC;
else if( modeFlag == "non-triggered" ) mode = RunMode::NonTriggered;

// after
if ( modeFlag == "mc" ) { mode = RunMode::MC; }
else if ( modeFlag == "non-triggered" ) { mode = RunMode::NonTriggered; }
```

**Alignment** (`macros/runAsymmetry.C`):
```cpp
// before
std::string input   = cl.getValue<std::string>( "input" );
std::string output  = cl.getValue<std::string>( "output" );

// after
std::string input = cl.getValue<std::string>( "input" );
std::string output = cl.getValue<std::string>( "output" );
```

**Placeholder comments given real meaning** (`src/RunAsymmetry.cxx`):
```cpp
// before
// ?
size_t trigConeIdx = nCones;

// after
// find the trigger cone
size_t trigConeIdx = nCones;
```

**Legalistic comment tightened** (`src/RunAsymmetry.cxx`):
```cpp
// before
// JetCorrector — chain per-cone L2Residual files onto the base JEC, but
// only for data (Triggered/NonTriggered); MC must never see residual
// corrections, so this check is the only thing that decides that, not
// what's in the TOML.

// after
// chain L2Residual files onto the base JEC for data only -- mode decides
// this, never the TOML
```

**Long function-header comment** (`src/CalibrationExtractor.cxx`, `ExtractAndFit`) — naming schema and JER SF backstory move to docs, one-line pointer stays:
```cpp
// after
// Runs the extraction loop for one eta mode. Output naming and the JER SF
// derivation are documented in docs/calibration-extraction.md.
```

**CMake defensive comments** (`CMakeLists.txt`):
```cpp
// before
# Must be before project(): CMSSW GCC has a broken annobin plugin and
# system ROOT on lxplus lacks libCling.so — both fixed with CMSSW GCC + CMSSW ROOT.

// after
# CMSSW GCC's annobin plugin is broken on cvmfs; strip it from the spec file.
```
