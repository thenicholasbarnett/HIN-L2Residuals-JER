# Style Guide

This is the reference for the cleanup pass described in `cleanup_plan.md`. It exists so cleanup has one target instead of drifting file to file. Style reference repo: `CMS_2026PbPb`.

## Hard rule: never touch vendored code

**Nothing under `external/` gets touched, ever — not formatting, not comments, not casts, nothing.** This includes all vendored JetMET/JetMETAnalysis and CMSSW `CondFormats/JetMETObjects` source (`external/jetmet/`, `external/jetmet_jer/`, `external/nlohmann/`, `external/toml.hpp`). This is not a style preference, it's a correctness boundary — vendored code stays byte-identical to upstream so it can be diffed/updated against its real source later. If a vendored file genuinely needs a local patch (as `JetResolutionObject.cc` already does for `STANDALONE` mode), that's a deliberate, clearly-commented deviation, not a cleanup pass.

## Control flow and spacing

```cpp
if ( condition ) { continue; }
if ( !ptr ) { return; }
for ( size_t i = 0; i < n; i++ ) { ... }
functionCall( arg1, arg2 );
```

Every `if`/`for`, including one-line bodies and guard clauses, uses `{ }` — no bare `if ( x ) return;`. Spacing is for readability inside an expression, not for lining up unrelated lines. No column alignment on declarations, assignments, or struct members — one line each, no padding to a common column.

```cpp
// prefer
std::string input = cl.getValue<std::string>( "input" );
std::string output = cl.getValue<std::string>( "output" );

// avoid
std::string input   = cl.getValue<std::string>( "input" );
std::string output  = cl.getValue<std::string>( "output" );
```

Casts stay C-style (`( Type )value`), matching ROOT/CMS convention and the rest of this codebase — not `static_cast`.

Variable names like `fi`, `fo`, `dcs`, `js` stay as-is. That's the author's actual working style, not a mistake to fix.

## Comments

Default is no comment. Casual, first-person voice is fine ("we skip this because...") but the target is brevity — most functions need zero comments, not a friendlier version of the same paragraph. One line is normal; more than three lines means it probably belongs in docs instead.

Keep a comment only for:
- A real ROOT ownership/lifetime trap.
- A CMS/JEC/JER text format detail that isn't obvious from the code.
- A build or environment quirk that can't be named away.
- Short section labels in long event loops: `// event loop`, `// trigger`, `// output`.

Avoid "must," "never," "always," "not stale," "physics-affecting choice" — legalistic hedging, not information.

**No physics explanations in source comments.** Don't write "we take the limit as alpha goes to zero because..." in code — that reasoning lives in docs. A short label naming *what* something is (e.g. `// alpha extrapolation`, `// third-jet pT`) is fine; a paragraph justifying *why* is not.

Debugging/incident history (annobin, RPATH, SCRAM path quirks, CMSSW compatibility fights) moves to a single docs file, not scattered inline comments or repeated across files.

## Before / after (real examples from this repo)

**Alignment** (`macros/runAsymmetry.C`):
```cpp
// before
std::string input   = cl.getValue<std::string>( "input" );
std::string output  = cl.getValue<std::string>( "output" );
std::string mode    = cl.getValue<std::string>( "mode" );

// after
std::string input = cl.getValue<std::string>( "input" );
std::string output = cl.getValue<std::string>( "output" );
std::string mode = cl.getValue<std::string>( "mode" );
```

**Guard clauses** (`src/CalibrationExtractor.cxx`):
```cpp
// before
if( !CanFit( h, controls.minEntries ) ) return r;

// after
if ( !CanFit( h, controls.minEntries ) ) { return r; }
```

**Long function-header comment** (`src/CalibrationExtractor.cxx`, `ExtractAndFit`) — the naming schema and JER SF backstory move to docs, a one-line pointer stays:
```cpp
// before (36 lines: naming schema + JER SF derivation history + Nicky's-call narration)

// after
// Runs the extraction loop for one eta mode. Output naming and the JER SF
// derivation are documented in docs/calibration-extraction.md.
```

**CMake defensive comments** (`CMakeLists.txt`):
```cpp
// before
# Must be before project(): CMSSW GCC has a broken annobin plugin and
# system ROOT on lxplus lacks libCling.so — both fixed with CMSSW GCC + CMSSW ROOT.
if(DEFINED ENV{CMSSW_BASE})
    # CMSSW GCC has broken annobin plugin specs (the .so is missing from cvmfs).
    # Fix: dump its specs, strip annobin entries, compile with -specs=<cleaned>.

// after
# CMSSW GCC's annobin plugin is broken on cvmfs; strip it from the spec file.
if(DEFINED ENV{CMSSW_BASE})
```
