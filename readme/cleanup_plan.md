# L2Residuals Cleanup Plan

**Hard rule, before anything else: `external/` (vendored JetMET/JetMETAnalysis/CMSSW/nlohmann/toml++ source, plus `JetCorrector.h`/`JSON_handler.h` -- code written by other people, relocated into `external/` on 2026-07-04 for exactly this reason) is never touched by this cleanup — not formatting, not comments, nothing. See `style.md`.**

This plan is for making the repository feel plainly authored and maintainable without erasing the real analysis history or replacing Nicky's code style with a generic formatter style.

Style decisions are resolved — see `style.md` and the "Resolved" answers under Questions For Nicky below.

The goal is not to make the code look "machine perfect." The goal is to make it look like a person owns it: direct comments, consistent local spacing, build logic that reads as designed, and docs that carry long explanations instead of burying them in executable files.

## Guiding Style

**Mechanical style (spacing, braces, alignment) is entirely `.clang-format`'s job — run it, don't hand-author it.** (2026-07-05: reconciled with Nicky — an earlier version of this doc hand-specified a spaced-paren convention, `if ( condition )`, as the target and assigned bringing the repo into line with it to Phase 3. That was never actually encoded in `.clang-format`, which is `BasedOnStyle: LLVM` + `ReflowComments: false` + `SortIncludes: false` + `InsertBraces: true` — real output is unspaced parens and clang-format-inserted braces, see below. `.clang-format` is the actual source of truth now; this doc's job is comments and things a formatter can't decide, not re-describing what `clang-format -style=file` already produces.)

```cpp
if (condition) {
  continue;
}
if (!ptr) {
  return;
}
for (size_t i = 0; i < n; i++) {
  ...
}
functionCall(arg1, arg2);
```

`InsertBraces: true` means every `if`/`for`/`while`, including one-line bodies and guard clauses, gets braces on format — no bare `if (x) return;` survives a pass. No column alignment on declarations/assignments — clang-format's LLVM base already defaults to no consecutive-alignment, so this falls out for free:

```cpp
std::string input = cl.getValue<std::string>("input");
std::string output = cl.getValue<std::string>("output");
bool useNorm = cl.getValue<bool>("norm", true);
```

Comments should mostly answer "why." Short section labels are fine when they help navigation. Long failure histories should move to docs unless they are needed right beside the code.

Comment brevity is a hard target. After cleanup, the default should be no comment. Add one only when it helps a reader avoid a real misunderstanding. A comment should usually be one short line. Multi-line comments are allowed for physics choices, ROOT ownership/lifetime traps, CMS text-file formats, and build quirks that cannot be made clear by naming. `ReflowComments: false` means clang-format never touches comment wording or wraps them — comment content is entirely hand-authored, always.

Style reference: `CMS_2026PbPb` is closer to the desired voice than the current L2Residuals build/Condor layer. It uses short section markers like `// looping over events`, direct comments like `// event objects`, and compact guards. The cleanup should borrow that restraint for comments; control-flow spacing is `.clang-format`'s call, not a voice choice.

## What Feels Off Today

- The top-level build files carry a lot of environment history inline.
- `condor/make_condor.sh` reads like a manual and an implementation at the same time.
- The same policy is repeated in README, macros, C++, tests, and Condor scripts.
- Many comments are correct but too defensive in tone: "must never," "no implicit default," "physics-affecting choice," "not stale," etc.
- Visual alignment spacing makes code look generated or over-polished.
- Some source files have human roughness, like `// ?`, `// inpt`, and compact one-line guards. Those are not the problem. They should either be clarified or left alone if they express the way the author actually works.
- Compared with `CMS_2026PbPb`, the current repo has too many "because this could be misunderstood" comments. The desired endpoint is closer to short navigation comments plus a few carefully chosen explanations.

## Comment Brevity Rules

Use these rules during every cleanup phase.

- [ ] Prefer deleting a comment over rewriting it.
- [ ] Keep a comment only if code, naming, or nearby docs cannot carry the meaning.
- [ ] Keep section comments short: `// event loop`, `// trigger`, `// output`.
- [ ] Avoid comments that restate the next line.
- [ ] Avoid long comments in headers unless they are API-critical.
- [ ] Move debugging history into docs.
- [ ] Move command examples into README or topic docs.
- [ ] Avoid "must/never/always" unless it is a real correctness rule.
- [ ] Avoid explaining the same policy in multiple files.
- [ ] If a comment needs more than three lines, first ask whether it belongs in docs.
- [ ] If a comment exists only to make generated-looking code feel safe, delete it.
- [ ] Preserve short human notes when they help orientation.

Target density:

- [ ] Event-loop source files may have short landmarks.
- [ ] Small headers should have almost no prose.
- [ ] Build files should have short headings and links/paths to deeper notes.
- [ ] Condor scripts should have usage plus sparse comments near non-obvious shell behavior.
- [ ] Tests should use assertion messages more than comments.

## Phase 0: Preserve The Current State

Do this before changing style or structure.

- [ ] Check the current branch and working tree.
- [ ] Note existing untracked files, especially `README_draft.md` and `readme/`.
- [ ] Decide whether cleanup happens on a new branch.
- [ ] Build once from the current state if the environment supports it.
- [ ] Run available tests once from the current state if the environment supports it.
- [ ] Save current command outputs in notes if build/test behavior is fragile.
- [ ] Do not combine behavior changes with style cleanup in the first pass.

Best time: first.

Reason: if a later cleanup breaks something, there needs to be a known baseline.

## Phase 1: Lock In The Human Style

This phase should be small and explicit. It prevents the cleanup from drifting into generic clang-format territory.

- [ ] Add a short style note, probably `readme/style.md` or `docs/style.md`.
- [ ] Document the intentional control-flow style: `if ( condition ) { continue; }`.
- [ ] Document spacing around function calls: `f( arg )`, ROOT-style casts if kept, and compact guards.
- [ ] Document that column alignment is not wanted for declarations or assignments.
- [ ] Document when one-line `if` statements are acceptable.
- [ ] Document that comments should be shorter in code and longer in docs.
- [ ] Document that comments should be absent by default.
- [ ] Document the `CMS_2026PbPb` style reference.
- [ ] Document that ROOT/CMS names should not be renamed just to look modern.
- [ ] Include 4-6 before/after examples from the repo.

Best time: before any broad edit.

Reason: this gives a target style and avoids arguing with future automated tooling.

## Phase 2: Calibration Edits With Nicky

Before sweeping the whole repo, edit one or two files by hand to learn the preferred voice.

Suggested files:

- [ ] `src/RunAsymmetry.cxx`
- [ ] `include/Dijet.h`
- [ ] `include/JetSelection.h`
- [ ] Optional: one macro such as `macros/runTextFile.C`

Small tasks:

- [ ] Replace `// ?` in `src/RunAsymmetry.cxx` with a real explanation.
- [ ] Replace `// inpt` with either `// input file` or no comment.
- [ ] Decide whether comments like `// structures`, `// sort`, `// trig` should stay as navigation marks.
- [ ] Compare the result against `CMS_2026PbPb` for comment density.
- [ ] Delete at least one comment that is merely decorative.
- [ ] Remove column alignment in a few variable declarations.
- [ ] Keep one-line guard clauses where they make the event loop easier to scan.
- [ ] Rewrite one long comment in Nicky's preferred voice.
- [ ] Review the diff together before touching more files.

Best time: immediately after the style note.

Reason: the broad cleanup should be based on an example that feels right to the actual author.

## Phase 3: Mechanical Spacing Cleanup

This phase should avoid behavior changes.

Small tasks:

- [ ] Find aligned declarations, assignments, struct members, arrays, and comments.
- [ ] Clean obvious alignment in `macros/*.C`.
- [ ] Clean obvious alignment in `src/*.cxx`.
- [ ] Clean obvious alignment in `include/*.h`.
- [ ] Clean obvious alignment in `include/plotting/*.h`.
- [ ] Clean obvious alignment in `tests/*.cxx`.
- [ ] Leave tables or arrays aligned only when the alignment genuinely improves reading.
- [ ] **`external/` is off-limits, full stop — see style.md's Hard Rule. Not "avoid bulk formatting," do not touch at all.**
- [ ] Avoid changing generated or third-party files.
- [ ] Re-read the diff file by file.

Examples to fix:

```cpp
bool useNorm        = ...
std::string trig    = ...
std::string noTrig  = ...
```

```cpp
TCanvas* c     = nullptr;
TPad*    main  = nullptr;
TPad*    ratio = nullptr;
```

Best time: after calibration edits.

Reason: this is broad but low risk if it stays mechanical.

## Phase 4: Comment Triage

Sort comments into four categories.

Keep:

- [ ] Comments that explain physics decisions.
- [ ] Comments that explain ROOT ownership/lifetime pitfalls.
- [ ] Comments that explain non-obvious CMS/JEC/JER formats.
- [ ] Comments that mark sections in long ROOT workflows.
- [ ] Comments that explain why a simple-looking change would be wrong.
- [ ] Short comments that match the `CMS_2026PbPb` style and help scan a long event loop.

Shorten:

- [ ] Legalistic comments with repeated "must," "never," and "always."
- [ ] Multi-line explanations that can become one or two lines.
- [ ] Comments that include three examples when one is enough.
- [ ] Comments that narrate the next line of code.
- [ ] Comments that would be acceptable in docs but too loud in source.

Move to docs:

- [ ] lxplus/CMSSW environment history.
- [ ] annobin failure details.
- [ ] ROOT/libCling mismatch history.
- [ ] RPATH/RUNPATH explanation.
- [ ] SCRAM package path constraints.
- [ ] Condor sandbox packaging rationale.
- [ ] Long CLI philosophy.

Delete:

- [ ] Comments that restate function names.
- [ ] Comments that are stale.
- [ ] Comments that explain removed behavior.
- [ ] Comments that are only there to sound safe.
- [ ] Comments whose only purpose is to make code seem carefully justified.

Best time: after spacing cleanup, before structural build changes.

Reason: comment cleanup changes the repo's voice without touching behavior.

## Phase 5: CMake Cleanup, Behavior-Preserving

First pass: make the current `CMakeLists.txt` less strange without moving logic.

Small tasks:

- [ ] Replace the long opening failure text with a shorter message.
- [ ] Collapse repeated CMSSW comments into short headings.
- [ ] Remove visual alignment in `set()` and declaration-like blocks.
- [ ] Make comments describe the compatibility issue, not the whole debugging story.
- [ ] Keep all current logic in place.
- [ ] Check that local CMake configure still works.
- [ ] Check that local CMake build still works.

Best time: after comment triage.

Reason: this gives a readable diff and keeps build behavior stable.

Second pass: split compatibility modules.

Small tasks:

- [ ] Create `cmake/CMSSWCompat.cmake`.
- [ ] Move hostname/cmsenv checks if still wanted.
- [ ] Move annobin spec cleanup.
- [ ] Move CMSSW ROOT discovery/RPATH pieces.
- [ ] Create `cmake/RootSetup.cmake` only if it makes the top-level clearer.
- [ ] Keep top-level CMake focused on project, targets, tests, and install/build outputs.
- [ ] Reconfigure from a clean build directory.
- [ ] Build.
- [ ] Run tests.

Best time: after the behavior-preserving CMake cleanup.

Reason: moving CMake logic has real risk, so it should be a separate commit.

## Phase 6: SCRAM/CMake Boundary Cleanup

Small tasks:

- [ ] Decide whether SCRAM support is first-class or compatibility-only.
- [ ] Keep the `bin/*.cc` wrappers if SCRAM needs them.
- [ ] Add a short explanation of why `bin/runAsymmetry.cc` includes `../macros/runAsymmetry.C`.
- [ ] Remove duplicated include path prose from README if it belongs in SCRAM notes.
- [ ] Review `BuildFile.xml` for hard-coded `Analysis/L2Residuals` paths.
- [ ] Decide whether package-path independence is worth solving now.
- [ ] If not solving it, document it once and tersely.

Best time: after CMake cleanup, before Condor refactor.

Reason: Condor depends on which build output is considered canonical.

## Phase 7: Condor Script Voice Cleanup

First pass: comments only, no logic changes.

Small tasks:

- [ ] Replace the giant header with concise usage.
- [ ] Move long examples to README or `readme/condor.md`.
- [ ] Keep the important `-output`, `-alltxt`, `-filelists`, `-config`, `-tag`, and `-nosubmit` behavior visible.
- [ ] Shorten the explanation of `[condor.filelist_modes]`.
- [ ] Shorten the `analysis_config.toml` transfer explanation.
- [ ] Shorten the CMSSW source explanation.
- [ ] Shorten the data staging explanation.
- [ ] Keep comments near tricky shell code.
- [ ] Use `CMS_2026PbPb/executable/Condor/make_condor.sh` as the brevity reference, not as a feature reference.
- [ ] Check `bash -n condor/make_condor.sh`.

Best time: after build boundary cleanup.

Reason: comment-only cleanup is low risk and will already make the script feel much less synthetic.

Second pass: split into helpers.

Small tasks:

- [ ] Extract `usage`.
- [ ] Extract `parse_args`.
- [ ] Extract `require_config`.
- [ ] Extract `lookup_filelist_mode`.
- [ ] Extract `require_cmssw`.
- [ ] Extract `choose_binary`.
- [ ] Extract `resolve_filelists`.
- [ ] Extract `normalize_output_dir`.
- [ ] Extract `prepare_submission_sandbox`.
- [ ] Extract `write_submit_header`.
- [ ] Extract `append_job_to_submit_file`.
- [ ] Extract `submit_filelist`.
- [ ] Keep the final main flow short and readable.
- [ ] Check `bash -n condor/make_condor.sh`.
- [ ] Run with `-nosubmit` on a tiny filelist if possible.

Best time: after comment-only Condor cleanup.

Reason: splitting shell code is a behavior risk. It deserves its own pass.

## Phase 8: Macro Front Door Cleanup

Small tasks:

- [ ] Shorten the top comments in `macros/runAsymmetry.C`.
- [ ] Shorten the top comments in `macros/runCalibration.C`.
- [ ] Shorten the top comments in `macros/runTextFile.C`.
- [ ] Shorten the top comments in `macros/runResponse.C`.
- [ ] Shorten the top comments in `macros/runPlotting.C`.
- [ ] Keep compiled usage strings accurate.
- [ ] Keep interpreted ROOT usage, but make it compact.
- [ ] Do not repeat the full `CommandLine` philosophy in every macro.
- [ ] Put the full CLI explanation in README or `readme/cli.md`.

Best time: after config policy is documented once.

Reason: macros are user-facing, so they need clarity, not essay headers.

## Phase 9: README And Docs Consolidation

**Cancelled 2026-07-04 — Nicky's call ("basically forget and remove").**
He's assembling his own README from the `readme/README_v*.md` snapshots
instead. Left below for history only; not part of the active plan.

Small tasks:

- [ ] Decide whether root `README.md` should be a quick start or a full manual.
- [ ] Move detailed build notes to `readme/build.md`.
- [ ] Move Condor details to `readme/condor.md`.
- [ ] Move style choices to `readme/style.md`.
- [ ] Move CLI/config policy to `readme/cli.md` if README is too long.
- [ ] Keep root README focused on what the repo does and how to run the main workflow.
- [ ] Remove repeated explanations from root README after they move.
- [ ] Keep physics workflow explanations where users will actually look for them.

Best time: after code comments are shortened.

Reason: moved comments need a home before they disappear from code.

## Phase 10: Source Code Voice Cleanup

Small tasks:

- [ ] `src/RunAsymmetry.cxx`: clarify section labels and the trigger cone lookup.
- [ ] `src/AnalysisConfig.cxx`: shorten config policy comments.
- [ ] `src/TextFileWriter.cxx`: preserve real JEC/JER format notes, shorten repeated policy.
- [ ] `src/CalibrationExtractor.cxx`: preserve fit/statistics notes, shorten section banners if desired.
- [ ] `src/ResponseExtractor.cxx`: review for repeated extraction policy.
- [ ] `include/TextFileWriter.h`: trim API comment to what callers need; move deep JER explanation to docs or `.cxx`.
- [ ] `include/Utilities.h`: shorten the long `FoldEtaAxis` entry-count comment but preserve the ROOT pitfall.
- [ ] `include/plotting/Utilities.h`: remove column alignment in small structs.
- [ ] `include/JetSelection.h`: decide whether the bottom "How to use" block belongs in docs.
- [ ] `tests/*.cxx`: keep descriptive test names but shorten comments that restate assertions.
- [ ] For each edited source file, ask: would this comment density look out of place in `CMS_2026PbPb`?

Best time: after style calibration and docs consolidation.

Reason: source cleanup should be informed by what documentation now owns.

## Phase 11: Tests And Safety Checks

Run after each risky phase.

Small tasks:

- [ ] `bash -n condor/make_condor.sh`
- [ ] `bash -n condor/runtime_wrapper.sh`
- [ ] `cmake -S . -B build`
- [ ] `cmake --build build`
- [ ] `ctest --output-on-failure`
- [ ] If on lxplus/CMSSW: `scram b -j4`
- [ ] If on lxplus/CMSSW: dry-run Condor with `-nosubmit`
- [ ] If ROOT is available: load each macro in batch mode.
- [ ] If possible: run one tiny end-to-end workflow fixture.

Best time: continuously, but especially after CMake, SCRAM, and Condor changes.

Reason: the cleanup is intended to be behavioral no-op unless explicitly stated.

## Phase 12: Commit Message History (optional, high-risk, not scheduled)

Not started. Written up 2026-07-04 at Nicky's request so the steps exist if
this is ever picked up, without committing to doing it now.

The problem: of 176 commits, the bulk read as uniform Capitalized-Imperative-
Mood ("Rework JES/JER response schema to eta_reco-binned, 3-variant
comparison"), a stylistic tell that assistance was involved in writing them.
The earliest ~15 commits ("making readme a working version", "initial repo
structure") are Nicky's own unassisted style and don't need touching.

What a full rewrite would take:

- [ ] `git filter-repo --message-callback` (or an interactive `rebase -i`
      with `reword` on ~160 commits) to rewrite message text repo-wide.
      `filter-repo` is the only realistic option at this commit count;
      `rebase -i` with 160 rewords by hand isn't practical.
- [ ] Every annotated tag (`v0.1.0` through `v0.6.0`) points at a commit
      hash that changes once its message (or any ancestor's) is rewritten.
      Each tag needs deleting and recreating against the new hash, locally
      and on the remote (`git push origin :refs/tags/vX.Y.Z` then
      re-push).
- [ ] Force-push `main` to `origin` afterward (`git push --force-with-lease`).
      Any other clone (lxplus included) is now diverged and needs the same
      diff-then-reset treatment as the [[project_lxplus_stale_clone_2026-07-03]]
      incident, not a plain `git pull`.
- [ ] Anyone else with a clone or fork loses the ability to fast-forward;
      out of scope here since this repo has no other known clones besides
      lxplus, but worth confirming before doing this for real.
- [ ] Decide up front what the new message text should say per commit —
      "vary tense/tone" isn't mechanical the way em-dash stripping is, so
      this is a real per-commit judgment pass, not a script.

Given the prior `filter-repo` rewrite (the AI-attribution scrub, see
[[session_2026-07-01b]]) already caused a real divergent-history incident
on lxplus, doing this again should be a deliberate, scheduled action with
time set aside to fix any clone that drifts, not a quick follow-on to a
comment cleanup pass.

## Suggested Commit Order

1. `docs: add local style guide`
2. `style: calibrate source cleanup on selected files`
3. `style: remove visual alignment spacing`
4. `docs: add build and condor notes`
5. `style: shorten comments without behavior changes`
6. `build: simplify CMake compatibility comments`
7. `build: move CMSSW compatibility into CMake module`
8. `condor: shorten submission script docs`
9. `condor: split submission script into helpers`
10. `docs: consolidate CLI and config policy`

## Questions For Nicky — Resolved 2026-07-04

See `style.md` for the full write-up. Answers, for reference:

- Voice: casual, mostly objective, first-person okay where it fits, but the emphasis is brevity — most functions need zero comments.
- `CMS_2026PbPb` confirmed as the style reference, explicitly preferred over any prior repo.
- Terse fragments over complete sentences in code comments.
- Final source: almost no comments outside section labels and real pitfalls — confirmed.
- Section labels (`// event loop`, `// trigger`, `// jet ID`) stay as landmarks.
- **No physics explanations in source comments at all** — a short label naming *what* (e.g. `// alpha extrapolation`) is fine; prose justifying *why* moves to docs, always.
- ROOT macro headers: minimal, but keep one compact usage example — confirmed.
- One-line `if`/`for` bodies, including guard clauses, always use `{ }` — no bare one-liners.
- Casts stay C-style (`( Type )value`), matching existing code and ROOT/CMS convention.
- Variable names (`fi`, `fo`, `dcs`, `js`) are preserved — that's the author's real style, not an error.
- "must/never/always" avoided unless it's a real correctness rule — confirmed.
- All debugging/incident history moves into one docs file, short reference style rather than narrative.
- Tests: style/comments only, no behavior changes — confirmed.
- **Vendored/external code (`external/`) is completely excluded from cleanup, permanently — not a phase-3 nicety, a hard rule.** See `style.md`'s "Hard rule" section.

## A Good Way To Accomplish This

Start with a short pairing pass on one source file. `src/RunAsymmetry.cxx` is probably the best because it includes both your human style and some comments that need translation from uncertainty into meaning.

Then make one small style guide file and use it as the reference for every later sweep. The guide should explicitly say that comment brevity is the goal, not just "better comments." After that, do broad mechanical spacing cleanup separately from comment cleanup. Only after the voice feels right should the build and Condor structure change.

The main discipline is to keep each commit boring:

- spacing only
- comments only
- docs only
- CMake move only
- Condor refactor only

That keeps review honest and makes it easy to back out any cleanup that accidentally changes behavior.
