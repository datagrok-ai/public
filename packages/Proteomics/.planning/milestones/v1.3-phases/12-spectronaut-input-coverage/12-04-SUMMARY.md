---
phase: 12-spectronaut-input-coverage
plan: 04
subsystem: parsers
tags: [spectronaut, streaming, proteomics, ux, gap-closure]

# Dependency graph
requires:
  - phase: 12-spectronaut-input-coverage
    provides: "12-02 — parseSpectronautStream + handleFields closure + finalizeSpectronaut shared tail; 12-03 — committed duckdb golden sidecar + stream-vs-text per-cell parity tests (the byte-identical-output proof this plan must keep green)"
provides:
  - "src/parsers/spectronaut-parser.ts — handleFields returns a discriminated LineOutcome ('kept' | 'malformed' | 'filtered'); two counters (malformed, filtered) replace the single skipped counter; the completion + in-progress messages reference ONLY the genuine-malformed count; by-design-filtered rows (CON__/REV__, q>threshold) are SILENT, matching the text path"
  - "src/tests/spectronaut-parser.ts — 2 new regression tests inside category('Spectronaut') (26 total): a decoy-only + q>threshold-only fixture emits NO 'malformed' message; a too-few-fields row still surfaces the genuine 'skipped N malformed line(s)' message"
affects: ["Phase 12 UAT Gap 1 closed — the false 'skipped ~millions of malformed line(s)' corruption signal on a provably-correct multi-GB precursor import is removed; the standing grok test regression net now locks the malformed-vs-filtered distinction against future re-collapse"]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Discriminated string-literal union (LineOutcome) replaces a boolean return to carry a malformed-vs-by-design-filtered distinction without altering any duckdb-parity predicate"
    - "User-facing-contract test: per-line counters are function-local and not returned, so the regression test spies grok.shell.info (with finally-restore so a thrown assertion cannot leak the stub — T-12-18) and asserts the completion-message contract directly"

key-files:
  created: []
  modified:
    - packages/Proteomics/src/parsers/spectronaut-parser.ts
    - packages/Proteomics/src/tests/spectronaut-parser.ts

key-decisions:
  - "filtered counter is incremented but never surfaced (added `void filtered;` to satisfy TS strict noUnusedLocals while preserving the plan's locked 'filtered is SILENT' UX decision) — by-design drops match pivotSpectronaut's identical silent drops"
  - "Each caller site (main while-loop + trailing-partial-line flush) switches on the LineOutcome discriminant; the trailing-flush branch was restructured from a chained `else if` into an explicit `if (!headerParsed) {...} else { switch }` so the discriminant switch is unambiguous (no predicate / control-flow change)"
  - "TDD-gate model: this plan's per-task tdd=\"true\" is satisfied as the phase-level refactor(...)+test(...) pair the 12-03 SUMMARY documented — Task 1 is a zero-output-change categorization refactor (Plan-03 golden + per-cell parity tests stay green proving byte-identical aggregation), Task 2 is the regression net pinning the new message contract; no RED fail-fast trigger applies"
  - "The machine-readable failed-count gate's regex has no match on a fully-green grok test run on this server (release/1.27.3 prints 'Tests passed.' + 'Exiting with code 0', no 'Failed: <n>' token); honored the plan's intent via the deterministic clean-pass signal — the identical resolution Plan 01/02/03 SUMMARYs documented"

patterns-established:
  - "A categorization/UX fix to a streaming aggregator is provably output-neutral when the prior plan's committed golden + per-cell parity tests stay green and the diff is constrained to return values + counters + message strings (predicates explicitly forbidden to change, verified by the diff)"

requirements-completed: [R1]

# Metrics
duration: 5min
completed: 2026-05-15
---

# Phase 12 Plan 04: Spectronaut Malformed-vs-Filtered Categorization Summary

**`handleFields` now returns a discriminated `'kept' | 'malformed' | 'filtered'` outcome and the streaming Spectronaut import surfaces "malformed line(s)" ONLY for genuinely unparseable rows — CON__/REV__ decoy and numeric-q>threshold drops are tracked separately and SILENT, exactly as the text path (`pivotSpectronaut`) already is — closing Phase 12 UAT Gap 1's false corruption signal on a provably-correct multi-GB precursor import, with two regression tests pinning the new behavior and the committed Plan-03 duckdb golden / per-cell parity tests proving aggregation output is byte-identical.**

## Performance

- **Duration:** ~5 min
- **Started:** 2026-05-15T18:45:56Z
- **Completed:** 2026-05-15T18:50:47Z
- **Tasks:** 2
- **Files modified:** 2 (0 created, 2 modified)

## Accomplishments

- **Task 1 (`refactor`):** Added a `LineOutcome = 'kept' | 'malformed' | 'filtered'` type alias next to `AggRow`. Changed the `handleFields` closure return type from `boolean` to `LineOutcome`, mapping each existing branch per the locked categorization **without touching any predicate**: `f.length < expectedFields` → `'malformed'`, `!protein` → `'malformed'`, `CON__`/`REV__` → `'filtered'`, `q !== null && q > qValueThreshold` → `'filtered'`, `qty === null && q === null` → `'malformed'`, aggregated → `'kept'`. Replaced the single `let skipped` counter with `let malformed` + `let filtered`. Updated both caller sites (main while-loop and trailing-partial-line flush) to `switch` on the discriminant. The in-progress `pi.update` message and the completion `grok.shell.info` now reference ONLY the `malformed` count (verbatim wording `Spectronaut import: skipped ${malformed} malformed line(s)`, guarded by `malformed > 0`); `filtered` is tracked but never surfaced. Diff confirmed: zero predicate-expression lines changed; aggregation/pivot/finalize/text-path code untouched.
- **Task 2 (`test`):** Added a `streamCapturingInfo` helper (spies `grok.shell.info`, restores the original in a `finally`) and two tests inside `category('Spectronaut')`: `by-design-filtered rows are not counted malformed` (a CON__ + REV__ + numeric-q-0.05 fixture with no truncation produces the correct 2-protein DataFrame and emits ZERO `grok.shell.info` matching `/malformed/i`) and `truncated line is counted malformed` (a valid fixture plus one `file\tCondA\t1` truncated line still surfaces a `/skipped \d+ malformed line\(s\)/` message). Added `import * as grok from 'datagrok-api/grok';`. `makeLongFormatTsv`, `streamTsv`, and the 24 prior tests are byte-unchanged (86 insertions, 0 deletions).
- Verified GREEN on localhost: `npx tsc --noEmit` 0 errors after each task; `grok test --category "Spectronaut" --host localhost` → **26 passed** (24 prior + 2 new); full Proteomics package suite fully green (`Parsers 15, Experiment Setup 3, Normalization 8, Imputation 9, DE 6, Generic 11, QC 7, Enrichment 7, Enrichment Viz 8, Spectronaut 26, SpectronautCandidates 18, SpectronautCandidates E2E 1, FragPipe 11, FragPipe E2E 1` — `Tests passed.`, exit 0). The Plan-03 `stream path matches text path` and `streaming output equals duckdb golden` tests stayed green — proof Task 1 left aggregation output byte-identical.

## Task Commits

Each task was committed atomically (intended file only):

1. **Task 1: split spectronaut streaming skipped counter into malformed vs by-design-filtered** — `755b0bfc63` (refactor)
2. **Task 2: pin spectronaut by-design-filtered rows are not reported as malformed** — `a0c5610e5b` (test)

**Plan metadata:** committed with this SUMMARY (docs) — worktree mode: SUMMARY only (STATE.md / ROADMAP.md owned by the orchestrator post-wave).

## Files Created/Modified

- `packages/Proteomics/src/parsers/spectronaut-parser.ts` — `+LineOutcome` type alias; `handleFields` returns the discriminated outcome (predicates unchanged); `skipped` → `malformed` + `filtered`; both caller sites switch on the discriminant; in-progress + completion messages reference only `malformed`; `void filtered;` so the unused-but-intentionally-silent counter passes TS strict. Aggregation/pivot/finalize/text-path untouched (60 insertions, 18 deletions; no predicate line changed).
- `packages/Proteomics/src/tests/spectronaut-parser.ts` — `+import * as grok`; `+streamCapturingInfo` helper; `+2` tests inside `category('Spectronaut')`. `makeLongFormatTsv`, `streamTsv`, and the 24 prior tests byte-unchanged (86 insertions, 0 deletions).

## Decisions Made

- **`void filtered;` for the silent-but-tracked counter.** The plan mandates `filtered` be tracked yet NEVER surfaced. TS strict (`noUnusedLocals`) would flag a counter that is only written, never read. `void filtered;` documents the intent (tracked for completeness, silent by locked UX decision) and keeps `tsc` clean without weakening the "filtered is silent" contract or referencing it in any message.
- **Trailing-flush control-flow restructure (no behavior change).** The original trailing-partial-line flush used a chained `if (!headerParsed) parseHeader(line); else if (handleFields(...)) rowCount++; else skipped++;`. To switch on the 3-member discriminant unambiguously, this became `if (!headerParsed) { parseHeader(line); } else { switch (handleFields(...)) { ... } }`. The header-guard predicate and the `parseHeader` call are byte-identical; only the post-guard increment changed from boolean to switch — same as the main loop.
- **Clean-pass gate fallback.** The plan's `<automated>` gate parses `/[Ff]ailed[:\s]+(\d+)/`; a fully-green run on this server (release/1.27.3) emits no such token, printing `Tests passed.` + `Exiting with code 0`. Kept the regex as primary; on no-match, evaluated the deterministic clean-pass signal (`Tests passed.` AND `Exiting with code 0` AND no `[1-9]\d* failed`) — the identical resolution Plan 01/02/03 SUMMARYs documented for this runner. No code/test change.

## Deviations from Plan

None — plan executed exactly as written. The two restructure/`void` decisions above are mechanical consequences of TS strict mode and the 3-member switch (the plan explicitly anticipated "Keep the surrounding control flow ... exactly as-is" for the predicate/aggregation; the header-guard predicate and all `handleFields` predicates ARE byte-identical). No Rule 1–4 trigger, no architectural change, no scope creep, no auth gate. Scope held to exactly `src/parsers/spectronaut-parser.ts` and `src/tests/spectronaut-parser.ts` (obs-A / obs-B and `experiment-setup.ts` / `normalization.ts` NOT touched, per the plan's scope guard).

## TDD Gate Compliance

This plan's frontmatter sets `tdd="true"` per task, but the plan body explicitly structures it (and the `<done>` clauses mandate) a code task committed as `refactor(...)` followed by a regression-test task committed as `test(...)` — the identical phase-level model the 12-03 SUMMARY documented for this phase set. The RED→GREEN cycle for the streaming **production code** (`parseSpectronautStream` / `handleFields`) was Plan 02 (Wave 1). Plan 04 Task 1 is a categorization + user-message refactor with **zero aggregation-output change** — the Plan-03 committed `streaming output equals duckdb golden` and `stream path matches text path` per-cell parity tests stayed green, which is the byte-identical-output proof. Task 2 is the regression net pinning the new message contract; both its tests passed GREEN immediately against the Task-1 implementation. No RED fail-fast trigger occurred (Task 2 asserts an already-implemented behavior — the intended flow for a gap-closure regression-net task). The phase-level RED/GREEN gate is satisfied across the wave (Plan 02 `feat(...)` + this plan's `refactor(...)` + `test(...)`).

## Issues Encountered

- This freshly-created git worktree had no `node_modules` (same as Plan 01/02/03). Resolved per the plan's `<conventions>` by symlinking the sibling `proteomics` worktree's `node_modules` — verified the sibling is at the **exact same base commit** `4d04e1f12c6f16c380cb23f2e74352214494f3ce` and this plan changes zero dependencies (only the two intended source/test files; package.json/lockfile untouched). The symlink was added to the common git-dir `info/exclude` (`packages/Proteomics/node_modules`), never `git add`ed, confirmed absent from every commit (`git diff --name-only 4d04e1f..HEAD` showed only the two intended files) and `git check-ignore`-excluded.

## Threat Surface Scan

The plan's `<threat_model>` (T-12-15..T-12-18) covers exactly this surface and the implemented work mitigates as specified: T-12-15 (re-collapse / re-label regression) — Task 2's `by-design-filtered rows are not counted malformed` spies `grok.shell.info` and fails loudly if any "malformed" message fires for a purely-filtered fixture, `truncated line is counted malformed` independently locks the genuine category still surfaces, both in the standing `grok test` net; T-12-16 (predicate tamper) — the Plan-03 duckdb-golden + per-cell parity tests stayed green and the Task-1 diff has zero predicate-expression changes (verified); T-12-17 (info disclosure) — accepted, the message remains a bounded integer + fixed string, no row content/path/PII; T-12-18 (leaked spy stub) — `streamCapturingInfo` restores the original `grok.shell.info` in a `finally`, and the full-category green gate would catch a leaked stub. No new security-relevant surface: no network endpoints, no auth, no schema changes, no runtime file writes — only client-side categorization of synthetic fixtures in the test runner. CSV/TSV formula-injection out of scope (the message is a plain integer count via `grok.shell.info`; test assets consumed only by `grok test`). No threat flags.

## Self-Check: PASSED

- Files exist on disk: `packages/Proteomics/src/parsers/spectronaut-parser.ts`, `packages/Proteomics/src/tests/spectronaut-parser.ts` (both modified, present).
- Both task commits present in git history: `755b0bfc63` (refactor), `a0c5610e5b` (test).
- `npx tsc --noEmit` 0 errors; `grok test --category "Spectronaut" --host localhost` 26 passed; full Proteomics package suite green on localhost (`Tests passed.`, exit 0); Plan-03 golden / per-cell parity tests green (byte-identical aggregation confirmed).

---
*Phase: 12-spectronaut-input-coverage*
*Completed: 2026-05-15*
