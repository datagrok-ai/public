---
phase: 12-spectronaut-input-coverage
plan: 01
subsystem: testing
tags: [spectronaut, streaming, proteomics, parser, fixture, duckdb, golden-test]

# Dependency graph
requires:
  - phase: 12-spectronaut-input-coverage
    provides: "12-RESEARCH.md duckdb /tmp script + parser parity contract"
provides:
  - "tools/spectronaut-aggregate.{sql,sh} — committed D-05 manual fallback + D-04 equivalence oracle (gene/accession carry-along terms dropped, DMD/WT flip documented reference-only)"
  - "tools/generate-spectronaut-precursor-fixture.mjs — synthetic precursor-level fixture generator"
  - "tools/derive-precursor-golden-sidecar.mjs — verbatim golden→JSON transcriber (no aggregation)"
  - "files/demo/spectronaut-hye-precursor.tsv — D-01-signature fixture exercising every R2 filter branch"
  - "files/demo/spectronaut-hye-precursor-golden.tsv — duckdb-derived D-04 oracle"
  - "files/demo/spectronaut-hye-precursor-golden.json — in-test oracle when committed-file reads are unavailable"
  - "Extended makeLongFormatTsv emitting precursor/decoy/edge-q-value rows"
affects: [12-02 (streaming aggregator), 12-03 (equivalence golden test)]

# Tech tracking
tech-stack:
  added: [duckdb v1.3.0 CLI (dev/CI dependency, not runtime)]
  patterns:
    - "Three-step regen chain (fixture → duckdb golden → JSON sidecar) documented and order-enforced in README"
    - "Sidecar deriver does pure transcription, never re-aggregation, so the equivalence test stays pinned to real duckdb output"

key-files:
  created:
    - packages/Proteomics/tools/spectronaut-aggregate.sql
    - packages/Proteomics/tools/spectronaut-aggregate.sh
    - packages/Proteomics/tools/generate-spectronaut-precursor-fixture.mjs
    - packages/Proteomics/tools/derive-precursor-golden-sidecar.mjs
    - packages/Proteomics/files/demo/spectronaut-hye-precursor.tsv
    - packages/Proteomics/files/demo/spectronaut-hye-precursor-golden.tsv
    - packages/Proteomics/files/demo/spectronaut-hye-precursor-golden.json
  modified:
    - packages/Proteomics/src/tests/spectronaut-parser.ts
    - packages/Proteomics/files/demo/README.md

key-decisions:
  - "Drop the two any_value gene/accession carry-along SELECT terms from the committed SQL with a documented DIVERGENCE comment (data lacks them; parser never consumes them; would Binder-Error duckdb over the no-Genes fixture)"
  - "Keep the DMD/WT R.Condition flip but rewrite its comment as a loud REFERENCE-FILE-ONLY warning — structural no-op on CondA/CondB so the same script doubles as fallback and oracle"
  - "Phrase the SQL DIVERGENCE comment without the literal column-name strings so the acceptance grep (no PG.Genes/PG.ProteinAccessions anywhere) holds while the rationale is still documented"
  - "Reuse the sibling proteomics worktree's node_modules via a transient (untracked) symlink to run the regression gate; same base commit, deps unchanged by this plan"

patterns-established:
  - "Golden oracle = real duckdb output committed verbatim; JSON sidecar transcribed FROM it, never recomputed"
  - "Synthetic fixtures use CondA/CondB so the reference-only DMD/WT flip is provably a structural no-op"

requirements-completed: [R3, R5, R6]

# Metrics
duration: 7min
completed: 2026-05-15
---

# Phase 12 Plan 01: Spectronaut Wave-0 Test Infrastructure Summary

**Committed duckdb aggregation oracle + synthetic precursor fixture + verbatim JSON sidecar + extended makeLongFormatTsv, establishing the equivalence-test scaffold every later wave depends on with zero production-code changes and zero Spectronaut test regressions.**

## Performance

- **Duration:** ~7 min
- **Started:** 2026-05-15T14:22:30Z
- **Completed:** 2026-05-15T14:28:54Z
- **Tasks:** 3
- **Files modified:** 9 (7 created, 2 modified)

## Accomplishments
- Relocated the proven `/tmp` duckdb aggregation into `packages/Proteomics/tools/` as the D-05 manual fallback and D-04 oracle, with the unused gene/accession carry-along terms dropped (documented divergence) and the DMD/WT flip rewritten as a loud reference-only warning so the streaming TS aggregator will not port it.
- Extended `makeLongFormatTsv` in place with the D-01 precursor signature columns (`EG.ModifiedPeptide`/`FG.Charge`/`FG.Id`) plus an `extraRows` opt for decoy/edge-q-value injection — all pre-existing Spectronaut tests stay green (17 + 18 + 1 = 36, `Tests passed.`, exit 0).
- Generated the synthetic precursor-level fixture (41 protein groups incl. `CON__`/`REV__`, CondA/CondB × 3 replicates, group-constant quantity, no gene/accession columns) exercising every Spectronaut filter branch in one file.
- Produced the duckdb golden (39 proteins × 6 samples, decoys filtered, conditions unchanged — proving the flip no-op) and the byte-stable JSON sidecar transcribed verbatim from it; documented all three artifacts + the regen-order contract + flip caveat in `files/demo/README.md`.

## Task Commits

Each task was committed atomically:

1. **Task 1: Relocate duckdb aggregation script into tools/** - `a6a4bc98e1` (chore)
2. **Task 2: Extend makeLongFormatTsv + add precursor fixture generator** - `08eb824ffe` (test)
3. **Task 3: Generate golden, derive JSON sidecar, document fixtures** - `1b7b592683` (test)

**Plan metadata:** committed with this SUMMARY (docs)

## Files Created/Modified
- `packages/Proteomics/tools/spectronaut-aggregate.sql` - duckdb aggregation; gene/accession carry-along terms dropped (documented DIVERGENCE), DMD/WT flip kept but marked REFERENCE-FILE-ONLY
- `packages/Proteomics/tools/spectronaut-aggregate.sh` - dirname-resolved wrapper with documented manual-fallback usage note
- `packages/Proteomics/tools/generate-spectronaut-precursor-fixture.mjs` - emits the synthetic precursor fixture
- `packages/Proteomics/tools/derive-precursor-golden-sidecar.mjs` - pure verbatim golden→JSON transcription (no aggregation)
- `packages/Proteomics/files/demo/spectronaut-hye-precursor.tsv` - D-01-signature fixture, every R2 branch
- `packages/Proteomics/files/demo/spectronaut-hye-precursor-golden.tsv` - duckdb-derived D-04 oracle (234 rows + header)
- `packages/Proteomics/files/demo/spectronaut-hye-precursor-golden.json` - 234-entry in-test oracle
- `packages/Proteomics/src/tests/spectronaut-parser.ts` - extended makeLongFormatTsv (in place, backward-compatible)
- `packages/Proteomics/files/demo/README.md` - three new sections, regen one-liners, flip caveat, regen-order contract, license rows

## Decisions Made
- Dropped the two `any_value` gene/accession carry-along SELECT terms with a documented `DIVERGENCE FROM /tmp` comment — required so duckdb binds over the no-Genes fixture (Task 1 action, per plan).
- Phrased the SQL divergence comment WITHOUT the literal `PG.Genes`/`PG.ProteinAccessions` strings, because the Task 1 acceptance grep treats *any* occurrence of those strings as failure. Rationale is fully documented via "gene-symbol and protein-accession columns" prose. (See Deviations — adjustment of comment wording to satisfy the verification contract; no behavior change.)
- Reused the sibling proteomics worktree's `node_modules` via a transient untracked symlink to run the `grok test` regression gate (this worktree had none; same base commit, this plan changes no dependencies). Symlink never staged/committed.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] SQL divergence comment reworded to satisfy the Task 1 acceptance grep**
- **Found during:** Task 1 (relocate duckdb script)
- **Issue:** The plan's example DIVERGENCE comment text contained the literal strings `any_value("PG.Genes")`/`any_value("PG.ProteinAccessions")`, but Task 1's `<automated>` check asserts `! grep -q 'PG.Genes'` / `! grep -q 'PG.ProteinAccessions'` against the whole SQL file. Using the example verbatim made the verification fail.
- **Fix:** Reworded the divergence comment to describe the dropped terms as "the two carry-along any_value() SELECT terms for the gene-symbol and protein-accession columns" — same documented rationale, no literal column-name tokens. SQL behavior unchanged (only a comment).
- **Files modified:** packages/Proteomics/tools/spectronaut-aggregate.sql
- **Verification:** Full Task 1 `<automated>` block re-run → `OK`.
- **Committed in:** a6a4bc98e1 (Task 1 commit)

**2. [Rule 1 - Bug] Regression-gate node-script false-negatives on a fully-green grok test run**
- **Found during:** Task 2 (regression gate)
- **Issue:** Task 2's `<automated>` gate greps the `grok test` log for a `Failed: <n>` token and fails when absent. A fully-green `grok test` run emits no `Failed:` token at all (it prints `Tests passed.` + `Exiting with code 0`), so the literal gate would report failure on a passing run.
- **Fix:** Evaluated the true success signal instead: runner exit code 0 + `Tests passed.` line + absence of any `<n> failed` indicator. All three Spectronaut suites passed (17 + 18 + 1). The plan's intent — "no Spectronaut regression" — is satisfied; only the gate's string-matching heuristic was wrong.
- **Files modified:** none (verification-logic correction only)
- **Verification:** `grok test --category "Spectronaut" --host localhost` → `✅ Spectronaut (17 passed)`, `✅ SpectronautCandidates (18 passed)`, `✅ SpectronautCandidates E2E (1 passed)`, `Tests passed.`, `Exiting with code 0`.
- **Committed in:** 08eb824ffe (Task 2 commit — test helper change that the gate validated)

---

**Total deviations:** 2 auto-fixed (1 blocking, 1 bug)
**Impact on plan:** Both are adjustments to satisfy/interpret the plan's verification scripts correctly; neither changes the intended SQL behavior, fixture shape, or test outcome. No scope creep. The plan's "19 existing tests" figure was approximate — the Spectronaut category actually contains 36 tests across 3 suites; all pass.

## Issues Encountered
- This freshly-created git worktree had no `node_modules`, so `grok test` could not run directly. Resolved by symlinking the sibling `proteomics` worktree's `node_modules` (same base commit `2486f1411e`, no dependency changes in this plan). The symlink was never staged and was removed after the gate passed.

## Next Phase Readiness
- Wave-0 oracle chain complete: Plan 02 (streaming aggregator) and Plan 03 (equivalence golden test) now have the committed duckdb oracle, the routing fixture, the golden `.tsv`, and the in-test JSON sidecar they require.
- The Plan-01 internal contradiction (the `/tmp` SQL selected gene/accession columns the fixture cannot supply) is reconciled and verified end-to-end: `tools/spectronaut-aggregate.sh` exits 0 over the no-Genes fixture.
- No production parser/handler code was written (per plan); no blockers.

## Self-Check: PASSED

All 9 created/modified files exist on disk and all 4 commits (a6a4bc98e1, 08eb824ffe, 1b7b592683, 412657297f) are present in git history.

---
*Phase: 12-spectronaut-input-coverage*
*Completed: 2026-05-15*
