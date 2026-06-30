---
phase: 12-spectronaut-input-coverage
plan: 03
subsystem: testing
tags: [spectronaut, streaming, proteomics, test, regression, golden-test, duckdb-parity]

# Dependency graph
requires:
  - phase: 12-spectronaut-input-coverage
    provides: "12-01 — precursor fixture + duckdb golden + JSON sidecar + extended makeLongFormatTsv; 12-02 — parseSpectronautStream + exported sniffIsPrecursor + finalizeSpectronaut shared tail"
provides:
  - "src/tests/spectronaut-parser.ts — 7 new streaming tests (R1 smoke, D-01 routing, stream-vs-text per-cell, streaming-vs-duckdb-golden, filter parity, tag/group, tools-presence) inside the existing category('Spectronaut') suite; 17 prior tests untouched"
  - "webpack.config.js asset/source rule for .sql/.sh — inlines the committed tools/spectronaut-aggregate.{sql,sh} into the test bundle as the R5 drift guard (tools/ is not runtime-readable via _package.files)"
  - "src/global.d.ts — ambient module declarations for the .sql/.sh raw-text imports (keeps tsc clean)"
affects: ["12 phase regression net — locks streaming↔text↔duckdb-golden parity and D-01 routing against future drift"]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Per-cell numeric equivalence (within 1e-3) keyed by protein id between parseSpectronautStream and parseSpectronautText — catches the documented first-encountered-vs-max divergence a structural-only check would miss"
    - "Golden test loads the committed JSON sidecar via _package.files (the demo-file convention) and asserts exact (protein × sample) key-set equality + per-key quantity within 1e-3 — never re-aggregates, no conditional fallback"
    - "webpack asset/source to pin non-deployed tools/ content into the test bundle as a deletion/rename/content-change drift guard"

key-files:
  created:
    - packages/Proteomics/src/global.d.ts
  modified:
    - packages/Proteomics/src/tests/spectronaut-parser.ts
    - packages/Proteomics/webpack.config.js

key-decisions:
  - "tools/ is not deployed under files/ and is unreachable via _package.files at runtime; resolved (Rule 3) by a narrowly-scoped webpack asset/source rule (.sql/.sh) + ambient module decls so the committed fallback content is pinned into the test bundle — the R5 drift guard the plan intends, with the assertion logic kept in src/tests/spectronaut-parser.ts"
  - "Golden sidecar read via _package.files.readAsText('demo/...') + JSON.parse rather than a static JSON import (tsconfig has no resolveJsonModule; matches the existing spectronaut-candidates-e2e / fragpipe-e2e committed-file convention)"
  - "The machine-readable failed-count gate's regex has no match on a fully-green grok test run (which prints 'Tests passed.' + 'Exiting with code 0' and no 'Failed: <n>' token); honored the plan's intent via the clean-pass signal (Tests passed. + exit 0 + no '<n> failed') — the identical resolution Plan 01 & Plan 02 SUMMARYs documented"

patterns-established:
  - "Streaming regression net: every R1/R2/R3/R5/R6 behavior has a named grok test; result gate is deterministic (clean-pass signal, not brittle substring)"

requirements-completed: [R1, R2, R3, R5, R6]

# Metrics
duration: 6min
completed: 2026-05-15
---

# Phase 12 Plan 03: Spectronaut Streaming Regression Net Summary

**Seven new `grok test` cases lock the Plan-02 streaming path against the Plan-01 fixture/golden/sidecar — an R1 smoke test, the D-01 `sniffIsPrecursor` routing branch in isolation, per-(protein × raw sample) numeric equivalence (within 1e-3) to the text path, exact key-set + per-key-quantity equivalence to the verbatim duckdb golden sidecar, every R2 filter edge, the tag/group contract, and a committed-fallback drift guard — with the 17 prior Spectronaut tests untouched and the full package suite green.**

## Performance

- **Duration:** ~6 min
- **Started:** 2026-05-15T14:52:06Z
- **Completed:** 2026-05-15T14:57:56Z
- **Tasks:** 2
- **Files modified:** 3 (1 created, 2 modified)

## Accomplishments
- Added `parseSpectronautStream` + `sniffIsPrecursor` imports and a `streamTsv` / `readDemoFile` / `rawSampleCols` / `proteinRowIndex` helper set to `src/tests/spectronaut-parser.ts`; the 17 prior Spectronaut tests and `makeLongFormatTsv` are byte-unchanged (RESEARCH Open Q1 — they stay locked on `parseSpectronautText`).
- **Task 1 (5 tests):** `streams precursor fixture` (R1 smoke — the exact 12-VALIDATION.md name, streams the committed `files/demo/spectronaut-hye-precursor.tsv`); `sniffIsPrecursor routes by header` (true for a precursor-signature header, false for a PG-only header — the D-01 routing branch in isolation); `stream path matches text path` (structural equality + per-(protein × raw sample column) value equality within 1e-3, catching the first-encountered-vs-max divergence); `streaming filter parity` (CON__/REV__/all-q>0.01 excluded, non-numeric `Profiled`/empty-string q-value included); `streaming tag set and groups`.
- **Task 2 (2 tests):** `streaming output equals duckdb golden` (loads the Plan-01 JSON sidecar via `_package.files`, asserts per-key quantity within 1e-3 AND exact (protein × sample) key-set equality vs the verbatim duckdb-derived sidecar — no conditional fallback, reads fail loudly); `duckdb fallback tooling is committed` (asserts `tools/spectronaut-aggregate.sql` contains `max(TRY_CAST` and `.sh` references `spectronaut-aggregate.sql`).
- Verified GREEN on localhost at every gate: Spectronaut category 24 passed (17 prior + 7 new) after each task, full Proteomics package suite fully green (`Parsers 15, Experiment Setup 3, Normalization 8, Imputation 9, DE 6, Generic 11, QC 7, Enrichment 7, Enrichment Viz 8, Spectronaut 24, SpectronautCandidates 18, SpectronautCandidates E2E 1, FragPipe 11, FragPipe E2E 1` — `Tests passed.`, exit 0). `npx tsc --noEmit` clean (0 errors).

## Task Commits

Each task was committed atomically:

1. **Task 1: streaming smoke / sniffIsPrecursor routing / stream-vs-text per-cell / filter-parity / tag-set tests** - `013236a70b` (test)
2. **Task 2: streaming-vs-duckdb-golden (sidecar) + tools-file-presence tests** - `203cf570f8` (test)

**Plan metadata:** committed with this SUMMARY (docs) — worktree mode: SUMMARY only (STATE.md/ROADMAP.md owned by the orchestrator post-wave).

## Files Created/Modified
- `packages/Proteomics/src/tests/spectronaut-parser.ts` — +`parseSpectronautStream`/`sniffIsPrecursor`/`_package` imports, raw `.sql`/`.sh` imports, 4 helpers, 7 new tests inside `category('Spectronaut')`; the 17 prior tests and `makeLongFormatTsv` untouched.
- `packages/Proteomics/webpack.config.js` — added an `asset/source` module rule for `.sql`/`.sh` so the committed `tools/spectronaut-aggregate.{sql,sh}` are inlined into the test bundle (Rule 3 — see Deviations).
- `packages/Proteomics/src/global.d.ts` — **new**; ambient `declare module '*.sql'`/`'*.sh'` so the raw-text imports type-check (Rule 3 — see Deviations).

## Decisions Made
- **Golden sidecar read path:** `_package.files.readAsText('demo/spectronaut-hye-precursor-golden.json')` + `JSON.parse` instead of a static JSON import — `tsconfig.json` has no `resolveJsonModule`, and this matches the existing `spectronaut-candidates-e2e` / `fragpipe-e2e` committed-file convention exactly.
- **Golden key split:** key form is `${protein}${condition}_${replicate}`; split via `/^(.*?)((?:Cond[A-Za-z]+)_\d+)$/` so a multi-segment protein id stays intact and the condition/replicate suffix is isolated.
- **Key-set equality is built from non-null streamed cells** (matching how the duckdb golden / sidecar is a sparse 39-protein × 6-sample grid, 234 entries) so the streamed set equals the sidecar set exactly.
- **No conditional fallback in the golden test** (the plan's removed "if committed-file reads unavailable…" clause): the sidecar is committed precisely so the comparison is always available; a read failure FAILS the test, it does not degrade it. Verified `! grep -q "if committed-file reads"`.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] `tools/` is not runtime-readable; pinned via webpack `asset/source` + ambient decls**
- **Found during:** Task 2 (tools-file-presence test).
- **Issue:** The plan's `duckdb fallback tooling is committed` test must read `tools/spectronaut-aggregate.{sql,sh}` "via the package files API … same helper convention as the demo-file reads." But `_package.files` maps only the deployed `files/` directory; `tools/` is **not** deployed to the runtime and is unreachable via `_package.files` in the Puppeteer test browser. A direct `_package.files.readAsText('../tools/...')` (or any runtime read) cannot resolve it.
- **Fix:** Added a narrowly-scoped `asset/source` webpack rule for `.sql`/`.sh` so the committed `tools/spectronaut-aggregate.{sql,sh}` are inlined into the **test bundle** as raw strings, and a `src/global.d.ts` with `declare module '*.sql'`/`'*.sh'` so `tsc` type-checks the imports. The assertion logic stays in `src/tests/spectronaut-parser.ts`. This satisfies the plan's R5 *intent* exactly — a deletion/rename/content change of the real `tools/` files fails the build + `grok test` (the drift guard T-12-11 wants) — and is the standard Datagrok pattern for bundling non-code assets (same precedent as the WASM rule's glue imports). The inline `declare module` blocks were first tried in the test file but TS2664 (wildcard module augmentation is illegal inside an ES module) forced the ambient `.d.ts`.
- **Files modified:** `packages/Proteomics/webpack.config.js`, `packages/Proteomics/src/global.d.ts` (new).
- **Verification:** `npx tsc --noEmit` 0 errors; `grok test --category "Spectronaut" --host localhost` → `duckdb fallback tooling is committed` passes; full suite green.
- **Committed in:** `203cf570f8` (Task 2 commit).

**2. [Rule 1 - Bug] Machine-readable failed-count gate false-negatives on a clean run**
- **Found during:** Both task verifies.
- **Issue:** The plan's `<automated>` gate parses `/[Ff]ailed[:\s]+(\d+)/` from the `grok test` log and fails when absent. A fully-green `grok test` run on this server (release/1.27.3) emits **no** `Failed: <n>` token — it prints `Tests passed.` + `Exiting with code 0`. The literal regex gate would therefore report failure on a passing run.
- **Fix:** Kept the regex as the primary path; when it does not match, evaluated the deterministic clean-pass signal instead — `Tests passed.` AND `Exiting with code 0` AND no `[1-9]\d* failed` indicator. This is the **identical** resolution Plan 01 and Plan 02 SUMMARYs documented for the same gate on the same runner; the plan's intent ("failed-count == 0, deterministic, not substring-`grep` of pass/fail") is fully satisfied — the clean-pass signal is itself deterministic and is *not* a substring `grep` for the words "pass"/"fail".
- **Files modified:** none (verification-logic interpretation only; no code/test change).
- **Verification:** Gate evaluated GREEN for Task 1 (`/tmp/12-03-t1.log`), Task 2 (`/tmp/12-03-t2.log`), and the full package suite (`/tmp/12-03-full.log`).
- **Committed in:** n/a (no file change; documented here per the plan's gate semantics).

---

**Total deviations:** 2 auto-fixed (1 blocking, 1 bug).
**Impact on plan:** Deviation 1 adds the minimum build support (one webpack rule + one ambient `.d.ts`) required for the plan's own declared `tools/`-presence test to be possible at all — the assertion and behavior contract are unchanged. Deviation 2 is a verification-interpretation correction (no code change) that matches the precedent both prior plans in this phase set. No scope creep; no behavior, fixture, or production-code change. The plan's "19 existing tests" figure was approximate (the Spectronaut category itself has 17; 24 after this plan; the full Spectronaut grouping across 3 suites is 24 + 18 + 1) — every prior test stays green, consistent with Plan 01 & 02 SUMMARYs.

## TDD Gate Compliance

This plan's frontmatter is `type: tdd`, but the plan body explicitly defines the model: *"TDD plan: tests are written to assert the behavior contract, then run against the already-built Plan 02 implementation."* The RED→GREEN cycle for the streaming **production code** (`parseSpectronautStream`, `sniffIsPrecursor`, `finalizeSpectronaut`) was Plan 02 (Wave 1). Plan 03 is the regression net: both commits are `test(...)` commits that passed GREEN immediately against the already-merged Plan-02 implementation — no `feat(...)` commit is expected or appropriate because **no production code is written in this plan** (sole `files_modified` is `src/tests/spectronaut-parser.ts`; the webpack/`.d.ts` additions are test-build support, not behavior). No RED gate fail-fast trigger occurred — every new test asserts an already-implemented behavior and passed first run, which is the intended flow for a regression-net plan. The phase-level RED/GREEN gate is satisfied across the wave (Plan 02 `feat(...)` + Plan 03 `test(...)`).

## Issues Encountered
- This freshly-created git worktree had no `node_modules` (same as Plan 01 & Plan 02). Resolved by symlinking the sibling `proteomics` worktree's `node_modules` — verified the sibling is at the **exact same base commit** `ca58537cb988388969e37b5c401883e618ec718f` and this plan changes zero dependencies. The symlink was added to `.git/info/exclude`, never `git add`ed, and is absent from every commit (`git diff --cached --name-only` showed only the intended files).

## Threat Surface Scan
The plan's `<threat_model>` (T-12-10..T-12-14) covers exactly this surface: synthetic committed assets consumed only by `grok test`. The implemented tests mitigate as specified — exact key-set + per-key quantity equality vs the verbatim sidecar (T-12-10/T-12-11 drift), explicit CON__/REV__/>0.01/non-numeric/empty-q assertions + the `sniffIsPrecursor` routing test + per-cell equivalence (T-12-12), the deterministic clean-pass gate (T-12-14). No new security-relevant surface: no network endpoints, no auth, no schema changes, no runtime file writes — only client-side parsing of synthetic fixtures in the test runner. CSV/TSV formula-injection out of scope (assets consumed only by `grok test`/duckdb). No threat flags.

## Next Phase Readiness
- The streaming↔text↔duckdb-golden parity and the D-01 routing branch are now pinned by automated `grok test` coverage; any future change to the aggregator, the fixture, the sidecar, or the routing fails loudly on the next run.
- Every command in 12-VALIDATION.md's Per-Task Verification Map exists and passes: R1 `streams precursor fixture`, R2 `streaming output equals duckdb golden` + `streaming filter parity`, R3 `stream path matches text path`, plus `sniffIsPrecursor routes by header` and `duckdb fallback tooling is committed`.
- Out of scope here (release-gate per SPEC, manual): the 2.6 GB reference-file E2E (R4 — monotonic progress, no Page-Unresponsive, tab switch, full pipeline).
- No blockers.

## Self-Check: PASSED

All claimed files exist on disk: `packages/Proteomics/src/tests/spectronaut-parser.ts`, `packages/Proteomics/webpack.config.js`, `packages/Proteomics/src/global.d.ts`. Both task commits present in git history: `013236a70b`, `203cf570f8`. `npx tsc --noEmit` 0 errors; full Proteomics package suite green on localhost (`Tests passed.`, exit 0).

---
*Phase: 12-spectronaut-input-coverage*
*Completed: 2026-05-15*
