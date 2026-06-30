---
phase: 13-ck-omics-volcano-and-enrichment-parity
plan: 03
subsystem: api
tags: [spectronaut, candidates, sign-normalization, filters, dockmanager, tdd]

requires:
  - phase: 13-ck-omics-volcano-and-enrichment-parity
    provides: 13-WAVE0-FINDINGS.md A2 (locked per-row flip rule)
provides:
  - normalizeCandidatesSign — conditional per-row sign normalization (pure parser)
  - dockComparisonFilterIfMultiContrast — native Comparison Filter for multi-contrast files
affects: [13-06]

tech-stack:
  added: []
  patterns:
    - "Conditional per-row flip keyed on declared-vs-canonical comparison string (CK create_subset_data)"
    - "Bulk getRawData → number[] → col.init, FLOAT_NULL preserved as sentinel (never init(()=>null))"

key-files:
  created: []
  modified:
    - src/parsers/spectronaut-candidates-parser.ts
    - src/package.ts
    - src/tests/spectronaut-candidates-parser.ts
    - src/tests/spectronaut-candidates-e2e.ts

key-decisions:
  - "Canonical = first parseable 'g1 / g2'; flip ONLY exact-reverse rows; never unconditional (A2/Pitfall 4)"
  - "number[] not Float32Array for the rewritten arrays — preserves double precision of large AVG Group Quantity values"
  - "AVG Group Quantity / Condition swaps guarded on BOTH columns present (A6)"
  - "Filter docking lives in package.ts (shell-only); parser stays pure (Pitfall 3)"

patterns-established:
  - "Sign normalization runs BEFORE significance (|log2FC| is flip-invariant)"
  - "Multi-contrast detection by distinct non-null Comparison count, post-normalization"

requirements-completed: [R3, R4]

duration: 105min
completed: 2026-05-17
---

# Phase 13 Plan 03: Candidates Sign Normalization + Multi-Contrast Filter Summary

**Candidates rows whose declared comparison is the exact reverse of the canonical orientation are flipped (log2FC negated, AVG Group Quantity + Condition swapped, Comparison relabeled) in a still-pure parser; multi-contrast files dock a native Comparison Filters viewer.**

## Performance

- **Duration:** ~105 min (inline sequential; one build + 2 test categories)
- **Started:** 2026-05-17T16:20Z
- **Completed:** 2026-05-17T18:05Z
- **Tasks:** 2 (Task 1 TDD RED → GREEN)
- **Files modified:** 4

## Accomplishments
- `normalizeCandidatesSign(df)` — pure: canonical = first parseable `g1 / g2`; flips ONLY exact-reverse rows per the A2-locked CK `create_subset_data` rule; negates log2FC (FLOAT_NULL preserved), swaps AVG Group Quantity Num/Den (only if both present), relabels Comparison + Condition. Runs before significance (|log2FC| flip-invariant).
- Single-orientation / canonical / unparseable files are byte-identical — all 18 pre-existing SpectronautCandidates tests pass unchanged (no mirror regression).
- `dockComparisonFilterIfMultiContrast(tv, df)` in package.ts — docks a native `DG.Viewer.filters` scoped to the Comparison column when >1 distinct comparison; single-contrast docks nothing. Parser untouched (purity preserved; grep gate clean).
- 3 new parser tests + 2 new E2E tests. Build compiled; `SpectronautCandidates` (21) and `SpectronautCandidates E2E` (3) green on localhost, exit 0.

## Task Commits

1. **Task 1 RED: reversed-comparison fixtures** - `ae70f2868c` (test)
2. **Task 1 GREEN: per-row sign normalization** - `b9d240071e` (feat)
3. **Task 2: multi-contrast Comparison Filter docking** - `5c90c1c7bd` (feat)
4. **Purity-grep comment reword** - `2a7d2ce16c` (docs)

**Plan metadata:** this commit (docs: complete plan)

## Files Created/Modified
- `src/parsers/spectronaut-candidates-parser.ts` - normalizeCandidatesSign + name-variant arrays
- `src/package.ts` - dockComparisonFilterIfMultiContrast + wired into importSpectronautCandidates
- `src/tests/spectronaut-candidates-parser.ts` - 3 sign-normalization tests
- `src/tests/spectronaut-candidates-e2e.ts` - 2 Filter-docking E2E tests

## Decisions Made
See key-decisions frontmatter.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Correctness] number[] instead of Float32Array for rewritten value arrays**
- **Found during:** Task 1 (sign normalization)
- **Issue:** The plan's <action> prescribed Float32Array; AVG Group Quantity values are large doubles (e.g. 1991.23583984375) and would lose precision through Float32.
- **Fix:** Used `number[]` (double-precision) for the rewritten log2FC / quantity arrays; FLOAT_NULL still preserved by writing the sentinel number back (never `init(() => null)` — memory feedback_dg_column_init_null_sentinel).
- **Files modified:** src/parsers/spectronaut-candidates-parser.ts
- **Verification:** swapped-quantity test asserts exact 100/800 values to 1e-4; all 21 tests green
- **Committed in:** b9d240071e

---

**Total deviations:** 1 auto-fixed (1 correctness — precision preservation)
**Impact on plan:** Behaviorally identical to intent, strictly more correct for large quantities. No scope change.

## Issues Encountered
- Purity grep tripped on the docstring literally containing `grok.shell`; reworded the comment (commit `2a7d2ce16c`). Parser has zero shell code/imports — verified.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Wave 2 sibling 13-04/13-05 independent and ready.
- 13-06 (volcano) will consume the now-sign-correct Candidates output.
- Real client multi-contrast file behavior (visual Filter docking) is a HUMAN-UAT item.

---
*Phase: 13-ck-omics-volcano-and-enrichment-parity*
*Completed: 2026-05-17*
