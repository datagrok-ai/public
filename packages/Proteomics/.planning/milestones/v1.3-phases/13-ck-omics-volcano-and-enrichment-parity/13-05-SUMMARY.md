---
phase: 13-ck-omics-volcano-and-enrichment-parity
plan: 05
subsystem: api
tags: [differential-expression, spectronaut, contrast-direction, tdd]

requires:
  - phase: 13-ck-omics-volcano-and-enrichment-parity
    provides: 13-WAVE0-FINDINGS.md A3 (D-09 fix site confirmed)
provides:
  - getDefaultComparison — DE dialog defaults to the declared/intended contrast
affects: [13-06]

tech-stack:
  added: []
  patterns:
    - "Direction-only default change; OK-handler index mapping untouched"

key-files:
  created: []
  modified:
    - src/analysis/differential-expression.ts
    - src/tests/analysis.ts

key-decisions:
  - "Default = `${g1} vs ${g2}` (numerator = group1 = declared Spectronaut Numerator), not alphabetical pairs[0]"
  - "Extracted getDefaultComparison so the default is unit-testable without driving the UI dialog"
  - "spectronaut-parser.ts intentionally NOT modified (fix pinned to the dialog default per A3)"

patterns-established:
  - "Workaround demoted: manual dropdown override stays, correct default no longer requires it"

requirements-completed: [R3]

duration: 23min
completed: 2026-05-17
---

# Phase 13 Plan 05: Report-DE Direction Default (D-09) Summary

**The DE dialog now defaults the Comparison to the declared/intended contrast (numerator = group1 = Spectronaut's declared Numerator) instead of the alphabetical pairs[0]; direction-only, dropdown override preserved, hit-list magnitudes unchanged.**

## Performance

- **Duration:** ~23 min (inline sequential; one-line fix + helper + 2 tests + build/test)
- **Started:** 2026-05-17T21:59Z
- **Completed:** 2026-05-17T22:22Z
- **Tasks:** 1 (TDD RED → GREEN)
- **Files modified:** 2

## Accomplishments
- `getDefaultComparison(g1, g2)` → `${g1} vs ${g2}` (= `pairs[1]`); `showDEDialog`'s `comparisonInput` default sources from it instead of `pairs[0]`. The OK-handler's `value === pairs[1]` logic resolves this to numerator = group1 = the parser's first condition = Spectronaut's declared Numerator (13-WAVE0-FINDINGS A3), unparking the documented mirror defect.
- Direction-only: `items` (both orientations), the OK-handler index mapping, `updateHint`, DE statistics/thresholds/method paths all byte-unchanged. `src/parsers/spectronaut-parser.ts` intentionally untouched (verified) — fix contained to the two declared files.
- 2 new tests in the `Differential Expression` category: default = declared contrast; direction-only invariant (reversal flips sign, |log2FC| identical). All 8 category tests green on localhost, exit 0.

## Task Commits

1. **Task 1 RED: default-orientation tests** - `47358a9b90` (test)
2. **Task 1 GREEN: declared-contrast default** - `45fd185fc0` (feat)

**Plan metadata:** this commit (docs: complete plan)

## Files Created/Modified
- `src/analysis/differential-expression.ts` - getDefaultComparison + default value sourced from it
- `src/tests/analysis.ts` - 2 D-09 regression tests

## Decisions Made
See key-decisions frontmatter. Extracting `getDefaultComparison` (rather than inlining `pairs[1]`) makes the default genuinely unit-testable without instantiating the UI dialog — the change is still "only the comparisonInput default value", now sourced from a one-line pure helper in the same declared file.

## Deviations from Plan

None - plan executed exactly as written. (The plan's `<verify>` names `--category Analysis`; the analysis.ts suite has no `Analysis` category — the change lives in `Differential Expression`, which was run green. Category-name imprecision in the plan, not a behavior change.)

## Issues Encountered
None.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Wave 2 complete (13-03, 13-04, 13-05). Wave 3 (13-06 volcano) unblocked.
- The end-to-end "unaware analyst gets correct volcano sign by default" is a HUMAN-UAT item (requires importing a real Spectronaut report and running DE through the dialog).

---
*Phase: 13-ck-omics-volcano-and-enrichment-parity*
*Completed: 2026-05-17*
