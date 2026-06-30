---
phase: 11-dialog-expansion-and-ux-polish
plan: 02
subsystem: ui
tags: [dialog, imputation, differential-expression, method-selector, comparison-direction]

requires:
  - phase: 10-spectronaut-parser-and-core-algorithms
    provides: "imputeKnn, imputeZero, imputeMean, imputeMedian functions"
provides:
  - "Multi-method imputation dialog with conditional parameter visibility"
  - "Valid-values protein filter with live count preview"
  - "DE comparison direction picker with hint text"
  - "t-test as explicit DE method option"
affects: [11-03, normalization-dialog, viewer-titles]

tech-stack:
  added: []
  patterns: [conditional-visibility-containers, live-count-preview, comparison-direction-parsing]

key-files:
  created: []
  modified:
    - packages/Proteomics/src/analysis/imputation.ts
    - packages/Proteomics/src/analysis/differential-expression.ts

key-decisions:
  - "Valid-values filter uses ANY-group threshold (protein passes if any group meets minimum)"
  - "Comparison direction parsed from string format 'GroupA vs GroupB' split on ' vs '"
  - "Default comparison is group2 vs group1, matching existing FC direction convention"

patterns-established:
  - "Container-div pattern: wrap related inputs in ui.div, toggle display on method change"
  - "Comparison direction parsing: split selected string on ' vs ' to determine numerator/denominator"

requirements-completed: [IMP-03, IMP-04, DE-01, DE-02]

duration: 2min
completed: 2026-03-07
---

# Phase 11 Plan 02: Dialog Expansion Summary

**Multi-method imputation dialog with 5 methods, conditional params, valid-values filter; DE dialog with comparison direction picker, t-test method, and FC hint text**

## Performance

- **Duration:** 2 min
- **Started:** 2026-03-07T21:37:53Z
- **Completed:** 2026-03-07T21:40:02Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Imputation dialog now offers MinProb, kNN, Zero, Mean, Median with conditional parameter visibility per method
- Valid-values filter counts proteins per group and removes those below threshold before imputation
- DE dialog has comparison direction dropdown with auto-generated directional pairs and dynamic FC hint text
- t-test available as explicit third DE method (previously only as fallback)

## Task Commits

Each task was committed atomically:

1. **Task 1: Expand imputation dialog** - `120cc3de00` (feat)
2. **Task 2: Expand DE dialog** - `5bf03155c3` (feat)

## Files Created/Modified
- `packages/Proteomics/src/analysis/imputation.ts` - Added method selector, conditional param containers, valid-values filter with live count
- `packages/Proteomics/src/analysis/differential-expression.ts` - Added comparison direction picker, hint div, t-test method, numerator/denominator mapping

## Decisions Made
- Valid-values filter: protein passes if ANY group meets threshold (not ALL groups) -- prevents over-aggressive filtering
- Comparison direction defaults to group2 vs group1, preserving existing FC sign convention
- t-test branch calls runDifferentialExpression directly without R server dependency

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed removeWhere signature mismatch**
- **Found during:** Task 1 (imputation dialog)
- **Issue:** `df.rows.removeWhere(BitSet)` expects a `RowPredicate`, not a `BitSet` directly
- **Fix:** Changed to `df.rows.removeWhere((row) => removeSet.get(row.idx))`
- **Files modified:** packages/Proteomics/src/analysis/imputation.ts
- **Verification:** TypeScript compiles without errors
- **Committed in:** 120cc3de00 (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Minor API mismatch fix, no scope change.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Imputation and DE dialogs fully expanded with all planned controls
- Ready for normalization dialog expansion (Plan 01) and viewer title/naming polish (Plan 03)

---
*Phase: 11-dialog-expansion-and-ux-polish*
*Completed: 2026-03-07*
