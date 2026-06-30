---
phase: 02-analysis-pipeline
plan: 02
subsystem: analysis
tags: [welch-t-test, fdr, benjamini-hochberg, differential-expression, statistics]

requires:
  - phase: 01-data-import-and-foundation
    provides: MaxQuant parser with intensity columns and semantic types
provides:
  - runDifferentialExpression() with Welch's t-test and BH FDR correction
  - showDEDialog() with prerequisite validation and threshold inputs
  - log2FC, p-value, adj.p-value, and significant result columns
affects: [03-visualization]

tech-stack:
  added: [@datagrok-libraries/statistics (tTest, fdrcorrection)]
  patterns: [client-side statistical testing, BH FDR correction, null-safe protein analysis]

key-files:
  created: []
  modified: [packages/Proteomics/src/analysis/differential-expression.ts]

key-decisions:
  - "Client-side Welch's t-test instead of R limma (per user decision, limma deferred to Phase 4)"
  - "Proteins with <2 replicates per group get null p-values rather than throwing errors"

patterns-established:
  - "Per-protein statistical testing with null-safe value extraction across sample columns"
  - "FDR correction applied only to testable proteins, untestable get FLOAT_NULL"
  - "State tag proteomics.de_complete prevents double-application"

requirements-completed: [ANLY-03, ANLY-05]

duration: 2min
completed: 2026-02-28
---

# Phase 02 Plan 02: Differential Expression Summary

**Client-side Welch's t-test with Benjamini-Hochberg FDR correction producing log2FC, p-value, adj.p-value, and significance columns**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-28T21:21:57Z
- **Completed:** 2026-02-28T21:24:00Z
- **Tasks:** 1
- **Files modified:** 1

## Accomplishments
- Replaced limma R stub with client-side Welch's t-test per-protein loop
- BH FDR correction via fdrcorrection from @datagrok-libraries/statistics
- Four result columns: log2FC, p-value, adj.p-value (float), significant (boolean)
- Semantic types assigned: Proteomics-Log2FC and Proteomics-PValue
- Graceful null handling for proteins with insufficient replicates
- Dialog validates group annotation exists and prevents re-running DE

## Task Commits

Each task was committed atomically:

1. **Task 1: Implement Welch's t-test DE with FDR correction and significance column** - `eaa5f0a653` (feat)

## Files Created/Modified
- `packages/Proteomics/src/analysis/differential-expression.ts` - Client-side DE with Welch's t-test, BH FDR, significance column, and dialog

## Decisions Made
- Used client-side Welch's t-test per user decision (R limma deferred to Phase 4)
- Proteins with <2 replicates per group get DG.FLOAT_NULL instead of throwing errors
- Significance requires both |log2FC| >= threshold AND adj.p-value <= threshold

## Deviations from Plan

None - plan executed exactly as written. The experiment-setup.ts dependency already existed in the scaffold (no stub creation needed).

## Issues Encountered
- Minor TypeScript strict null check on `col.get()` return type (nullable float) -- resolved with type assertion after null guard check

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- DE computation ready for volcano plot visualization in Phase 3
- Result columns have correct semantic types for viewer detection
- proteomics.de_complete tag enables downstream pipeline awareness

---
*Phase: 02-analysis-pipeline*
*Completed: 2026-02-28*

## Self-Check: PASSED
