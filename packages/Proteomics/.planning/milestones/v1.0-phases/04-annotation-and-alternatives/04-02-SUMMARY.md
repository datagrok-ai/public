---
phase: 04-annotation-and-alternatives
plan: 02
subsystem: analysis
tags: [deqms, limma, differential-expression, r-script, proteomics]

requires:
  - phase: 02-analysis-pipeline
    provides: "limma DE pipeline, runLimmaDE function, buildExpressionDf helper"
provides:
  - "DEqMS R script (deqms_de.R) with peptide-count-weighted variance estimation"
  - "Method selector (limma/DEqMS) in DE dialog"
  - "Peptide count column picker with auto-detection"
  - "Three-level fallback chain: DEqMS -> limma -> client-side t-test"
affects: [visualization, analysis-pipeline]

tech-stack:
  added: [DEqMS]
  patterns: [server-side-r-fallback-chain, conditional-dialog-inputs]

key-files:
  created:
    - packages/Proteomics/scripts/deqms_de.R
  modified:
    - packages/Proteomics/src/analysis/differential-expression.ts

key-decisions:
  - "DEqMS R script uses limma as prerequisite base then applies spectraCounteBayes"
  - "Peptide counts passed as single-column DataFrame (Datagrok R input constraint)"
  - "Missing peptide counts default to 1 (conservative: no count weighting advantage)"
  - "ui.input.column value takes undefined not null for optional default"

patterns-established:
  - "Three-level R fallback: specialized -> general -> client-side"
  - "Conditional dialog input visibility via onChanged + style.display toggle"
  - "proteomics.de_method tag for tracking which DE method was used"

requirements-completed: [ANLY-04]

duration: 3min
completed: 2026-03-02
---

# Phase 04 Plan 02: DEqMS Alternative DE Method Summary

**DEqMS peptide-count-weighted DE with method selector dialog and three-level fallback chain**

## Performance

- **Duration:** 3 min
- **Started:** 2026-03-02T03:27:27Z
- **Completed:** 2026-03-02T03:31:03Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Created DEqMS R script with spectraCounteBayes variance estimation and limma fallback
- Enhanced DE dialog with method selector (limma/DEqMS) and conditional peptide count column picker
- Implemented three-level fallback: DEqMS -> limma -> client-side t-test
- Output columns identical to limma contract so existing viewers work without changes

## Task Commits

Each task was committed atomically:

1. **Task 1: Create DEqMS R script** - `e02700a10a` (feat)
2. **Task 2: Add method selector and peptide count picker** - `77aeb1a529` (feat)

## Files Created/Modified
- `packages/Proteomics/scripts/deqms_de.R` - DEqMS R script with peptide-count-weighted DE and limma fallback
- `packages/Proteomics/src/analysis/differential-expression.ts` - Method selector, peptide count picker, runDeqmsDE function, three-level fallback

## Decisions Made
- DEqMS R script uses limma as prerequisite base then applies spectraCounteBayes -- keeps code DRY
- Peptide counts passed as single-column DataFrame because Datagrok R inputs support dataframe but not raw column vectors
- Missing peptide counts default to 1 (conservative approach: no count weighting advantage for unknown proteins)
- ui.input.column value parameter takes undefined not null for optional default column

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed null vs undefined type mismatch for column input default**
- **Found during:** Task 2 (DE dialog modification)
- **Issue:** Plan used `?? null` for default peptide column, but `ui.input.column` value expects `undefined` not `null`
- **Fix:** Changed `?? null` to `?? undefined`
- **Files modified:** packages/Proteomics/src/analysis/differential-expression.ts
- **Verification:** TypeScript compiles cleanly with `npx tsc --noEmit`
- **Committed in:** 77aeb1a529 (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Minor type fix for API compatibility. No scope creep.

## Issues Encountered
- Linter reverted file changes during first commit attempt (git lock file conflict). Resolved by removing stale lock file and re-applying changes via full file write.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- DEqMS method available as alternative to limma for peptide-count-aware DE
- All existing viewers (volcano, heatmap, PCA) work with DEqMS output without changes
- R environment on server needs DEqMS Bioconductor package for full functionality (graceful fallback to limma if absent)

---
*Phase: 04-annotation-and-alternatives*
*Completed: 2026-03-02*
