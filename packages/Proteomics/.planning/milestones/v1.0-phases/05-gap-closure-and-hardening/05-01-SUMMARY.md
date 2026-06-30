---
phase: 05-gap-closure-and-hardening
plan: 01
subsystem: visualization
tags: [heatmap, dendrogram, clustering, detectors, differential-expression, r-scripts]

requires:
  - phase: 03-visualization
    provides: Expression heatmap grid with z-score normalization
  - phase: 04-annotation-and-alternatives
    provides: R-based DE scripts (limma, DEqMS) and three-level fallback

provides:
  - Corrected PascalCase R function call names (LimmaDE, DeqmsDE)
  - Extended detector regex for log2-wrapped intensity columns
  - Filter-isolated heatmap via DataFrame clone
  - Dendrogram tree visual alongside heatmap grid via IDendrogramService

affects: []

tech-stack:
  added: []
  patterns:
    - "DataFrame clone for viewer filter isolation (df.clone(filter))"
    - "TreeHelper + DendrogramService for hierarchical clustering with visual tree injection"

key-files:
  created: []
  modified:
    - packages/Proteomics/src/analysis/differential-expression.ts
    - packages/Proteomics/detectors.js
    - packages/Proteomics/src/viewers/heatmap.ts

key-decisions:
  - "DistanceMetric.Euclidean enum instead of string literal for type safety"
  - "Z-score statistics computed from all original rows, applied only to cloned top-N rows"

patterns-established:
  - "Clone-based filter isolation: clone DataFrame with BitSet filter before creating Grid, preventing filter side-effects on other viewers"
  - "Dendrogram service pattern: getTreeHelper + getDendrogramService with graceful try/catch fallback"

requirements-completed: [VIZ-02]

duration: 3min
completed: 2026-03-03
---

# Phase 5 Plan 1: Gap Closure and Hardening Summary

**Fixed R script call names, log2 detector regex, heatmap filter isolation via DataFrame clone, and dendrogram tree injection via IDendrogramService**

## Performance

- **Duration:** 3 min
- **Started:** 2026-03-03T14:15:18Z
- **Completed:** 2026-03-03T14:18:13Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- Corrected PascalCase mismatch in R function calls (limmaDE -> LimmaDE, deqmsDE -> DeqmsDE) so server-side DE actually executes
- Extended detectIntensity to recognize log2-wrapped column names (e.g. 'log2(LFQ intensity Sample1)') for semType persistence across save/reload
- Isolated heatmap filter from shared DataFrame by cloning with top-N BitSet filter -- volcano plot retains all proteins
- Replaced hand-rolled distance matrix and _heatmap_sort_order hack with TreeHelper.calcDistanceMatrix + hierarchicalClusteringByDistance + dendrogramService.injectTreeForGrid

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix R script call names, detector regex, and heatmap filter isolation** - `dc7f89f0ed` (fix)
2. **Task 2: Render dendrogram tree visual alongside heatmap grid** - `35c8d1f238` (feat)

## Files Created/Modified
- `packages/Proteomics/src/analysis/differential-expression.ts` - PascalCase R function call names + improved fallback warning messages
- `packages/Proteomics/detectors.js` - Extended detectIntensity with log2-wrapped column name matching
- `packages/Proteomics/src/viewers/heatmap.ts` - Cloned DataFrame for filter isolation + TreeHelper/DendrogramService clustering with visual tree injection

## Decisions Made
- Used DistanceMetric.Euclidean enum value instead of string literal 'euclidean' -- TypeScript type system requires the enum for calcDistanceMatrix parameter
- Z-score statistics (mean, std) computed across ALL rows of original DataFrame for statistical correctness, then z-score values computed on cloned top-N rows only

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] DistanceMetric enum required instead of string literal**
- **Found during:** Task 2 (Dendrogram tree visual)
- **Issue:** TypeScript compilation error -- calcDistanceMatrix expects DistanceMetric enum, not string 'euclidean'
- **Fix:** Imported DistanceMetric from @datagrok-libraries/bio/src/trees/consts and used DistanceMetric.Euclidean
- **Files modified:** packages/Proteomics/src/viewers/heatmap.ts
- **Verification:** npx tsc --noEmit passes cleanly
- **Committed in:** 35c8d1f238 (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Type-level fix required for compilation. No scope change.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All v1.0 milestone audit risks (RISK-01, RISK-02, RISK-03) resolved
- Heatmap viewer fully hardened with dendrogram visual, filter isolation, and proper R script integration
- Package ready for final verification and release

---
*Phase: 05-gap-closure-and-hardening*
*Completed: 2026-03-03*
