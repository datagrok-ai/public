---
phase: 07-qc-dashboard
plan: 01
subsystem: visualization
tags: [scatter-plot, box-plot, corr-plot, dock-manager, ma-plot, cv, missingness]

# Dependency graph
requires:
  - phase: 01-data-import-and-foundation
    provides: "intensity columns with SEMTYPE.INTENSITY and log2 prefix"
  - phase: 02-analysis-pipeline
    provides: "group annotations via getGroups(), DE results for conditional MA coloring"
provides:
  - "QC computation functions: computeMA, computeCV, computeLoessTrend, createMissingnessMatrix, unpivotIntensities, computeMissingBarData, getIntensityColumns"
  - "Dashboard orchestration function: openQcDashboard"
affects: [07-02, menu-wiring]

# Tech tracking
tech-stack:
  added: []
  patterns: [multi-viewer-dock-layout, computed-columns-on-main-df, separate-df-for-different-granularity, moving-average-trend-line]

key-files:
  created:
    - packages/Proteomics/src/viewers/qc-computations.ts
    - packages/Proteomics/src/viewers/qc-dashboard.ts
  modified: []

key-decisions:
  - "Moving-average approximation for loess trend line -- sliding window on sorted A values, sufficient for QC bias detection"
  - "Separate DataFrames for missingness heatmap, bar chart, and box plot -- cross-selection with protein-level viewers not possible (accepted per CONTEXT.md)"
  - "MA trend shown as separate small scatter plot docked below MA plot rather than overlay"
  - "CV plot defaults to group1 replicates"

patterns-established:
  - "DockManager layout: capture DockNode return values from dock() for nested docking"
  - "Column cleanup on re-run: remove and recreate computed columns to prevent accumulation"
  - "ensureFreshFloat helper: removes column if exists then adds fresh float column"

requirements-completed: [QC-01, QC-02, QC-03, QC-04, QC-05, QC-06]

# Metrics
duration: 3min
completed: 2026-03-06
---

# Phase 7 Plan 1: QC Dashboard Summary

**Pure computation utilities (MA, CV, loess trend, missingness, unpivot) plus dashboard orchestration docking 7 viewers (MA, trend, CV, correlation, missing heatmap, missing bar, box plot) in tiled layout**

## Performance

- **Duration:** 3 min
- **Started:** 2026-03-06T21:49:25Z
- **Completed:** 2026-03-06T21:52:25Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- 7 exported computation functions covering all QC metrics (MA, CV, loess trend, missingness matrix, unpivot, missing bar data, intensity column detection)
- Dashboard orchestration function creating and docking 7 viewers in a tiled layout
- Conditional MA plot coloring by direction when DE results are available
- Null/edge case guards throughout (DG.FLOAT_NULL for insufficient data)

## Task Commits

Each task was committed atomically:

1. **Task 1: Create QC computation utilities** - `a0d61eeb79` (feat)
2. **Task 2: Create QC dashboard orchestration with viewer layout** - `37adbdc05c` (feat)

## Files Created/Modified
- `packages/Proteomics/src/viewers/qc-computations.ts` - Pure computation functions for MA, CV, loess trend, missingness, unpivot, missing bar data
- `packages/Proteomics/src/viewers/qc-dashboard.ts` - Dashboard orchestration with viewer creation and dock layout management

## Decisions Made
- Used moving-average sliding window as loess approximation -- simple, efficient, sufficient for showing intensity-dependent bias in QC overview
- MA trend line displayed as a separate small scatter plot docked below the MA plot since Datagrok scatter plot does not support native overlay of a second series
- CV plot defaults to group1 -- both groups' CV columns are computed and available for later switching
- missGrid, barDf, and longDf use separate DataFrames registered via grok.shell.addTable() -- cross-selection with protein-level viewers is not possible (accepted limitation)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Fixed DockManager dock() refNode type mismatch**
- **Found during:** Task 2 (QC dashboard orchestration)
- **Issue:** Plan suggested passing viewer objects as refNode parameter to dock(), but dock() requires DockNode (returned by previous dock() calls)
- **Fix:** Captured DockNode return values from each dock() call and used them as refNode for subsequent docking operations
- **Files modified:** packages/Proteomics/src/viewers/qc-dashboard.ts
- **Verification:** TypeScript compiles without errors
- **Committed in:** 37adbdc05c (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Type-safety fix required for correct DockManager API usage. No scope creep.

## Issues Encountered
None beyond the dock() type mismatch resolved above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- openQcDashboard exported and ready to be wired into package.ts menu (Plan 07-02)
- All computation functions exported for unit testing (Plan 07-02)
- Computed columns pattern established for QC metrics on main DataFrame

---
*Phase: 07-qc-dashboard*
*Completed: 2026-03-06*
