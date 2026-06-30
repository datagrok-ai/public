---
phase: 07-qc-dashboard
plan: 03
subsystem: ui
tags: [datagrok, qc-dashboard, dg-viewer, dock-manager, dataframe]

requires:
  - phase: 07-qc-dashboard (plan 01)
    provides: QC computation functions and dashboard layout
  - phase: 07-qc-dashboard (plan 02)
    provides: Menu wiring for QC dashboard
provides:
  - Fixed QC dashboard with all 7 viewers rendering correctly using DG.Viewer factory methods
affects: []

tech-stack:
  added: []
  patterns:
    - "DG.Viewer factory methods (barChart, boxPlot, grid) for auxiliary DataFrames instead of tv.addViewer"
    - "Avoid grok.shell.addTable() for DataFrames used only as viewer data sources"

key-files:
  created: []
  modified:
    - packages/Proteomics/src/viewers/qc-dashboard.ts

key-decisions:
  - "Used DG.Viewer.grid() instead of missDf.plot.grid() for proper Viewer instance that docks correctly"
  - "Used 'as any' cast for barChart/boxPlot options to work around strict Partial<ISettings> typings"

patterns-established:
  - "Auxiliary DataFrame pattern: create DF locally, use DG.Viewer.factory(df, opts), dock directly -- never addTable"

requirements-completed: [QC-01, QC-03, QC-05]

duration: 1min
completed: 2026-03-07
---

# Phase 7 Plan 3: QC Dashboard UAT Gap Closure Summary

**Fixed 4 UAT failures by replacing grok.shell.addTable + tv.addViewer with DG.Viewer factory methods for auxiliary DataFrames**

## Performance

- **Duration:** 1 min
- **Started:** 2026-03-07T00:25:12Z
- **Completed:** 2026-03-07T00:26:16Z
- **Tasks:** 2
- **Files modified:** 1

## Accomplishments
- Eliminated circular JSON serialization error caused by grok.shell.addTable() on auxiliary DataFrames
- Fixed missingness heatmap using DG.Viewer.grid(missDf) instead of missDf.plot.grid()
- Fixed missing values bar chart using DG.Viewer.barChart(barDf) bound to correct DataFrame
- Fixed intensity box plot using DG.Viewer.boxPlot(longDf) bound to correct DataFrame
- Package builds and bundles successfully with webpack

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix auxiliary DataFrame viewer creation in QC dashboard** - `9af4184248` (fix)
2. **Task 2: Build package to verify webpack bundle** - no commit (verification-only task, no file changes)

## Files Created/Modified
- `packages/Proteomics/src/viewers/qc-dashboard.ts` - Replaced broken addTable+addViewer pattern with DG.Viewer factory methods for 3 auxiliary DataFrame viewers

## Decisions Made
- Used DG.Viewer.grid() instead of df.plot.grid() to get a proper Viewer instance compatible with DockManager
- Used `as any` type cast on factory method options because the strict Partial<ISettings> interfaces don't include all viewer-specific properties like splitColumnName/valueColumnName
- Removed allowColSelection=false on missGrid (cosmetic, not blocking) and dropped isHeatmap/isGrid props that are not available on Grid returned by DG.Viewer.grid()

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- QC Dashboard is fully functional with all 7 viewers docked in tiled layout
- Ready for Phase 8 (Gene ID Mapping and Enrichment Analysis)

---
*Phase: 07-qc-dashboard*
*Completed: 2026-03-07*
