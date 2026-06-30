---
phase: 03-visualization
plan: 01
subsystem: visualization
tags: [scatter-plot, volcano, pca, eigendecomposition, formula-lines]

# Dependency graph
requires:
  - phase: 02-analysis-pipeline
    provides: differential expression results (log2FC, adj.p-value, significant columns)
provides:
  - createVolcanoPlot factory returning configured ScatterPlotViewer with threshold lines
  - ensureNegLog10Column and ensureDirectionColumn helper functions
  - computePCA client-side Jacobi eigendecomposition returning sample-level DataFrame
  - createPcaPlot factory returning ScatterPlotViewer with group coloring and sample labels
affects: [03-02 integration, heatmap viewer, workflow orchestration]

# Tech tracking
tech-stack:
  added: []
  patterns: [factory function returning configured built-in viewer, Jacobi rotation eigendecomposition]

key-files:
  created:
    - packages/Proteomics/src/viewers/volcano.ts
    - packages/Proteomics/src/analysis/pca.ts
    - packages/Proteomics/src/viewers/pca-plot.ts
  modified: []

key-decisions:
  - "Use built-in ScatterPlotViewer factory pattern instead of custom JsViewer for both volcano and PCA"
  - "PCA creates a separate sample-level DataFrame (not added to protein-level DataFrame) since rows represent different entities"
  - "Confidence ellipses use annotation regions API with graceful fallback if API unavailable in installed version"
  - "Direction column with categorical color coding for three-color volcano scheme instead of boolean significant column"

patterns-established:
  - "Viewer factory: function returns DG.ScatterPlotViewer configured with formula lines and color coding"
  - "Derived column pattern: ensureXxxColumn checks existence before computing, returns column name"

requirements-completed: [VIZ-01, VIZ-03]

# Metrics
duration: 6min
completed: 2026-02-28
---

# Phase 3 Plan 1: Volcano & PCA Scatter Plots Summary

**Volcano plot factory with three-color direction scheme and formula-line thresholds, plus client-side Jacobi PCA returning sample-level ScatterPlotViewer with group coloring**

## Performance

- **Duration:** 6 min
- **Started:** 2026-03-01T03:01:06Z
- **Completed:** 2026-03-01T03:07:31Z
- **Tasks:** 2
- **Files modified:** 4 (3 created, 1 deleted)

## Accomplishments
- Volcano plot factory with three formula lines (1 horizontal p-value threshold, 2 vertical FC thresholds)
- Three-color direction column (red=up, blue=down, gray=NS) with categorical color coding
- Client-side PCA via Jacobi eigendecomposition on sample covariance matrix
- PCA creates new sample-level DataFrame with PC column names including variance explained percentages
- 95% confidence ellipse computation per group (best-effort via annotation regions API)
- Deleted old VolcanoViewer JsViewer stub

## Task Commits

Each task was committed atomically:

1. **Task 1: Create volcano plot factory and -log10 column helper** - `0775f12853` (feat)
2. **Task 2: Implement client-side PCA and create PCA plot factory** - `4bb903c385` (feat)

## Files Created/Modified
- `packages/Proteomics/src/viewers/volcano.ts` - Volcano plot factory with ensureNegLog10Column, ensureDirectionColumn, createVolcanoPlot
- `packages/Proteomics/src/analysis/pca.ts` - Client-side PCA with Jacobi eigendecomposition, computePCA returning sample-level DataFrame
- `packages/Proteomics/src/viewers/pca-plot.ts` - PCA plot factory with group coloring, sample labels, and confidence ellipses
- `packages/Proteomics/src/viewers/volcano-viewer.ts` - DELETED (JsViewer stub replaced by scatter plot factory)

## Decisions Made
- Used categorical color coding on direction column (`col.meta.colors.setCategorical`) for three-color volcano scheme
- PCA creates a separate sample-level DataFrame since protein-level and sample-level data have different row counts -- linked selection between volcano and PCA is not applicable (standard in proteomics tools)
- Confidence ellipses use `(sp.meta as any).annotationRegions?.add()` with graceful fallback since the installed datagrok-api type declarations don't include the annotationRegions API yet
- Gene name labels use `displayLabels = 'Auto'` mode which naturally labels isolated extreme points (most significant proteins)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] AreaAnnotationRegion type not in installed datagrok-api**
- **Found during:** Task 2 (PCA plot factory)
- **Issue:** `AreaAnnotationRegion` interface and `annotationRegions` on ViewerMetaHelper exist in source js-api but not in the npm-installed datagrok-api type declarations
- **Fix:** Defined local `AreaAnnotationRegionDef` interface and used `(sp.meta as any).annotationRegions?.add()` with optional chaining for graceful degradation
- **Files modified:** `packages/Proteomics/src/viewers/pca-plot.ts`
- **Verification:** TypeScript compiles clean (only expected package.ts import error remains)
- **Committed in:** `4bb903c385` (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Minor type workaround for annotation regions API. No scope change. Ellipse rendering is best-effort and documented.

## Issues Encountered
- package.ts still imports deleted `volcano-viewer.ts` -- this is expected and will be resolved in Plan 03-02 which updates the package integration

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Volcano and PCA plot factories ready for integration in Plan 03-02
- package.ts import of deleted volcano-viewer.ts must be updated in Plan 03-02
- Heatmap viewer (Plan 03-02) can follow same factory pattern established here

---
*Phase: 03-visualization*
*Completed: 2026-02-28*
