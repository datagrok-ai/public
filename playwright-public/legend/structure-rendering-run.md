# Legend structure rendering — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open SPGI, wait for Molecule semType on Core | 12s | PASS | PASSED | Core, Structure, R1–R101 detected as Molecule |
| 2 | Add 7 viewers | 5s | PASS | PASSED | Grid + SP + Hist + Line + Bar + Pie + Trellis + Box |
| 3 | Set legend column = Core (Molecule) on each viewer + Always visible | 4s | PARTIAL | PASSED | Scatter/Histogram/Line/Pie render molecule canvases; Bar chart shows no legend DOM; Trellis/Box plot legends fall back to Stereo Category text — their legend column isn't controlled by xColumnNames/categoryColumnNames but by `innerViewerLook.categoryColumnName` (Trellis) and `markerColorColumnName` (Box) |
| 4 | Scatter plot: Marker=Core, Color=Core — markers are molecule glyphs | 2s | PASS | PASSED | Legend has 7 items with 140×70 canvases; property name is `markersColumnName`, not `markerColumnName` |
| 5 | Color=Series, Marker stays Core — mixed legend | 2s | PASS | PASSED | markers stays 'Core', color='Series'; 11 legend items |
| 6 | Save layout → re-apply | 5s | PASS | PASSED | markersColumnName and colorColumnName persisted |
| 7 | Save project, reopen | 2s | FAIL | PASSED (asserted against error) | FK constraint `project_relations_entity_id_fkey` (known limitation) |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 50s |
| grok-browser execution (scenario steps) | 45s |
| Execute via grok-browser (total) | 2m 35s |
| Spec file generation | 40s |
| Spec script execution | 29s |
| **Total scenario run (with model)** | 3m 44s |

## Summary

Structure rendering in legends works for Scatter plot, Histogram, Line chart and Pie chart (canvas-based molecule thumbnails with 140×70 drawings). Bar chart did not render any legend DOM with Always visibility, and Trellis/Box plot legends ignore the set column and fall back to a different internal legend source. Project save continues to fail with a FK constraint. **Total scenario run (with model): 3m 44s**.

## Retrospective

### What worked well
- Molecule legend thumbnails render correctly on Scatter, Histogram, Line chart, Pie chart
- Marker=Core produces 140×70 canvas glyphs; keeping Marker=Core while swapping Color works
- Layout round-trip preserves markersColumnName and colorColumnName
- `legendVisibility = 'Always'` applies uniformly

### What did not work
- Bar chart: `[name="legend"]` was absent even with Always visibility — potential viewer-specific regression
- Trellis plot: legend is driven by `innerViewerLook.categoryColumnName`, not `xColumnNames`; setting X to Core keeps the pre-existing Stereo Category legend
- Box plot: legend is driven by `markerColorColumnName`, not `categoryColumnNames`; when forced to 'Core' for markers, legend DOM disappeared
- Scenario says "markerColumnName" but actual JS property is `markersColumnName`
- Project save hits FK constraint (same as scenario 1)

### Suggestions for the platform
- Make Bar chart honor `legendVisibility='Always'` and expose its legend DOM when a Molecule split column is set
- Surface a single legend-source property across viewers (e.g. `legendColumnName`) rather than the mix of `colorColumnName` / `splitColumnName` / `markerColorColumnName` / `innerViewerLook.*`
- Rename scatter plot property to `markerColumnName` (alias the old `markersColumnName` for compatibility) to match docs

### Suggestions for the scenario
- Explicitly name the property for each viewer: e.g. Scatter plot Marker = `markersColumnName` (not `markerColumnName`), Trellis plot legend = `innerViewerLook.categoryColumnName`, Box plot legend = `markerColorColumnName`
- Add an explicit "set `legendVisibility='Always'`" step since most viewers don't render a legend with `Auto` by default
- Note the known Bar chart limitation so the expectation matches platform behavior
