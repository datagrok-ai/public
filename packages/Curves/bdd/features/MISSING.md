# Curves features — what is missing

## Resolved

- **Multi curve viewer readings** (2026-10-18): `getWidgetStatus()` reports `curves`, `legend shown`, `legend corner`, `legend rows`, `legend row N`, `legend shortened rows`, `legend hidden` and a `legend row N` / `shortened legend row N` / `legend more` area per legend row, with `isRenderPending` and `onRendered`; `multi-curve-viewer` and the `curves-legend` guide claim them.
- **The axes of the Multi curve viewer** (GROK-20661): `getWidgetStatus()` lays the chart out as the render does
  (`FitChartCellRenderer.layOutChart`) and reports `x range`, `y range` and `points outside the axes`; the feature
  claims every point inside them with all five columns, the package test a cell's declared range not narrowing
  the overlay.

## Still open

Nothing planned.

## Won't do

- **Curve cell readings in the grid and in Context Panel > Chart**: what a curve cell draws is the renderer's (ruling of 2026-10-06); the package tests it (legend-tests, panel-renderer-tests, converter-tests).
- **Points as hit areas in a cell** (outlier toggling, curves-statistics.md scenario 4): renderer drawing again, covered by the package tests, not a feature (ruling of 2026-10-06).
