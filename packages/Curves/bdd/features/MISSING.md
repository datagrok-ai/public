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

## Data to Curves

Feature: `data/data-to-curves.feature` (30 Oct 2026, local stand).

### Questions (Data to Curves)

- **No help page** covers the Data to Curves dialog; the feature claims what `data-to-curves.ts` does.
- **The column pickers report no value.** `should have value` on `Concentration input` reads nothing (a Dart
  column input publishes no value), and `contain text "concentration"` would match the label itself; the
  feature names the column `concentration_uM` so that the text claim cannot pass on the label.

### Findings 30 Oct 2026 (Data to Curves)

- **Data > Curves > Data to Curves is missing from the table view's Data menu.** Steps: open any table; open
  Data > Curves. Measured: the group holds Calculate MSR... only. Expected: Data to Curves, which the function
  declares (`'top-menu': 'Data | Curves | Data to Curves'`, Curves `src/package.ts`). Cause: the table view's
  menu adds a package function only when it takes a dataframe input (core `xamgle/lib/src/views/table_view.dart`
  `fillMenu`), and `dataToCurvesTopMenu` takes none (it declares `requiresTableView`). PowerPack's Compare Columns
  is missing the same way (`PowerPack/bdd/features/compare-columns/MISSING.md`). The feature's first scenario is a
  `@known-failure` on this; the other opens the dialog by calling `Curves:dataToCurvesTopMenu`.
