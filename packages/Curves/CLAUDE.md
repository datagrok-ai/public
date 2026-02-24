# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**Curves** (`@datagrok/curves`) is a Datagrok plugin for **fitted curves** — dose-response curves, sigmoid fits, and general curve fitting. It provides in-grid cell rendering of charts stored as JSON/XML in cells, interactive editing (outlier toggling, property panels), automatic curve fitting, a "Data to Curves" pipeline that converts well-level assay data into fitted curve summaries, a multi-curve overlay viewer, and statistics extraction (IC50, AUC, R², etc.).

Category: **Visualizations**. Top menu: `Data | Curves | ...`.

## Build Commands

```bash
npm install
npm run build              # grok api && grok check --soft && webpack
npm run test               # grok test
npm run build-all          # Builds js-api → utils → statistics → this package
npm run link-all           # Links datagrok-api, @datagrok-libraries/utils, @datagrok-libraries/statistics
```

## Key Dependencies

- `@datagrok-libraries/statistics` — **core fitting engine**: defines `IFitChartData`, `IFitSeries`, `IFitChartOptions`, `IFitPoint`, `FitStatistics`, `FitCurve`, `FitFunction` classes, Nelder-Mead optimizer, sigmoid/linear/log-linear/exponential/4PL fit functions, confidence intervals, box-plot statistics, and `FitConstants` (tags, sizing, thresholds)
- `@datagrok-libraries/utils` — `Viewport` (world↔screen coordinate transform), `StringUtils`
- `@datagrok-libraries/test` — test framework (`category`, `test`, `expect`)
- `jstat` — statistical distributions
- `wu` — lazy iteration
- `exceljs` — Excel file support
- `dayjs` — date manipulation
- `rxjs` — reactive subscriptions (debounce, merge)

### Data Model (from `@datagrok-libraries/statistics`)

Cell values are JSON strings conforming to `IFitChartData`:
```
IFitChartData {
  chartOptions?: IFitChartOptions   // axes bounds, log scale, title, axis names, showStatistics
  seriesOptions?: IFitSeriesOptions // default series options (fit function, colors, markers)
  series?: IFitSeries[]             // array of data series
}
```
Each `IFitSeries` has `points: IFitPoint[]` (x, y, outlier, color, marker, stdev) plus fit configuration (fitFunction, parameters, colors, showFitLine, clickToToggle, droplines, etc.).

Legacy XML format (3DX `<chart>...</chart>`) is auto-detected and converted to JSON via `fit-parser.ts`.

Options cascade: **DataFrame tags → Column tags → Cell JSON → Series defaults**, merged at render time.

## Architecture

### Entry Point — `src/package.ts`

`PackageFunctions` class registers all platform-visible functions via `@grok.decorators`:

**Initialization:**
- `_initCurves` — `@init`: registers `FitGridCellHandler` as an `ObjectHandler` for grid cell property panels

**Registered Functions:**

| Function | Decorator / Meta | Purpose |
|---|---|---|
| `curveFitDemo` | `@demo`, `demoPath: Curves \| Curve Fitting` | Loads demo CSV and opens table view |
| `dataToCurves` | `@func` | Programmatic API: converts well-level data (concentration, readout, batch/compound/assay/run IDs) into a DataFrame with fitted curve JSON column. Supports parent-table join for pre-fit parameters, reported IC50s, etc. |
| `dataToCurvesTopMenu` | `@func`, `top-menu: Data \| Curves \| Data to Curves` | Opens the Data to Curves dialog UI |
| `addStatisticsColumn` | `@func`, `vectorFunc`, `role: transform` | Extracts a single fit statistic (e.g., IC50, AUC, R²) from a specific series into a new float column |
| `addAggrStatisticsColumn` | `@func`, `vectorFunc`, `role: transform` | Same but aggregates across all series using a DG.Stats aggregation (min, max, avg, med, etc.) |

**Exports**: `_package`, `Sync` (sequential promise runner), auto-generated wrappers from `package.g.ts`

### Semantic Type Detector — `detectors.js`

`CurvesPackageDetectors` detects columns with semantic type `fit`:
1. `detectXMLCurveChart` — matches strings starting with `<chart>` and ending with `</chart>`, sets tag `.fitChartFormat = '3dx'`
2. `detectFit` — matches JSON strings containing both `series` and `points` keywords

### Cell Renderer — `src/fit/fit-renderer.ts`

`FitChartCellRenderer` (extends `DG.GridCellRenderer`) — the main grid cell renderer for `fit` cells:
- **Caching**: three `DG.LruCache` instances: `parsedCurves` (parsed JSON), `fittedCurves` (fitted curve results), `curvesDataPoints` (extracted data points)
- **render()**: parses cell → gets chart data → renders axes, fit lines, points, confidence intervals, droplines, statistics, title, legend
- **onClick()**: toggles outlier status on points (if `clickToToggle` is enabled), updates dependent statistics columns
- **onDoubleClick()**: opens `inspectCurve` dialog (resizable chart editor)
- **onMouseMove()**: shows tooltip with x/y values on point hover; shows enlarged chart tooltip for small cells
- **Visibility thresholds**: axes, labels, title, legend, points, droplines each have minimum pixel-size thresholds from `FitConstants`

Key helper functions:
- `getOrCreateParsedChartData(cell)` — cached JSON parse with options merge (DF → Column → Cell cascade)
- `getOrCreateCachedFitCurve(series, ...)` — cached Nelder-Mead fit
- `substituteZeroes(data)` — replaces x=0 with calculated substitute when logX is enabled
- `setOutlier(gridCell, point, ...)` — toggles outlier, updates cell JSON, fires custom event, recalculates stats columns
- `inspectCurve(gridCell)` — opens a dialog with a `GridCellWidget` for interactive chart editing
- `layoutChart(rect, ...)` — computes viewport, xAxis, yAxis rectangles with margins
- `mergeChartOptions(...)` / `mergeSeries(...)` — combine options/series from multiple cells
- `mergeProperties(properties, source, target)` — cascading property merge

### Property Panel — `src/fit/fit-grid-cell-handler.ts`

`FitGridCellHandler` (extends `DG.ObjectHandler`) — renders the context panel when a fit cell is selected:
- **Options pane**: switch level (Dataframe/Column/Cell), series options (fit function, colors, markers, line style, confidence intervals), chart options (log scale, axes bounds, zeroes handling)
- **Chart pane**: enlarged `GridCellWidget` rendering
- **Fit pane**: displays per-series or aggregated statistics (R², AUC, IC50, slope, top, bottom) with "+" buttons to extract each stat as a new column

Also contains:
- `calculateSeriesStats()` — computes `FitStatistics` for a single series (uses cached fit curve)
- `getChartDataAggrStats()` — aggregates statistics across all series using DG.Stats
- `changeCurvesOptions()` — applies property changes at Dataframe/Column/Cell level
- `detectSettings()` — scans column to detect which options are custom per-cell vs default

### Rendering Utilities — `src/fit/render-utils.ts`

Low-level canvas drawing functions, all using `CanvasRenderingContext2D`:

| Function | What it renders |
|---|---|
| `renderPoints` | Data points (circles, markers) and error bars (stdev whiskers) |
| `renderFitLine` | Fitted curve line (with line style: solid/dotted/dashed/dashdotted) |
| `renderConfidenceIntervals` | Shaded confidence bands around fit line |
| `renderConnectDots` | Direct point-to-point line connections |
| `renderDroplines` | IC50 droplines (dashed lines from axis to curve) |
| `renderStatistics` | In-chart text labels for selected statistics |
| `renderTitle` | Chart title above the plot area |
| `renderAxesLabels` | X and Y axis name labels |
| `renderLegend` | Color-coded legend with series names, markers, and optional column labels |

Also: `getSeriesColor()`, candlestick/box-plot rendering (for replicate data visualization).

### XML Parser — `src/fit/fit-parser.ts`

`convertXMLToIFitChartData(xmlText)` — parses legacy 3DX XML format `<chart>` documents into `IFitChartData`. Extracts grid bounds, axis labels, series parameters (reorders from `[IC50, tan, max, min]` to `[max, tan, IC50, min]`), points, and outlier masks.

### Data to Curves — `src/fit/data-to-curves.ts`

The main data transformation pipeline:

**`dataToCurvesUI()`** — dialog with two forms:
- **Well Table Form**: table, assay, concentration, readout, batch ID, run ID, compound ID, target entity, outliers, additional columns (auto-detects columns by name heuristics)
- **Parent Table Form**: optional table with pre-fit parameters (Max/Hill/Inflection/Min), reported IC50, qualified IC50, experiment ID, qualifier, fit function selection (4PL dose-response or 4PL regression), use-prefit-params toggle
- **Join Editor**: columns to join well-level and parent-level data
- Supports dialog history for re-running with same settings

**`convertDataToCurves()`** — core logic:
1. Groups well-level rows by compound|assay|target|run (or custom join column)
2. Builds `IFitSeries` per curve with points sorted by x, outlier flags, per-run marker types, categorical colors
3. Assigns pre-fit parameters from parent table if available
4. Creates result DataFrame with: `Fitted Curve` column (JSON, semType=fit), `Compound|Assay|Target` grouping column, info columns (assay, batch, compound, target, run, reported IC50, etc.)
5. Computes aggregated stats: max percent inhibition, geomean IC50, standard deviation, total/reported result counts
6. Auto-adds calculated statistics columns (IC50, Max, Min, Hill, AUC) if no parent table
7. Sets up layout: trellis plot with MultiCurveViewer, selected curves grid, overlay viewer

### Multi-Curve Viewer — `src/fit/multi-curve-viewer.ts`

`MultiCurveViewer` (extends `DG.JsViewer`) — superimposes curves from multiple rows/columns on one chart:
- Properties: `curvesColumnNames` (column list), `legendColumnName`, `showSelectedRowsCurves`, `showCurrentRowCurve`, `showMouseOverRowCurve`, `mergeColumnSeries`, `showOutliers`, plus all `IFitChartOptions` properties
- Reacts to: current cell change, mouse-over row, selection change (debounced 50ms)
- Merges chart options from all visible cells, assigns categorical colors, adjusts opacity for large selections
- Supports trellis mode (uses filter instead of selection/current)
- `fromChartData(chartData)` — static factory for standalone use

### Fit Inspection — `src/fit/fit-inspection.ts`

`inspectSeries(series, fitFunctionName)` — creates a temporary single-row DataFrame with the series as a fit cell, opens it in the inspect dialog with click-to-toggle and IC50 droplines. Used for programmatic curve inspection.

### Demo — `src/fit/fit-demo.ts`

- `createSigmoidPoints(length, step, pointsPerX)` — generates random sigmoid data with noise
- `createDemoDataFrame(rowCount, chartsCount, chartsPerCell)` — builds a full demo DataFrame with multiple fit columns
- `curveDemo()` — loads `files/curves-demo.csv` and opens it

### Python Script — `scripts/calculateMSR.py`

`Calculate MSR` — computes Minimum Significant Ratio using mixed-effects linear models (statsmodels). Groups by compound, assay, target entity; uses log10(IC50) with run dates as random effects. Falls back to simple variance calculation for small sample sizes. Top menu: `Data | Curves | Calculate MSR...`.

## Tests — `src/tests/`

Test entry point: `src/package-test.ts` — imports all test files.

| File | What it tests |
|---|---|
| `fit-tests.ts` | Core fitting: `getSeriesFitFunction`, `getCurve`, `fitSeries` (sigmoid, linear, log-linear, exponential, polynomial), confidence intervals, statistics, benchmarks |
| `detector-tests.ts` | Semantic type detection for JSON fit data and XML 3DX format |
| `curves-cell-renderer-tests.ts` | Cell renderer creation and rendering performance benchmarks |
| `transform-tests.ts` | Viewport coordinate transformations (world↔screen) |

## Source Structure

```
src/
  package.ts              — Entry point, PackageFunctions, Sync helper
  package.g.ts            — Auto-generated function wrappers (DO NOT EDIT)
  package-api.ts          — Auto-generated API namespace (scripts, funcs)
  package-test.ts         — Test entry point
  fit/
    fit-renderer.ts       — FitChartCellRenderer, caching, chart layout, merging
    fit-grid-cell-handler.ts — FitGridCellHandler (property panel), stats calculation
    render-utils.ts       — Canvas drawing functions (points, lines, CIs, droplines, legend)
    fit-parser.ts         — XML 3DX → IFitChartData converter
    data-to-curves.ts     — Data to Curves pipeline (UI dialog + conversion logic)
    multi-curve-viewer.ts — MultiCurveViewer (overlay viewer)
    fit-inspection.ts     — inspectSeries helper
    fit-demo.ts           — Demo data generation and demo entry point
  tests/
    fit-tests.ts          — Fitting algorithm tests
    detector-tests.ts     — Semantic type detector tests
    curves-cell-renderer-tests.ts — Renderer benchmarks
    transform-tests.ts    — Viewport transform tests
scripts/
  calculateMSR.py         — MSR calculation (Python, statsmodels)
files/
  curves-demo.csv         — Demo data
detectors.js              — Semantic type detectors (XML and JSON fit formats)
```

## Quick Lookups

| Looking for... | Check first |
|---|---|
| Function/demo/init registration | `src/package.ts` |
| Cell renderer (render, click, tooltip) | `src/fit/fit-renderer.ts` (`FitChartCellRenderer`) |
| Property panel (options, stats, chart) | `src/fit/fit-grid-cell-handler.ts` (`FitGridCellHandler`) |
| Canvas drawing (points, lines, axes) | `src/fit/render-utils.ts` |
| Data to Curves dialog & conversion | `src/fit/data-to-curves.ts` |
| Multi-curve overlay viewer | `src/fit/multi-curve-viewer.ts` |
| XML 3DX format parsing | `src/fit/fit-parser.ts` |
| Inspect/edit a single curve | `src/fit/fit-inspection.ts` + `inspectCurve` in `fit-renderer.ts` |
| Fit data types (IFitChartData, etc.) | `@datagrok-libraries/statistics/src/fit/fit-curve.ts` |
| Fit functions & optimizer | `@datagrok-libraries/statistics/src/fit/fit-data.ts` |
| Constants (tags, sizes, thresholds) | `@datagrok-libraries/statistics/src/fit/const.ts` |
| New fit API (FitFunction classes) | `@datagrok-libraries/statistics/src/fit/new-fit-API.ts` |
| Viewport (coordinate transforms) | `@datagrok-libraries/utils/src/transform.ts` |
| Semantic type detection | `detectors.js` |
| MSR calculation (Python) | `scripts/calculateMSR.py` |
| Auto-generated function wrappers | `src/package.g.ts` / `src/package-api.ts` |
| Test entry point | `src/package-test.ts` |
| Demo data | `files/curves-demo.csv` |
