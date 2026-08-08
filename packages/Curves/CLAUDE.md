# CLAUDE.md

## Overview

**Curves** (`@datagrok/curves`) is a Datagrok plugin for **fitted curves** — dose-response curves, sigmoid fits, and general curve fitting. It provides in-grid cell rendering of charts stored in various formats in cells, interactive editing (outlier toggling, property panels), automatic curve fitting, a "Data to Curves" pipeline that converts well-level assay data into fitted curve summaries, a multi-curve overlay viewer, and statistics extraction (IC50, AUC, R², etc.).

Category: **Visualizations**. Top menu: `Data | Curves | ...`.

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

The **native format** for curve data in cells is JSON conforming to `IFitChartData`:
```
IFitChartData {
  chartOptions?: IFitChartOptions   // axes bounds, log scale, title, axis names, showStatistics
  seriesOptions?: IFitSeriesOptions // default series options (fit function, colors, markers)
  series?: IFitSeries[]             // array of data series
}
```
Each `IFitSeries` has `points: IFitPoint[]` (x, y, outlier, color, marker, stdev) plus fit configuration (fitFunction, parameters, colors, showFitLine, clickToToggle, droplines, etc.).

Labels are values the fit cannot produce - a plate's Z prime, an assay name, a compound id.
`chartOptions.labels` describes the whole plot and renders once; `series.labels` describes one curve
and renders with it. `chartOptions.showLabels` picks which names to draw, exactly as `showStatistics`
does for statistics, so it inherits the cascade and the persistence below for free. Unlike every other
option, `chartOptions.labels` merges **per key** across levels - labels are data, not a setting.

Options precedence, applied at render time in `fit/fit-chart-data.ts`:

**value the curve declares → dataframe explicit → column explicit → cell explicit**

Every level stores what it holds in a tag (`FitConstants.TAG_FIT`, dataframe and column) or in the
cell JSON, alongside an `explicit` list naming the options a user actually set there. An option that
is merely present fills a gap (`mergeProperties`); only an **explicitly set** one overrides the value
a series declares for itself. Setting a level clears the claims of the narrower levels below it, so a
Column change reaches every cell and a Dataframe change reaches every column - without deleting the
data's own values, which is how a setting survives a layout applied to a fresh table.

## Curve Converter Architecture

The Curves package supports **multiple curve data formats** through a dynamic converter system. The central format is the native IFitChartData JSON. Other formats (XML 3DX, compact dose-response, PZFX, etc.) are converted to native format on-the-fly via registered converters.

### How It Works

1. **Detection**: Semantic type detectors in `detectors.js` identify curve data in columns. Each detector:
   - Sets `col.semType = 'fit'`
   - For non-native formats, sets column tag `.%curve-format` to a friendly format name (e.g. `3dx`, `compact-dr`, `pzfx`)
   - If no format tag is set, the data is assumed to be in native IFitChartData JSON format

2. **Conversion**: When a cell value needs to be parsed, `parseCellValue()` in `curve-converter.ts`:
   - Checks the column's `.%curve-format` tag for the format name
   - If present, looks up the converter function by format name in the registry and applies it
   - If absent, parses the value directly as JSON
   - Results are cached in a value-based `DG.LruCache` keyed by `formatName||cellValue`

3. **Registration**: Converters are registered in `_initCurves()` via `registerCurveConverter(formatName, fn)`:
   - Local converters (in this package) are registered directly by format name to avoid circular dependencies
   - External converters (from other packages) are discovered via `initExternalConverters()` which finds DG functions tagged with `meta.role: curveConverter` and `meta.curveFormat: '<name>'`

### Key Constants

- `FitConstants.TAG_CURVE_FORMAT` = `'.%curve-format'` — column tag storing format name (e.g. `'3dx'`, `'compact-dr'`, `'pzfx'`)
- `FitConstants.FIT_SEM_TYPE` = `'fit'` — semantic type for fit curves
- `FitConstants.TAG_FIT_CHART_FORMAT` = `'.fitChartFormat'` — legacy tag for format identification
- `FitConstants.TAG_FIT` = `'.%fit'` — dataframe/column tag holding that level's options. The `.%`
  prefix is what makes the platform serialize it into a layout; `TAG_FIT_LEGACY` (`'.fit'`) is read
  and migrated once, then removed

### Supported Formats

| Format | Detector | Format Name | File |
|--------|----------|-------------|------|
| **Native JSON** (IFitChartData) | `detectFit` — looks for `series` + `points` | None (native) | — |
| **XML 3DX** (`<chart>...</chart>`) | `detectXMLCurveChart` | `3dx` | `fit/converters/xml-converter.ts` |
| **Compact Dose-Response** (`{"p":..., "poi":...}`) | `detectCompactDoseResponse` | `compact-dr` | `fit/converters/compact-dr-converter.ts` |
| **PZFX** (GraphPad Prism XY) | `detectPzfxCurveChart` | `pzfx` | `fit/converters/pzfx-converter.ts` |

### Adding a New Format

To add support for a new curve format, use the `/add-curve-format` skill. In summary:
1. Create a converter function in `src/fit/converters/` that converts the format to native JSON
2. Add a detector in `detectors.js` that identifies the format and sets the converter tag
3. Register the converter in `_initCurves()` in `package.ts`
4. Add tests in `src/tests/converter-tests.ts`

## Architecture

### Entry Point — `src/package.ts`

`PackageFunctions` class registers all platform-visible functions via `@grok.decorators`:

**Initialization:**
- `_initCurves` — `@init`: registers `FitGridCellHandler`, registers all local curve converters, discovers external converters

**Registered Functions:**

| Function | Decorator / Meta | Purpose |
|---|---|---|
| `curveFitDemo` | `@demo`, `demoPath: Curves \| Curve Fitting` | Loads demo CSV and opens table view |
| `dataToCurves` | `@func` | Programmatic API: converts well-level data into a DataFrame with fitted curve JSON column |
| `dataToCurvesTopMenu` | `@func`, `top-menu: Data \| Curves \| Data to Curves` | Opens the Data to Curves dialog UI |
| `addStatisticsColumn` | `@func`, `vectorFunc`, `role: transform` | Extracts a fit statistic from a specific series into a new float column |
| `addAggrStatisticsColumn` | `@func`, `vectorFunc`, `role: transform` | Aggregates statistics across all series |
| `convertXmlCurveToJsonFunc` | `@func`, `role: curveConverter` | Returns XML 3DX converter function |
| `convertCompactDrToJsonFunc` | `@func`, `role: curveConverter` | Returns compact dose-response converter function |
| `convertPzfxToJsonFunc` | `@func`, `role: curveConverter` | Returns PZFX converter function |

### Converter Infrastructure — `src/fit/curve-converter.ts`

Central module for the converter system:
- `registerCurveConverter(nqName, fn)` — registers a converter function
- `initExternalConverters()` — discovers converters from external packages via `DG.Func.find({meta: {role: 'curveConverter'}})`
- `parseCellValue(cellValue, column)` — **the single entry point** for parsing any cell value; applies converter if needed
- `convertCellValue(cellValue, column)` — converts using registered converter, with value-based LRU caching
- `isNativeFormat(column)` / `hasConverter(column)` — check whether a column uses native or converted format

### Semantic Type Detectors — `detectors.js`

`CurvesPackageDetectors` detects columns with semantic type `fit`:
1. `detectXMLCurveChart` — matches `<chart>...</chart>`, sets converter tag
2. `detectFit` — matches native JSON with `series` + `points`
3. `detectPzfxCurveChart` — matches PZFX XY table XML
4. `detectCompactDoseResponse` — matches compact DR JSON with `p` + `poi` fields

### Chart Data & Caches — `src/fit/fit-chart-data.ts`

A leaf module, so the renderer and the statistics do not import each other through it:
- **Caching**: three `DG.LruCache` instances: `parsedCurves` (parsed IFitChartData), `fittedCurves` (fitted curve results), `curvesDataPoints` (extracted data points)
- `getOrCreateParsedChartData(cell)` — cached parse with the precedence above applied
- `getColumnChartOptions` / `getDataFrameChartOptions` — read a level's options, migrating `.fit` onto `.%fit`
- `mergeProperties` (gap filling) and the explicit overrides that outrank a series' own value
- `mergeSeries`, `substituteZeroes`, `sanitizeCellValue`

### Cell Renderer — `src/fit/fit-renderer.ts`

`FitChartCellRenderer` (extends `DG.GridCellRenderer`) — the main grid cell renderer for `fit` cells:
- **render()**: parses cell via `parseCellValue()` → renders axes, fit lines, points, confidence intervals, droplines, statistics, title, legend
- **onClick()**: toggles outlier status on points (native format only), updates dependent statistics columns
- **onDoubleClick()**: opens `inspectCurve` dialog (resizable chart editor)
- **onMouseMove()**: shows tooltip with x/y values on point hover
- Outlier toggling is only available for native format columns (`isNativeFormat()`)

Key helper functions:
- `getOrCreateCachedFitCurve(series, ...)` — cached Nelder-Mead fit (in `fit-chart-data.ts`)
- `layoutChart(rect, ...)` and the `*Shown` predicates — geometry and what a cell of a given size has
  room for, in `fit-layout.ts`
- `setOutlier`, `inspectCurve`, `handleClick`, `handleMouseMove` — what a user does to a curve, in
  `fit-interaction.ts`. The renderer's `onClick`/`onMouseMove` delegate to them; `handleMouseMove`
  takes the renderer as the base `DG.GridCellRenderer` type so the two modules do not import each other

Parameters stored in a cell are **always data space**; `seriesInFitSpace()` hands the renderer a copy
in the optimizer's coordinates for the curve, the confidence band and the droplines, and everything
shown as a number goes through `toDataSpace()` instead. See "Fit statistics" below.

### Fit statistics

Statistics come from the typed fit, not positional parameter indices: `getSeriesFit()` → a `Fit`
subclass, read by name with `getStatistic(fit, name)`. Each fit function advertises its own set via
`statisticsProperties`, so a linear fit has no `ic50` and asking returns `undefined`, not `NaN`.
`toDataSpace()` converts log-fitted statistics back to data space exactly once — **aggregate before
calling it**, so an averaged IC50 stays a geometric mean.

Legacy names (`interceptX`, `top`, …) are persisted in `showStatistics`, in the `.statistics` column
tag and in recorded transforms; `LEGACY_FIT_STATISTICS` is append-only and resolved through an alias
table. `addStatisticsColumn` and the `.sourceColumn` sweep in `setOutlier` are the legacy path, kept
so existing projects keep working.

`curveStatistic` / `curveAggrStatistic` add statistics as **calculated columns** (`meta.vectorFunc` +
`action: join(table)`, called with `processed: false`). Three traps there fail silently — a parameter
must not be named after the package module global, the result column name must be stable, and on
recalculation the function gets a detached column. See `src/tests/calculated-columns-tests.ts`.

### Property Panel — `src/fit/fit-grid-cell-handler.ts`

`FitGridCellHandler` (extends `DG.ObjectHandler`) — renders the context panel when a fit cell is selected:
- **Options pane**: switch level (Dataframe/Column/Cell), series and chart options
- **Chart pane**: enlarged `GridCellWidget` rendering
- **Fit pane**: per-series or aggregated statistics with "+" buttons to extract stats as columns
- `addStatisticColumn()` runs `curveStatistic`/`curveAggrStatistic` and puts the result right after
  the curve column, in the dataframe and in the grid, which keep their order separately

### Options — `src/fit/fit-options.ts`

Writing an option at any level, kept apart from the panel that renders it:
- `changeCurvesOptions(gridCell, input, section, level)` — the single entry point
- `claim`/`unclaim` — record or drop "the user set this here", which is what outranks the data
- `chartPropertiesFor` / `seriesPropertiesFor` / `normalizeStatisticNames` — the property lists the
  panel binds to, built from the cell's own fit functions rather than a fixed legacy list
- `STATISTIC_AFFECTING_OPTIONS` gates the notification: only these refit dependent statistic columns,
  everything else notifies at dataframe level, which marks the table modified without recalculating

### Rendering Utilities — `src/fit/render-utils.ts`

Low-level canvas drawing functions using `CanvasRenderingContext2D`:

| Function | What it renders |
|---|---|
| `renderPoints` | Data points (circles, markers) and error bars |
| `renderFitLine` | Fitted curve line (solid/dotted/dashed/dashdotted) |
| `renderConfidenceIntervals` | Shaded confidence bands |
| `renderConnectDots` | Direct point-to-point line connections |
| `renderDroplines` | ICxx / ECxx droplines, resolved through `FitFunction.inverse` |
| `renderStatistics` | In-chart text labels for selected statistics |
| `renderTitle` | Chart title |
| `renderAxesLabels` | X and Y axis name labels |

### Legend — `src/fit/fit-legend.ts`

The legend annotates the plot from a corner, over a backdrop, taking at most 40% of the width and
half the height. `chooseLegendCorner()` picks the corner: `renderPoints` and `renderFitLine` report
where they put ink (`drawnAt`), each corner is scored by how much of it falls within a clearance ring,
ties go to the corner the ink keeps furthest away, and the top right one stays unless another is less
than half as crowded - otherwise the legend would hop about as rows are hovered or the grid scrolls.

A column that holds a single curve gets no row of its own; a name that two columns share is qualified
with the column, while the same name twice in one column is left alone (the prefix would be identical);
a name too long is ellipsized; rows past the corner collapse into `+N more`. `legendTooltip()` answers
the hover - the whole of a shortened name, or every row when the pointer is on `+N more` - reading the
layout the render recorded, since only the render knows which corner it chose. `isLegendVisible()` is
the single gate (the `showLegend` chart option and the size thresholds), used by the renderer and the
hover alike. Rows are measured on a canvas context of its own, since hit testing has none to draw on.

### Converters — `src/fit/converters/`

Each converter is a pure function: `(value: string) => string` (raw cell value → native JSON string).

| File | Function | Source Format |
|------|----------|---------------|
| `xml-converter.ts` | `convertXmlCurveToJson` | XML 3DX `<chart>` documents |
| `compact-dr-converter.ts` | `convertCompactDrToJson` | Compact dose-response JSON (`{p, poi, mxr, mnr, hs, ...}`) |
| `pzfx-converter.ts` | `convertPzfxToJson` | GraphPad Prism PZFX XY tables |

### XML Parser — `src/fit/fit-parser.ts`

`convertXMLToIFitChartData(xmlText)` — parses legacy 3DX XML format `<chart>` documents into `IFitChartData`. Used by the XML converter. Reorders parameters from `[IC50, tan, max, min]` to `[max, tan, IC50, min]`.

### Data to Curves — `src/fit/data-to-curves.ts`

The main data transformation pipeline. `convertDataToCurves()` groups well-level assay data by compound|assay|target|run and builds fitted curve JSON.

### Multi-Curve Viewer — `src/fit/multi-curve-viewer.ts`

`MultiCurveViewer` (extends `DG.JsViewer`) — superimposes curves from multiple rows/columns on one chart. Supports trellis mode, selection/current/mouseover tracking, and legend columns.

### PZFX Parser — `src/formats/pzfx/pzfx-parser.ts`

Parses GraphPad Prism `.pzfx` files. Also used by the PZFX converter and the file viewer/handler.

## Tests — `src/tests/`

Test entry point: `src/package-test.ts` — imports all test files.

| File | What it tests |
|---|---|
| `detector-tests.ts` | Semantic type detection for all formats (JSON, XML, compact DR, PZFX) |
| `converter-tests.ts` | All converters (XML, compact DR, PZFX), caching, format detection, rendering integration |
| `fit-tests.ts` | Core fitting algorithms (sigmoid, linear, log-linear, exponential, polynomial) |
| `curves-cell-renderer-tests.ts` | Cell renderer creation and rendering performance |
| `calculated-columns-tests.ts` | `curveStatistic`/`curveAggrStatistic` as calculated columns, statistic-space round trips, aggregation |
| `panel-renderer-tests.ts` | Property panel construction, option precedence between levels, notification gating |
| `transform-tests.ts` | Viewport coordinate transformations |
| `pzfx-tests.ts` | PZFX file parser |
| `legend-tests.ts` | Legend rows, the reserved strip, ellipsis and overflow |
| `multi-curve-viewer-tests.ts` | The viewer works on copies of the cached cell data |
| `curve-data.ts` | Shared sigmoid data and `renderedTexts()` used by the suites above (not a suite) |

## Source Structure

```
src/
  package.ts              — Entry point, PackageFunctions, converter registration
  package.g.ts            — Auto-generated function wrappers (DO NOT EDIT)
  package-api.ts          — Auto-generated API namespace (DO NOT EDIT)
  package-test.ts         — Test entry point
  fit/
    curve-converter.ts    — Converter registry, LRU cache, parseCellValue() entry point
    fit-chart-data.ts     — Options precedence, parsed-chart-data and fit caches (leaf module)
    fit-layout.ts         — Chart geometry, and what fits in a cell of a given size
    fit-legend.ts         — Legend rows, the strip they are given, and their drawing
    fit-renderer.ts       — FitChartCellRenderer: assembling a chart onto a canvas
    fit-interaction.ts    — Outlier toggle, tooltip, click handling, the chart editor dialog
    fit-grid-cell-handler.ts — FitGridCellHandler (property panel), stat column extraction
    fit-options.ts        — changeCurvesOptions, claims, panel property lists
    fit-statistics.ts     — Per-series and aggregated statistics calculation
    render-utils.ts       — Canvas drawing functions (points, lines, CIs, droplines, statistics)
    fit-parser.ts         — XML 3DX → IFitChartData parser (used by xml-converter)
    data-to-curves.ts     — Data to Curves pipeline (UI dialog + conversion logic)
    multi-curve-viewer.ts — MultiCurveViewer (overlay viewer)
    fit-inspection.ts     — inspectSeries helper
    fit-demo.ts           — Demo data generation and demo entry point
    converters/
      xml-converter.ts    — XML 3DX format converter
      compact-dr-converter.ts — Compact dose-response JSON converter
      pzfx-converter.ts   — PZFX (GraphPad Prism) format converter
  formats/
    pzfx/
      pzfx-parser.ts      — PZFX file parser
  tests/
    curve-data.ts         — Shared sigmoid test data
    detector-tests.ts     — Semantic type detector tests (all formats)
    converter-tests.ts    — Converter tests (all formats, caching, integration)
    fit-tests.ts          — Fitting algorithm tests
    curves-cell-renderer-tests.ts — Renderer benchmarks
    transform-tests.ts    — Viewport transform tests
    pzfx-tests.ts         — PZFX parser tests
scripts/
  calculateMSR.py         — MSR calculation (Python, statsmodels)
files/
  curves-demo.csv         — Demo data
detectors.js              — Semantic type detectors (all curve formats)
```

## Quick Lookups

| Looking for... | Check first |
|---|---|
| Adding a new curve format | Use `/add-curve-format` skill |
| Converter infrastructure | `src/fit/curve-converter.ts` |
| Existing converters | `src/fit/converters/` |
| Function/demo/init registration | `src/package.ts` |
| Cell renderer (drawing) | `src/fit/fit-renderer.ts` |
| Click, tooltip, outlier toggle, editor dialog | `src/fit/fit-interaction.ts` |
| Chart geometry, size thresholds | `src/fit/fit-layout.ts` |
| Legend rows, width and drawing | `src/fit/fit-legend.ts` |
| Options precedence, parse and fit caches | `src/fit/fit-chart-data.ts` |
| Property panel (options, stats, chart) | `src/fit/fit-grid-cell-handler.ts` |
| Writing an option at a level, claims | `src/fit/fit-options.ts` |
| Statistics calculation (per-series, aggregated) | `src/fit/fit-statistics.ts` |
| Canvas drawing (points, lines, axes, statistics, labels) | `src/fit/render-utils.ts` |
| Data to Curves dialog & conversion | `src/fit/data-to-curves.ts` |
| Multi-curve overlay viewer | `src/fit/multi-curve-viewer.ts` |
| XML 3DX format parsing | `src/fit/fit-parser.ts` |
| Fit data types (IFitChartData, etc.) | `@datagrok-libraries/statistics/src/fit/fit-curve.ts` |
| Fit functions, typed fits, optimizer | `.../fit/fit-engine.ts` |
| Series-level API (getSeriesFit, toDataSpace) | `.../fit/fit-data.ts` |
| Constants (tags, sizes, thresholds) | `@datagrok-libraries/statistics/src/fit/const.ts` |
| Semantic type detection | `detectors.js` |
| Test entry point | `src/package-test.ts` |
