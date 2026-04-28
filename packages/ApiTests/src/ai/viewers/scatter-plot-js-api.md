# ScatterPlotViewer JS API — features-subdir focus

**JS API source:** `public/js-api/src/viewer.ts:611` (ScatterPlotViewer)

## What we are testing

The Scatter Plot is the richest viewer in `core/client/d4/lib/src/viewers/scatterplot/features/` — formula lines, annotation regions, regression lines, drop lines, smart labels — with most of its surface exposed via JS through `ScatterPlotViewer` (`viewer.ts:611-669`) and the shared `ViewerMetaHelper` (`viewer.ts:769-779`). This scenario targets JS-only entry points: the `df.plot.scatter` factory shorthand, the helper-based `meta.formulaLines` / `meta.annotationRegions` add/items/clear API, the typed `zoom()` method, the seven Scatter-specific event Observables, and `getInfo()` / `viewBox`-family geometry accessors that only `ScatterPlotViewer` exposes.

## Why it is uncovered

`src/ai/viewers/` covers BarChart, Histogram, LineChart and PieChart; `src/grid/` and the rest of `src/` mention scatter only through `setProperty` smoke tests. Dart-side `core/client/xamgle/lib/src/tests/viewers_test.dart` covers `xColumnName` / `yColumnName` / `xAxisType` / `yAxisType` / a raw-string `formulaLines` round-trip and an OOM regression — but never the JS helper façade (`FormulaLinesHelper.add` / `.items` / `.clear`, same for `AnnotationRegionsHelper`), the `df.plot.scatter` shorthand, the `zoom()` JS method, the event Observable shape, or the `viewBox` / `xAxisBox` / `yAxisBox` getters. The Dart parametric prop iteration handles bare bool/enum/string round-trips, so this scenario stays clear of those.

## Preconditions

- Demog dataset reachable via `grok.data.demo.demog(N)`.
- ApiTests package built and published to `dev`.
- No specific server entities required — every test creates and disposes its own viewer / table view.

## Test cases

1. **df.plot.scatter shorthand returns typed ScatterPlotViewer** — call `df.plot.scatter({xColumnName, yColumnName, colorColumnName})`; expect `instanceof DG.ScatterPlotViewer`, `type === VIEWER.SCATTER_PLOT`, friendly-key aliasing reflected in `props` and `getOptions(true).look`.
2. **meta.formulaLines helper round-trip** — `viewer.meta.formulaLines.add({...formula line...})`; expect `items` returns the parsed array, raw `props['formulaLines']` is the JSON-stringified storage, `clear()` empties both.
3. **meta.annotationRegions helper round-trip** — same shape as case 2 against `viewer.meta.annotationRegions` (uses a separate `AnnotationRegionsHelper`); region constructed with `header` / `description` / `opacity` per the `AnnotationRegion` interface in `js-api/src/helpers.ts:37`.
4. **zoom() reflects xMin/xMax/yMin/yMax in look** — call `viewer.zoom(x1, y1, x2, y2)` on a `view.addViewer`-attached scatter plot; expect no throw and that subsequent `getOptions(true).look` exposes the four `xMin`/`xMax`/`yMin`/`yMax` keys (numeric values not asserted — pre-paint the platform leaves them `null`; documented fallback applied per Author's first-paint contingency).
5. **Scatter event Observables are subscribable** — assert `onZoomed`, `onResetView`, `onViewportChanged`, `onAfterDrawScene`, `onBeforeDrawScene`, `onPointClicked`, `onPointDoubleClicked` are rxjs `Observable`s with working `.subscribe(...).unsubscribe()`.
6. **view.addViewer typed instance + viewBox/getInfo geometry** — `tv.addViewer(VIEWER.SCATTER_PLOT, ...)` returns `DG.ScatterPlotViewer` present in `tv.viewers`; `viewBox` is a `DG.Rect`; `xAxisBox` / `yAxisBox` are non-null at first-paint (asymmetric `toJs` wrapping observed: `viewBox` lands as a `DG.Rect` instance, the axis-box getters return raw objects pre-paint); `getInfo()` exposes `canvas` and `overlay` keys. Width/height numeric assertions dropped per Author's first-paint contingency.

Negative cases: skipped — the scatter factory accepts unknown property keys silently (JS-side `setOptions` swallows them), `zoom()` is a void method that has no documented invalid input, and the helper APIs accept any JS object as a formula-line / region (validation lives Dart-side beyond this layer's contract).

## Out of scope

- Dart-side property defaults, regression line math, layout serialization — covered by `viewers_test.dart` and the parametric prop iteration.
- Pointer event payloads (`onPointClicked` data shape under a real click) — would require synthetic pointer dispatch; not safe in a headless ApiTests run.
- `screenToWorld` / `worldToScreen` / `pointToScreen` — round 8 surfaced a `NaN.floor()` failure on `LineChartViewer.screenToWorld(0, 0)` before first paint, same risk applies here.
