# LineChart JS API

**JS API source:** `public/js-api/src/viewer.ts:243` (`DG.Viewer.lineChart`), `public/js-api/src/viewer.ts:577` (`DG.LineChartViewer`)

## What we are testing

LineChart-specific JS API surface that the existing Dart parametric `@Prop` round-trip tests cannot easily pin down: factory friendly-key aliasing on `DG.Viewer.lineChart` / `df.plot.line`, the multi-element array properties (`yColumnNames`, `splitColumnNames`, `chartTypes`) that round-trip through `getOptions(true).look`, the LineChart-typed events (`onLineSelected` carrying `LineChartLineArgs`, `onZoomed`) as rxjs Observables, and the LineChart-typed `screenToWorld(...)` / `activeFrame` accessors on the JS wrapper.

## Why it is uncovered

`src/widgets/<viewer>` only covers a handful of BarChart prop setters; there is no LineChart coverage on the JS side. The Dart parametric tests in `core/client/xamgle/lib/src/tests/viewers_test.dart` round-trip every `@Prop` for every viewer, and dedicated Dart tests already exercise `xColumnName`/`yColumnName`/`formulaLines` and column-rename for LineChart — those are deliberately out of scope here. Round 5 (Histogram) covered the generic `getOptions(includeDefaults)` envelope, dataFrame swap, events-as-Observable shape and `close()` lifecycle as generic patterns; round 6 (PieChart) covered `view.addViewer` typed-instance attach and `setOptions` round-trip; round 7 (BarChart) covered `df.plot.bar` shorthand, `resetView()` and `onResetView`. This scenario picks LineChart-specific JS surface that remains uncovered after rounds 5–7.

## Preconditions

- A running Datagrok instance with the JS API loaded.
- The demog demo dataset is reachable through `grok.data.demo.demog(...)` (used by every viewer test on the JS side).
- No package-level state is created on the server; viewers built via `DG.Viewer.lineChart` / `df.plot.line` without `tv.addViewer` stay detached and are GC'd. The `view.addViewer` case opens a `TableView` and closes it via `tv.close()` in a `finally` block (round-5 lifecycle precaution: do not call `viewer.close()` on a never-attached viewer — it throws).

## Test cases

1. **factory friendly-key aliasing via DG.Viewer.lineChart** — call `DG.Viewer.lineChart(df, {x: 'age', yColumnNames: ['height'], split: 'race'})`. Expect `v instanceof DG.LineChartViewer`, `v.type === DG.VIEWER.LINE_CHART`, and the friendly key `x` to land on canonical `xColumnName`, `split` (deprecated singular) on `splitColumnName`, and `yColumnNames` array round-tripped. Verify the `look` block from `v.getOptions(true)` agrees. (`y: 'height'` is intentionally NOT used: `ILineChartSettings` does not declare a `y` field, and the Dart `setOptions` alias resolution maps `'y'` → `'yColumnName'` (not present on LineChart) → no plural fallback because `'y'` doesn't end with `'s'`.)
2. **factory via df.plot.line shorthand** — call `df.plot.line({xColumnName: 'age', yColumnNames: ['weight']})`. Expect a `DG.Viewer` with `type === DG.VIEWER.LINE_CHART`, `dataFrame === df`, and `xColumnName`/`yColumnNames` to round-trip — the shorthand is a thin wrapper over `DG.Viewer.lineChart`.
3. **yColumnNames and splitColumnNames array round-trip** — call `DG.Viewer.lineChart(df, {xColumnName: 'age', yColumnNames: ['height', 'weight'], splitColumnNames: ['race', 'sex']})`. Expect both arrays to round-trip as arrays through `v.props['yColumnNames']` and through `v.getOptions(true).look['yColumnNames']` (same for `splitColumnNames`). The Dart parametric tests do not assert array shape on multi-Y/multi-split.
4. **multiAxis and splineTension round-trip via setOptions** — start from `{xColumnName, yColumnNames: ['height', 'weight']}`, then `setOptions({multiAxis: true, splineTension: 0.5})`. Both are `@Prop`s on `LineChartLook`, so `setProperty` accepts them; verify `v.props.multiAxis === true`, `v.props.splineTension === 0.5`, the same in `getOptions(true).look`, and that `getProperties()` lists both. (`chartTypes` was originally proposed here but is **not** a `@Prop` on `LineChartLook` — it is a plain field that `setProperty` cannot set and `toMap()` does not serialize, so it doesn't round-trip through `setOptions`/`getOptions`.)
5. **onLineSelected and onZoomed are rxjs Observables (LineChart-typed events)** — build a detached LineChart, assert `onLineSelected`, `onZoomed`, and `onResetView` are Observable-shaped (`subscribe` is a function, returns a `Subscription` whose `unsubscribe` is a function). No synthesized canvas events — the wrapper presence is enough; the event payload type (`EventData<LineChartLineArgs>` / `null`) is enforced statically by TypeScript.
6. **view.addViewer attaches a typed LineChartViewer with activeFrame** — open a `TableView`, `tv.addViewer(DG.VIEWER.LINE_CHART, {xColumnName: 'age', yColumnNames: ['height']})`. Expect the returned viewer to be `instanceof DG.LineChartViewer`, present in `tv.viewers`, and expose `activeFrame` (non-null) materialisable via `DG.toJs(af) instanceof DG.DataFrame`. Close the table view in `finally`. **Platform note:** `LineChartViewer.activeFrame` getter (`viewer.ts:582`) returns the raw Dart handle directly — no `toJs` wrap on the JS side. The test materialises the wrapper via `DG.toJs(af)` and notes this as a platform inconsistency to surface. `screenToWorld(x, y)` was originally part of this case but **dropped** — it requires the viewer to have laid out (the Dart `LineChartCore.toWorld` reads viewport `Rect` whose initial values are `NaN` and only stabilise after first paint). Without a deterministic post-paint signal in headless ApiTests, the call throws `NaN.floor()`. Worth a JIRA: `screenToWorld` should either return `null` or wait for first layout instead of throwing.

Negative cases skipped: `DG.Viewer.lineChart` accepts arbitrary partial options without an observable failure mode (unknown keys are dropped silently by the Dart side).

## Out of scope

- `xColumnName`/`yColumnName`/`formulaLines` round-trip (covered by Dart parametric and dedicated Dart tests in `core/client/xamgle/lib/src/tests/viewers_test.dart`).
- Column-rename behaviour for LineChart (Dart side).
- `getOptions(includeDefaults)` JSON envelope shape, `getInfo`/`getProperties` general structure, `dataFrame` swap, generic events-as-Observable shape, `close()` lifecycle (round 5 / Histogram).
- `setOptions` generic round-trip pattern (round 6 / PieChart).
- `df.plot.<x>` factory generic + `resetView` / `onResetView` (round 7 / BarChart — both already exercised on a different viewer).
- Canvas event simulation, zoom gesture, line selection by user input — out of scope for headless ApiTests.
