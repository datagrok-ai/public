# Bar chart viewer — JS API surface

**JS API source:** `public/js-api/src/viewer.ts:223` (`DG.Viewer.barChart`),
`public/js-api/src/viewer.ts:682` (`BarChartViewer` class — `resetView`,
`onCategoryClicked`, `onCategoryHovered`, `onResetView`),
`public/js-api/src/dataframe/data-frame.ts:566` (`DataFramePlotHelper.bar`),
`public/js-api/src/interfaces/d4.ts:786` (`IBarChartSettings`).


## What we are testing

The Bar-chart-specific JS-API surface. Specifically:

- The typed factory `DG.Viewer.barChart(table, options)` and the
  `df.plot.bar({...})` shorthand on `DataFramePlotHelper`.
- The friendly→canonical key aliasing seeded by the `IBarChartSettings`
  interface — `value` → `valueColumnName`, `split` → `splitColumnName`,
  `stack` → `stackColumnName`, `color` → `colorColumnName`.
- The `BarChartViewer.resetView()` JS-only wrapper (calls
  `grok_BarChartViewer_ResetView`) — must not throw on a freshly-attached
  Bar chart and must trigger the `d4-bar-chart-reset-view` event.
- The three Bar-only rxjs Observables: `onCategoryClicked`,
  `onCategoryHovered`, `onResetView`. The first two are verified for
  Observable shape only (they fire on real mouse interactions which we
  cannot synthesize from a Puppeteer-backed integration test). `onResetView`
  is exercised end-to-end: subscribe → call `viewer.resetView()` → assert
  the event fires within a small timeout.

## Why it is uncovered

Coverage scan across `public/packages/ApiTests/src/` (excluding `src/ai/`):
the only existing BarChart usage is six simple property-setter tests in
`src/grid/viewer-set-property.ts` (`valueAggrType`, `showValueAxis`,
`maxCategoryWidth`, `linearColorScheme`, `barSortOrder`, `relativeValues`)
plus a single `onContextMenu` subscription in `src/shell/events.ts`. None
of `BarChartViewer.resetView`, `onCategoryClicked`, `onCategoryHovered`,
or `onResetView` is referenced anywhere outside `src/ai/`.

The Dart parametric viewer test in
`core/client/xamgle/lib/src/tests/viewers_test.dart` round-trips every
`@Prop` and the targeted Dart Bar-chart tests cover `splitColumnName`,
`aggrType`, `aggrColumnName`, `splitMap`, and legend persistence. The
JS-only events and the `resetView()` JS wrapper are not exercised by any
of those.

The round-5 Histogram and round-6 PieChart scenarios covered the generic
`getOptions`/`setOptions`/`dataFrame`-swap surface and the Observable
shape contract. We do not duplicate those generic checks here; instead we
focus on Bar-specific JS-only surface — the friendly-key factory aliasing,
the `resetView()` wrapper plus the matching `onResetView` round-trip, and
the Observable shape of the two click/hover events.

## Preconditions

- A running Datagrok instance with the JS API loaded.
- `grok.data.demo.demog(<rowCount>)` available.
- No external server-side state created or modified.

## Test cases

1. **factory typed via DG.Viewer.barChart** —
   `DG.Viewer.barChart(df, {value: 'age', valueAggrType: 'avg', split: 'race'})`
   returns a `Viewer` whose `type === DG.VIEWER.BAR_CHART`. The seeded
   friendly aliases (`value`, `split`) are reflected as their canonical
   `*ColumnName` forms in `viewer.props` and in `getOptions(true).look`.
   Note: unlike `DG.Viewer.histogram`/`DG.Viewer.pieChart`, the factory is
   typed `Viewer<IBarChartSettings>` (it goes through `Viewer.fromType`),
   but the runtime instance is still a `DG.BarChartViewer` via `toJs`.

2. **factory via df.plot.bar shorthand** — `df.plot.bar({value: 'height',
   split: 'race', stack: 'sex'})` returns a viewer whose `dataFrame === df`,
   `type === DG.VIEWER.BAR_CHART`, and whose seeded keys are reflected as
   `valueColumnName`, `splitColumnName`, `stackColumnName`. Mirrors the
   round-5 Histogram `df.plot.histogram` case.

3. **resetView() does not throw on attached BarChart** — attach via
   `tv.addViewer(DG.VIEWER.BAR_CHART, {value: 'age', split: 'race'})`, cast
   the result to `DG.BarChartViewer`, call `resetView()` once. The call
   returns void and does not throw. (`resetView` is a thin wrapper around
   `grok_BarChartViewer_ResetView` — the JS-observable contract is
   "completes synchronously without throwing".)

4. **onCategoryClicked and onCategoryHovered are rxjs Observables** —
   both event getters are non-null, expose `subscribe`, and the returned
   `Subscription` has `unsubscribe`. Shape check only; firing requires
   real canvas mouse events which the Dart-side viewer event tests cover.

5. **onResetView round-trip via resetView()** — attach a BarChart, subscribe
   to `viewer.onResetView`, call `viewer.resetView()`, assert the event
   fires within a small timeout (1500 ms). Uses `take(1)` + a
   `Promise.race` against a `DG.delay`-backed timeout so a missed fire
   fails clearly instead of hanging the test runner.

## Out of scope

- Pixel-level rendering of bars, labels, legend (Dart-side).
- Selection / filter behavior driven by `onClick: 'Filter'` — covered by
  the manual TestTrack scenarios in
  `public/packages/UsageAnalysis/files/TestTrack/Viewers/bar-chart-tests.md`.
- `@Prop` round-trips that the Dart parametric viewer test already covers
  (legend position, color int props, margins, axis font, etc.).
- The simple Bar-chart property setters already covered in
  `src/grid/viewer-set-property.ts` (`valueAggrType`, `showValueAxis`,
  `maxCategoryWidth`, `linearColorScheme`, `barSortOrder`, `relativeValues`).
- The `splitColumnName` / `splitMap` / `aggrType` / `aggrColumnName` surface
  already covered by targeted Dart Bar-chart tests.
- Synthesized canvas mouse events for `onCategoryClicked` and
  `onCategoryHovered` — the firing path is Dart-side and is exercised by
  the Dart event tests; the JS-side surface is the Observable contract.

Negative case observed in rounds 5–6 (and inherited here): calling
`close()` on a viewer created via `DG.Viewer.barChart(df, ...)` but never
attached to a `TableView` throws inside the Dart interop. Detached viewers
in cases 1–2 / 4 are simply allowed to go out of scope rather than
`close()`-ed; cases 3 and 5 attach via `tv.addViewer(...)` and clean up
with `tv.close()` in `finally`.
