# Pie chart viewer — JS API surface

**JS API source:** `public/js-api/src/viewer.ts:279` (`DG.Viewer.pieChart`),
`public/js-api/src/viewer.ts:696` (`PieChartViewer` class),
`public/js-api/src/interfaces/d4.ts:2855` (`IPieChartSettings`).


## What we are testing

The Pie chart specific JS-API surface. The typed factory
`DG.Viewer.pieChart(table, options)`, the `PieChartViewer` wrapper class, the
`onSegmentClicked` rxjs observable (Pie-only event), the typed
`IPieChartSettings` round-trip via `setOptions`/`getOptions`, the Pie-specific
properties metadata (`pieSortType`, `pieSortOrder`, `segmentAngleColumnName`,
`segmentLengthColumnName`, `includeNulls`), and the lifecycle when attached to
a `TableView`.

Pie has **no** `df.plot.pie(...)` helper (DataFramePlotHelper exposes scatter,
grid, tile, form, histogram, bar, heatMap, box, line, network — pie is absent),
so this scenario covers the `DG.Viewer.pieChart` factory and the named
`view.addViewer(VIEWER.PIE_CHART, ...)` path instead.

## Why it is uncovered

Coverage scan across `public/packages/ApiTests/src/` (excluding `src/ai/`):
zero references to `Viewer.pieChart`, `PieChartViewer`, `IPieChartSettings`,
or the `d4-pie-chart-on-segment-clicked` event id. Only
`src/dapi/layouts.ts:151` references the string `"PieChartLook"` while
asserting a layout JSON shape — a different surface. The Dart parametric
viewer test in `core/client/xamgle/lib/src/tests/viewers_test.dart` already
round-trips every `@Prop`, so this scenario stays away from `look[prop] = …`
round-trips and instead exercises Pie-specific JS-only surface — the typed
factory, the typed wrapper, the `onSegmentClicked` event, the property-bag
JSON shape for Pie's distinctive props, and the lifecycle.

The round-5 Histogram test covered the generic `getOptions`/`setOptions`/
`dataFrame`-swap surface on Histogram. We do not duplicate those generic
checks; we focus on Pie-specific differences instead.

## Preconditions

- A running Datagrok instance with the JS API loaded.
- `grok.data.demo.demog(<rowCount>)` available.
- No external server-side state created or modified.

## Test cases

1. **factory typed** — `DG.Viewer.pieChart(df, {category: 'race'})` returns a
   `DG.PieChartViewer` instance, `viewer.type === DG.VIEWER.PIE_CHART`, and
   the seeded `category` alias is reflected as `categoryColumnName` in
   `props` and in `getOptions(true).look` (mirrors BarChart's
   `value`→`valueColumnName` aliasing).

2. **pie-specific look round-trip via setOptions** — calling `setOptions({
   pieSortType: 'by category', pieSortOrder: 'desc', includeNulls: false,
   labelPosition: 'Outside' })` is observable on both `viewer.props` and on
   `getOptions(true).look`. These four keys are Pie-only and worth pinning.

3. **segment angle / length columns** — `setOptions({segmentAngleColumnName:
   'age', segmentAngleAggrType: 'avg', segmentLengthColumnName: 'weight',
   segmentLengthAggrType: 'max'})` is reflected in `props` and `look`.
   Confirms that Pie's two segment-coding column families round-trip
   independently. Also clears them with `''` and confirms the empty value
   sticks (Pie's default for both is the empty string).

4. **getProperties exposes pie-specific properties** — `viewer.getProperties()`
   returns a non-empty `DG.Property[]` that contains entries for
   `categoryColumnName`, `pieSortType`, `pieSortOrder`,
   `segmentAngleColumnName`, `segmentLengthColumnName`, and `includeNulls`.
   The `pieSortType` property's `choices` array contains both `'by value'`
   and `'by category'`.

5. **onSegmentClicked is an rxjs Observable** — `viewer.onSegmentClicked` is
   non-null, has a `subscribe` method, and the returned `Subscription`
   exposes `unsubscribe`. We do not synthesize a click — verifying the
   observable contract is the JS-only testable surface; the Dart-side
   behavior (firing on canvas click) is exercised by the Dart event tests.

6. **view.addViewer attaches a typed PieChartViewer** — using
   `tv.addViewer(DG.VIEWER.PIE_CHART, {category: 'race'})`, the returned
   value is a `Viewer` whose `type === 'Pie chart'` AND is an instance of
   `DG.PieChartViewer` (the typed JS wrapper, not just a generic Viewer).
   The viewer also appears in `tv.viewers` (cloned per read), and the entry
   pulled out of that array is itself a `DG.PieChartViewer`. `tv.close()`
   tears down shell state in `finally`.

Negative case observed during round 5 (and inherited here): calling
`close()` on a viewer created via `DG.Viewer.pieChart(df, ...)` but never
attached to a `TableView` throws inside the Dart interop. For this reason
the detached viewers used in cases 1-5 are simply allowed to go out of
scope rather than `close()`-ed. This asymmetry is platform-wide for
viewers and was already documented in the round-5 Histogram scenario.

Cases dropped on critic review (rationale): a `dataFrame` getter/setter
swap and a "close on attached viewer" lifecycle case were both removed —
both are Histogram-equivalent and already exercised by the round-5
Histogram scenario, so duplicating them on Pie adds no Pie-specific
coverage.

Other negative cases skipped: the JS API has no defined failure mode for
unknown look keys (Dart silently ignores) and no JS-observable failure
path for invalid sort type / aggregation strings (rendering shows an
"error" banner but `getOptions` round-trips silently).

## Out of scope

- Pixel-level rendering of pie segments, labels, legend (Dart-side).
- Selection / filter behavior driven by `onClick: 'Filter'` — covered by
  the manual Test Track scenarios in
  `public/packages/UsageAnalysis/files/TestTrack/Viewers/pie-chart.md`.
- `@Prop` round-trips that the Dart parametric viewer test already covers
  (e.g. `legendPosition`, color int props, margins).
- Trellis layout (`PieChartDescriptor.createTrellisViewer`).
- The `df.plot.pie` shorthand — it does not exist on `DataFramePlotHelper`.
