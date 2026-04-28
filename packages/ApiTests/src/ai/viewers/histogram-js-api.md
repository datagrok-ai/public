# Histogram viewer — JS API surface

**JS API source:** `public/js-api/src/viewer.ts:219` (`DG.Viewer.histogram`),
`public/js-api/src/viewer.ts:671` (`HistogramViewer` class),
`public/js-api/src/dataframe/data-frame.ts:565` (`df.plot.histogram`).


## What we are testing

The JS-API surface of the Histogram viewer. Specifically: the typed factory
`DG.Viewer.histogram(table, options)` and its sibling `df.plot.histogram(options)`,
the inherited `Viewer.getOptions()`/`setOptions()` JSON round-trip, the
`getInfo()` / `getProperties()` / `type` introspection helpers, the `dataFrame`
getter/setter, the `props` ObjectPropertyBag, the `close()` cleanup path, and
the four custom event streams declared on `HistogramViewer`
(`onBinsSelected`, `onLineSelected`, `onMouseOverBins`, `onMouseOverLine`).

## Why it is uncovered

Coverage scan across `public/packages/ApiTests/src/` (excluding `src/ai/`) finds
zero references to `Viewer.histogram`, `HistogramViewer`, `plot.histogram`, or
the `d4-histogram-*` event ids. Only `Stats.histogramsByCategories` is touched
by `src/stats/stats.ts:88` — a different API. The Dart parametric tests in
`core/client/xamgle/lib/src/tests/viewers_test.dart` already round-trip every
`@Prop` on every viewer, so this scenario deliberately stays away from
`look[prop] = …` round-trips and instead exercises the JS-only surface
(factories, JSON shapes, event observable shape, lifecycle).

## Preconditions

- A running Datagrok instance with the JS API loaded.
- `grok.data.demo.demog(<rowCount>)` available (used by all hand-written
  ApiTests).
- No external server-side state is created or modified.

## Test cases

1. **factory typed** — `DG.Viewer.histogram(df, {value: 'age', bins: 15})`
   returns an instance of `DG.HistogramViewer`, has `type === DG.VIEWER.HISTOGRAM`,
   and the seeded options are reflected in `props` and in
   `getOptions(true).look`.

2. **factory via DataFrame.plot.histogram** — `df.plot.histogram(...)` returns
   a `DG.Viewer` whose `type === DG.VIEWER.HISTOGRAM` and whose `dataFrame` is
   the same dataframe that produced it. Confirms the deprecated-but-still-public
   helper still wires through.

3. **getOptions excludes defaults by default** — after constructing with no
   options, `getOptions().look` only contains keys explicitly set or marked
   non-default; `getOptions(true).look` contains many more (>=10) keys
   including `bins`. Verifies the `includeDefaults` flag actually changes the
   shape.

4. **setOptions round-trip** — `setOptions({bins: 7, showXAxis: true,
   splitColumnName: 'race'})` is observable on both `viewer.props` and
   `getOptions(true).look`, and the JSON returned by `getOptions(true)`
   carries `id` and `type === 'Histogram'`.

5. **dataFrame getter/setter swap** — building a histogram on `demog`, then
   assigning a freshly-created small DataFrame with one numerical column,
   `viewer.dataFrame` reads back the new frame. Confirms the
   `grok_Viewer_Set_DataFrame` interop path on Histogram.

6. **getInfo / getProperties / descriptor shape** — `getInfo()` returns a
   non-null object, `getProperties()` returns a non-empty array of
   `DG.Property` instances that includes a `bins` property, and
   `descriptor.name === 'Histogram'`. Light shape checks only — exact key
   sets are not part of the public contract.

7. **event streams are rxjs Observables** — each of `onBinsSelected`,
   `onLineSelected`, `onMouseOverBins`, `onMouseOverLine` returns an
   rxjs `Observable` whose `subscribe(...)` returns a `Subscription` with an
   `unsubscribe()`. We do not synthesize histogram clicks in JS — verifying
   the observable contract (subscribe/unsubscribe without throwing) is the
   testable JS-only surface.

8. **close detaches without throwing on attached viewer** — attach a
   histogram to a real `TableView` (`tv.addViewer(VIEWER.HISTOGRAM, ...)`),
   then `viewer.close()` completes without throwing. `TableView` is closed
   in `finally` to clean up shell state.

Negative case observed during round 5: calling `close()` on a viewer that
was created via `DG.Viewer.histogram(df, ...)` but never attached to a
`TableView` throws `Cannot read properties of null (reading 'C9')` from
inside `grok_Viewer_Close` (Dart `_initJsApi` interop). For this reason
the other seven cases above do **not** call `close()` on the detached
viewer they create — they simply let it go out of scope, mirroring the
pattern in `src/grid/viewer-set-property.ts`. This asymmetry is a real
JS-observable behaviour of the platform; surfacing it via the run is the
value-add of this scenario file.

Other negative cases skipped: the API has no defined failure mode for
invalid property names (Dart silently ignores unknown look keys), and
there is no JS-observable failure path for setting `bins` to 0 or
splitting on a non-existent column — both round-trip silently. Documented
here per the spec template.

## Out of scope

- Pixel-level rendering of bins, axes, or legends (Dart-side; not JS-API).
- Filter-panel integration (`FilterGroup.updateOrAdd({type: HISTOGRAM, …})`)
  — covered indirectly by `src/grid/filterGroup.ts` for the bar/grid
  pathway; a dedicated histogram-filter test would belong in
  `src/ai/filters/...`, not here.
- `@Prop` round-trips that the Dart parametric viewer test already covers.
- Trellis layout (`HistogramDescriptor.createTrellisViewer`) — requires
  Trellis JS bindings that are out of scope for this single-viewer file.
