# PC Plot (Parallel Coordinates) JS API

**JS API source:** `public/js-api/src/viewer.ts:275` (`DG.Viewer.pcPlot`), `public/js-api/src/viewer.ts:704` (`DG.PcPlot`), `public/js-api/src/interfaces/d4.ts:2692` (`IPcPlotSettings`)

## What we are testing

The Parallel Coordinates (PC) plot viewer surfaced through the JS API: the
factory `DG.Viewer.pcPlot(df, options)`, the typed wrapper class `DG.PcPlot`
(note: not suffixed with `Viewer` like most others), the defining
`columnNames: Array<string>` setting that picks which columns become axes,
PC-specific Boolean toggles (`normalizeEachColumn`, `showCurrentLine`,
`showAllLines`), the PC-only event observables `onLineClicked` and
`onLineHovered`, and that `view.addViewer(VIEWER.PC_PLOT, ...)` returns a
`DG.PcPlot` instance.

## Why it is uncovered

There are no PC-plot-specific tests anywhere under
`public/packages/ApiTests/src/` (verified by greping for `PC_PLOT`,
`pcPlot`, `PcPlot`). The only existing coverage is the generic Dart
parametric prop iteration in `core/client/xamgle/lib/src/tests/viewers_test.dart`,
which round-trips primitive `@Prop` values across all viewers but does not
assert PC's `Array<string> columnNames` round-trip, the typed `DG.PcPlot`
wrapper, the factory shorthand, or the PC-only line events.

## Preconditions

- A running Datagrok server (no extra entities required).
- `grok.data.demo.demog()` available — used as the data source. Demog has
  numeric `age`, `height`, `weight` columns suitable for PC axes.
- No `df.plot.pcPlot` shorthand exists on `DataFrameViewerHelper` (verified
  in `public/js-api/src/dataframe/data-frame.ts`); we only test
  `DG.Viewer.pcPlot` and `view.addViewer(VIEWER.PC_PLOT, ...)`.

## Test cases

1. **factory `DG.Viewer.pcPlot` returns typed `DG.PcPlot` with PC settings** — call
   `DG.Viewer.pcPlot(df, {columnNames: ['age', 'height', 'weight']})`; assert
   `v instanceof DG.PcPlot`, `v.type === DG.VIEWER.PC_PLOT`,
   `v.dataFrame === df`, and `v.props['columnNames']` is an array containing
   each of the three names.
2. **`columnNames` array round-trip via `getOptions(true).look` and `setOptions`** —
   set `columnNames` to `['age', 'height']` initially, assert array shape and
   contents through both `props['columnNames']` and `getOptions(true).look.columnNames`;
   then call `setOptions({columnNames: ['age', 'weight']})` and assert the new
   array round-trips through both surfaces.
3. **PC-specific Boolean props round-trip via `setOptions`** — flip
   `normalizeEachColumn`, `showCurrentLine`, and `showAllLines` from their
   defaults via `setOptions`, assert each is reflected in `v.props[...]` and
   in `getOptions(true).look`. Also assert `getProperties()` introspection
   surfaces a `Property` for `columnNames` whose name matches.
4. **`onLineClicked` and `onLineHovered` are subscribable rxjs Observables** —
   construct via the factory, assert both event getters return objects with a
   `subscribe` function and that the returned `Subscription` has `unsubscribe`.
   Pure shape check; no canvas geometry, no first-paint pitfalls.
5. **`view.addViewer(VIEWER.PC_PLOT, ...)` attaches a typed `DG.PcPlot`** — open
   a TableView, attach via `view.addViewer`, assert returned viewer is a
   `DG.PcPlot`, walk `tv.viewers` and find the PC plot, confirm that one is
   also a `DG.PcPlot` instance. Close the TableView in `finally`.

## Out of scope

- Dart parametric prop round-trips (already pinned by `viewers_test.dart`).
- Canvas geometry (`viewBox`, axes positions, screen-to-world) — known
  first-paint pitfall in this pipeline (rounds 8/9), and PC plot has no
  `viewBox` getter on its wrapper anyway.
- In-chart range-slider filtering and Pick up/Apply — UI flows owned by the
  Playwright TestTrack scenario at
  `public/packages/UsageAnalysis/files/TestTrack/Viewers/pc-plot.md`.
- Color encoding rendering, density styles, and legend positioning — visual
  concerns out of scope for an API integration test.
