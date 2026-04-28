# BoxPlot JS API

**JS API source:** `public/js-api/src/viewer.ts:231` (`DG.Viewer.boxPlot`), `public/js-api/src/viewer.ts:713` (`DG.BoxPlot`), `public/js-api/src/interfaces/d4.ts:999` (`IBoxPlotSettings`), `public/js-api/src/dataframe/data-frame.ts:568` (`df.plot.box`)

## What we are testing

The Box Plot viewer surfaced through the JS API: the typed factory
`DG.Viewer.boxPlot(df, options)` returning `DG.BoxPlot` (note: bare `BoxPlot`,
no `Viewer` suffix — same convention as `DG.PcPlot` from round 12), the
shorthand `df.plot.box(...)` returning a base `Viewer` whose `type` still
matches `DG.VIEWER.BOX_PLOT`, friendly-key aliasing on the factory options
(`value`/`category1`/`category2` collapsing to `valueColumnName` /
`category1ColumnName` / `category2ColumnName`), the multi-key `category1` +
`category2` JSON envelope (uncovered combination), the BoxPlot-specific
`statisticsFormat` choices via `getProperties()` introspection, the BoxPlot-only
event observables (`onResetView`, `onAfterDrawScene`, `onBeforeDrawScene`,
`onPointClicked`), and the typed-instance check via
`view.addViewer(VIEWER.BOX_PLOT, ...)`.

## Why it is uncovered

There are no BoxPlot-specific tests anywhere under
`public/packages/ApiTests/src/` (verified by greping for `boxPlot`, `BoxPlot`,
`BOX_PLOT`). The Dart-side targeted feature test at
`core/client/d4/lib/src/viewers/box_plot/features/box_plot_t_test.dart` covers
Welch/Alexander-Govern p-value rendering and the `T`-key shortcut, and
`core/client/xamgle/lib/src/tests/viewers_test.dart` parametric-prop iteration
covers individual primitive `@Prop` round-trips (including
`categoryColumnName`, `categoryColumnName2`, `valueColumnName`, `showPValue`,
`statisticsFormat`). What remains uncovered: the typed JS wrapper `DG.BoxPlot`
factory shape, friendly-key aliasing, the *combined* `category1` + `category2`
JSON envelope, BoxPlot's four event observables, and the
`view.addViewer(VIEWER.BOX_PLOT)` typed-instance contract.

## Preconditions

- A running Datagrok server (no extra entities required).
- `grok.data.demo.demog()` available — its `age`, `sex`, `race` columns are
  suitable for value/category combinations.
- Lifecycle: detached factory-built viewers are never `close()`-d (round 12
  pinning); attached viewers are `close()`-d via `tv.close()` in `finally`.

## Test cases

1. **factory `DG.Viewer.boxPlot` returns typed `DG.BoxPlot`** — call
   `DG.Viewer.boxPlot(df, {valueColumnName: 'age', category1ColumnName: 'race'})`;
   assert `v instanceof DG.BoxPlot`, `v.type === DG.VIEWER.BOX_PLOT`,
   `v.dataFrame === df`, and `v.props['valueColumnName'] === 'age'` plus
   `v.props['category1ColumnName'] === 'race'`. Also assert that
   `df.plot.box(...)` shorthand returns a base `Viewer` (not necessarily
   `DG.BoxPlot`) whose `type` is `DG.VIEWER.BOX_PLOT`.
2. **friendly-key aliases `value`/`category1`/`category2` on factory options** —
   call `DG.Viewer.boxPlot(df, {value: 'age', category1: 'race', category2: 'sex'})`;
   assert that the friendly aliases collapse to canonical names: `props` and
   `getOptions(true).look` both surface `valueColumnName === 'age'`,
   `category1ColumnName === 'race'`, `category2ColumnName === 'sex'`.
3. **combined `category1ColumnName` + `category2ColumnName` JSON envelope round-trip** —
   construct with one category, then `setOptions({category1ColumnName: 'race',
   category2ColumnName: 'sex'})`; assert both are reflected together in
   `getOptions(true).look` (the two-key combined JSON envelope is what
   distinguishes this from the Dart parametric harness, which only round-trips
   each separately).
4. **`showStatistics` + `statisticsFormat` combined round-trip with `getProperties` choices** —
   set `{showStatistics: true, statisticsFormat: 'two digits after comma'}`;
   assert both reflected in `props` and `getOptions(true).look`. Then walk
   `getProperties()` and assert that the `statisticsFormat` `Property` exists
   and has a non-empty `choices` array including `'two digits after comma'`.
   (The author's first guess was a Unicode-bearing label `'count, mean ± stdev'`;
   the actual `BoxPlotLook.statisticsFormatChoices` is a numeric-format list:
   `'int'`, `'one/two/three/four digits after comma'`, `'scientific'`,
   `'2/3/4 significant digits'`, `'full precision'`.)
5. **BoxPlot-specific events are subscribable rxjs Observables** — construct via
   the factory (detached); assert each of `onResetView`, `onAfterDrawScene`,
   `onBeforeDrawScene`, `onPointClicked` returns an object with a `subscribe`
   function and that the returned `Subscription` has `unsubscribe`. Pure shape
   check, no canvas geometry, no first-paint pitfalls (rounds 8/9 lessons).
6. **`view.addViewer(VIEWER.BOX_PLOT, ...)` attaches a typed `DG.BoxPlot`** — open
   a TableView, attach via `tv.addViewer`, assert returned viewer is a
   `DG.BoxPlot`, walk `tv.viewers` and find the box plot, confirm that one is
   also a `DG.BoxPlot` instance. Close the TableView in `finally`.

## Out of scope

- Dart parametric prop round-trips (already pinned by `viewers_test.dart`).
- Welch's t-test / Alexander-Govern p-value rendering, `T`-key toggle,
  click-to-Wikipedia (already covered by the Dart targeted feature test
  at `core/client/d4/lib/src/viewers/box_plot/features/box_plot_t_test.dart` —
  off-limits, and heavy DOM/keyboard interaction is not feasible in headless
  ApiTests anyway).
- `resetView()` JS method — calling it triggers a canvas redraw on a
  freshly-attached viewer; rounds 8/9 first-paint pitfalls warn against this.
- Canvas geometry (`viewport: Rect`) — same first-paint pitfall.
- Violin plot KDE rendering, layout save/restore, project zoom preservation —
  Playwright UI scenarios at
  `public/packages/UsageAnalysis/files/TestTrack/Viewers/box-plot.md`.
