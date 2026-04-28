# BarChart extras

**JS API source:** `public/js-api/src/viewer.ts:682` (`BarChartViewer`), `public/js-api/src/interfaces/d4.ts:786` (`IBarChartSettings`)

## What we are testing

A second pass at `DG.BarChartViewer` and `IBarChartSettings` that pins down behaviour the round-7 BarChart suite, the `src/grid/viewer-set-property.ts` parametric round-trips, and the Dart-side targeted BarChart tests do not already cover. Specifically: the friendly-key `color` alias on the typed factory, the JSON-envelope shape when `splitColumnName` and `stackColumnName` are set together, the `barSortType` + `barSortOrder` pair round-trip, choices introspection on the Bar-specific `barSortType` enum, `includeNulls` boolean round-trip via `setOptions`, and the `onClick` `RowGroupAction` enum value persisting as a string in the look. All are checked through the public JS API only (`getOptions(true).look`, `props`, `getProperties()`).

## Why it is uncovered

Round 7 (`bar-chart-js-api.ts`) covers `valueAggrType`, `value`/`split`/`stack` factory aliases, `resetView`, and the three event Observables. `src/grid/viewer-set-property.ts` covers six simple props (`valueAggrType`, `showValueAxis`, `maxCategoryWidth`, `linearColorScheme`, `barSortOrder`, `relativeValues`) — but never `barSortType`, `includeNulls`, `onClick`, or the `color` factory alias. Round 6 PieChart pinned `pieSortType` choices; round 10 Histogram pinned `colorAggrType` choices — neither touches Bar's `barSortType` choices. Dart-side `viewers_test.dart` covers `splitColumnName + splitMap` together but not `splitColumnName + stackColumnName` in the JS look envelope.

## Preconditions

- `grok.data.demo.demog(N)` is available (no server-side state created).
- All cases create a detached viewer via `DG.Viewer.barChart(...)` or `df.plot.bar(...)`. No table view is opened, so no canvas geometry is touched.

## Test cases

1. **Color friendly-key alias on factory** — `DG.Viewer.barChart(df, {value: 'age', split: 'race', color: 'height'})` resolves the `color` short key to `colorColumnName` in both `props` and `getOptions(true).look`.
2. **Split + stack JSON-envelope shape** — set both `splitColumnName: 'race'` and `stackColumnName: 'sex'` via the factory; both must appear in `props` AND in `getOptions(true).look` together.
3. **`barSortType` + `barSortOrder` pair round-trip** — set `{barSortType: 'by value', barSortOrder: 'asc'}` via `setOptions`; both must persist together in `look` and `props`.
4. **`getProperties` choices on `barSortType`** — `getProperties()` must surface the `barSortType` property with a non-empty `choices` list containing the two known sort modes (`by category`, `by value`).
5. **`includeNulls` boolean round-trip** — toggle `includeNulls` to `false` via `setOptions`; the value must persist as a boolean in both `look` and `props`. (Default is `true`; this also exercises the false→true defaulting.)
6. **`onClick` `RowGroupAction` enum round-trip** — set `onClick: 'Filter'` via `setOptions`; the JSON envelope must store the string `'Filter'` in `look` and the wrapper must surface it via `props`. Repeat with `'Select'` to confirm the channel is not pinning a single value.

Negative case is intentionally skipped: setting unknown enum strings on these props is a "validate-on-render" path with no defined synchronous failure mode in the JS API. Already pinned on Histogram (round 10) for `colorAggrType` introspection; not duplicated here.

## Out of scope

- Rendering / canvas geometry / hit-testing — never-painted detached viewer, by design (see round 8/9 first-paint pitfalls).
- Property-grid editor UI — Dart-side concern.
- Layout persistence across `saveLayout`/`applyLayout` — already covered Dart-side in `viewers_test.dart`.
- Color-coding semantics (which rows feed which color bin) — covered by Dart unit tests; we only pin the JSON envelope.
