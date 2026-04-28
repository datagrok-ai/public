# Histogram viewer — extras (round 10)

**JS API source:**
- `public/js-api/src/viewer.ts:219` (`DG.Viewer.histogram`)
- `public/js-api/src/viewer.ts:671` (`HistogramViewer` class)
- `public/js-api/src/interfaces/d4.ts:1733` (`IHistogramSettings`)
- `public/js-api/src/entities/property.ts:121` (`Property` class — `choices`, `columnTypeFilter`)


## What we are testing

Histogram-specific JS-API surface that round 5 (`histogram-js-api.{md,ts}`) did
not exercise and that the Dart parametric viewer harness in
`core/client/xamgle/lib/src/tests/viewers_test.dart` (lines 3–99) cannot reach
through its boolean-toggle / int-midpoint / column-iterator strategy. The
focus is on JSON shapes round-tripped through `setOptions` /
`getOptions(true).look` and on the contents of `getProperties()` exposed to JS
(the parametric Dart harness only iterates choices/min-max — it never reads
`Property.choices` or `Property.columnTypeFilter` back out and never sets
array-typed or formula-typed look values).

## Why it is uncovered

Round 5 (`histogram-js-api.ts`) covers: typed factory, `df.plot.histogram`,
`getOptions` defaults flag, a single `setOptions` round-trip on
`bins`/`showXAxis`/`splitColumnName`, dataFrame swap, descriptor/getInfo
shape, the four event observables, and `close()` on attached. It does
**not** touch:

- clearing a `*ColumnName` look field (assignment of empty string)
- formula-string look values such as `filter`
- array-typed look values such as `linearColorScheme`
- `Property.columnTypeFilter` introspection on the seven `*ColumnName` props
- `Property.choices` introspection on `colorAggrType`
- string-with-CSV-shape look values such as `aggTooltipColumns`

A repo-wide grep for `HistogramViewer` outside `src/ai/` returns zero hits, so
no other ApiTests file covers any of these.

## Preconditions

- A running Datagrok instance with the JS API loaded.
- `grok.data.demo.demog(...)` available (used by every existing histogram
  ApiTest) — provides numerical columns `age`, `height`, `weight` and
  categorical columns `race`, `sex`, `disease`.
- No external server-side state created or modified — every viewer is
  detached and discarded; no `addViewer` calls.

## Test cases

1. **clear splitColumnName via setOptions empty string** — set
   `splitColumnName: 'race'` then `splitColumnName: ''`. Assert that both
   `viewer.props['splitColumnName']` and `getOptions(true).look.splitColumnName`
   are equal to `''` (or `null`/`undefined` — accept any falsy round-trip,
   document the actual shape) after the clear. Asymmetric to round 5 which
   only sets a non-empty value.

2. **filter formula round-trip** — `setOptions({filter: '${age} > 20'})`
   leaves `getOptions(true).look.filter === '${age} > 20'` (string identity)
   and `viewer.props['filter']` matches. Confirms formula strings flow through
   the JSON envelope unmolested. Bin-count delta is **not** asserted (no
   first-paint signal on a detached viewer; round 8/9 lessons).

3. **linearColorScheme number-array round-trip** — set
   `linearColorScheme: [0xFF112233, 0xFFAABBCC]` via `setOptions`. Assert the
   look JSON returns an `Array.isArray(look.linearColorScheme) === true` of
   length 2 with element-wise equality on the integer values. The Dart
   parametric harness only handles BOOL / numeric-with-min-max / column /
   string-with-choices — array-typed properties are never round-tripped.

4. **getProperties columnTypeFilter introspection** — call `getProperties()`
   on a fresh histogram, find the entries whose name is `valueColumnName`,
   `splitColumnName`, `colorColumnName`. Assert each one's `propertyType ===
   DG.TYPE.COLUMN` (or `'column'`) and that `columnTypeFilter` is a non-null
   string ('numerical' for value/color, 'categorical' for split). Pure JS-side
   introspection — Dart parametric reads `prop.columnTypeFilter` only to
   *select* a column, never asserts its identity.

5. **getProperties choices introspection on colorAggrType** — find the
   property with name `colorAggrType` in `getProperties()`. Assert
   `Array.isArray(prop.choices)` is `true` and `prop.choices.length > 0`.
   We don't pin the exact set (it grows with new aggregations) — just that
   the JS getter returns a usable list. The Dart parametric harness only
   *iterates* choices; it never asserts the JS shape.

6. **aggTooltipColumns string round-trip** — `setOptions({aggTooltipColumns:
   'avg(age), count(race)'})`. Confirm the exact string echoes through
   `getOptions(true).look.aggTooltipColumns` and `viewer.props
   ['aggTooltipColumns']`. Acts as a regression check on string-typed,
   non-choice look properties (the Dart parametric harness skips strings
   without `choices`).

## Out of scope

- Any boolean-toggle, choice-iteration, or column-iteration round-trip — the
  Dart parametric harness already exercises every `@Prop` getter/setter on
  Histogram across those three axes (`viewers_test.dart:68-98`).
- Pixel-level rendering / first-paint geometry — round 8 / round 9 already
  proved this is unreliable on a detached viewer.
- Filter-panel integration (`FilterGroup.updateOrAdd({type: HISTOGRAM, ...})`)
  — covered by `src/grid/filterGroup.ts` for the underlying mechanism;
  histogram-specific filter routing is a separate feature area.
- `close()` on either attached or detached viewer — already the subject of
  round 5's last two test cases (and a known platform bug surfaced there).
- Trellis layout (`createTrellisViewer`) — out of scope for a single-viewer
  file.

## Negative cases

The histogram look-property API has no defined JS-observable failure mode
for any of the inputs above: unknown keys are silently dropped, malformed
formulas in `filter` are accepted at JSON-set time and only surface during
rendering (not exposed as a JS event), and arrays of the wrong length /
wrong type for `linearColorScheme` are silently coerced. Per the spec
template's negative-case skip clause, this scenario records the skip
rather than synthesising a meaningless failing case.
