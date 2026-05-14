---
name: datagrok-filtering
description: Filter rows of a Datagrok DataFrame inside a datagrok-exec block — by predicate, range, equals/contains/in-set, regex, substructure (SMILES / SMARTS / molblock auto-detect), or by combining masks. Also covers clearing, inverting, the show-only-filtered vs destructive-drop split, and the filter event lifecycle (onRowsFiltering / onFilterChanged / onRowsFiltered). Use whenever the user says "filter", "show only", "hide rows where", "find rows that", "clear the filter", "invert the filter", or asks for the filtered subset as a new table. Does NOT cover selection (separate skill) or generic DataFrame cloning (datagrok-df-and-columns).
---

# datagrok-filtering

Use the `grokky.*` filter helpers inside a `datagrok-exec` block. They wrap
`df.filter` (a `DG.BitSet`), the `FilterGroup` UI, and the Chem substructure
path so Claude doesn't have to remember that `setAll(true)` clears the filter,
that `df.filter.init` zeros the buffer before applying its predicate, that
`Chem:substructureFilter` reads `state.molBlock` (and silently matches zero
rows if you hand it raw SMILES), or that `df.rows.removeWhereIdx` deletes
rows whereas `df.filter` only hides them.

## What this skill covers

Row-shape filtering on `DG.DataFrame`: reading the current filter state,
writing a predicate-based filter, single-column range/equals/in-set/regex
filters via the UI filter group, substructure filtering on molecule columns,
combining masks (AND / OR / AND NOT), clearing and inverting the filter,
non-destructive "give me only the filtered rows" via `df.clone(df.filter)`,
destructive `removeWhereIdx`, and subscribing to the filter event lifecycle.

**Out of scope.** Selection state (`df.selection`) lives in the
`datagrok-selection` skill. Generic DataFrame cloning, column work, and
"give me the visible rows as a new DF without filters at all" live in
`datagrok-df-and-columns` (`grokky.cloneDf`).

## Quick reference

| Helper                                                | One-liner                                                                |
|-------------------------------------------------------|--------------------------------------------------------------------------|
| `grokky.filterRows(target, col, criteria, opts?)`     | Single-column range / values / substructure via `FilterGroup` UI.        |
| `grokky.filterByPredicate(df, pred)`                  | Keep rows where `pred(i)` is true. Buffer-direct, single notification.   |
| `grokky.filterSubstructure(target, col, query, opts?)`| SMILES / SMARTS / molblock — auto-detect, auto-switch above 50k rows.    |
| `grokky.combineFilters(df, op, ...preds)`             | AND / OR several predicates in a single pass.                            |
| `grokky.clearFilter(dfOrView)`                        | Show every row. Polymorphic: also disables UI filter group when given a view. |
| `grokky.invertFilter(df)`                             | Flip the current mask.                                                   |
| `grokky.filteredDf(df, opts?)`                        | Non-destructive — clone of currently-passing rows.                       |
| `grokky.dropRows(df, pred)`                           | **Destructive** — removes rows where `pred(i)` is true.                  |

Globals available inside every `datagrok-exec` block: `grok`, `ui`, `DG`,
`view`, `t` (the current `DG.DataFrame`, when the view is a TableView),
`grokky`.

## The mental model

`df.filter` is a `DG.BitSet` — one bit per row in the DataFrame. The
**polarity is the #1 footgun in this skill**:

| Bit value | Meaning                              |
|-----------|--------------------------------------|
| `true`    | row **passes** the filter (visible)  |
| `false`   | row is **filtered out** (hidden)     |

So "clear the filter" (show every row) is `df.filter.setAll(true)`. Reaching
for `setAll(false)` hides everything. Always `setAll(true)` — or just call
`grokky.clearFilter(df)`.

Two other consequences of the BitSet model:

- All bitwise mutators (`and`, `or`, `xor`, `andNot`, `invert`, `setAll`,
  `init`, `copyFrom`) are **in-place** and return `this`. They do **not**
  return a fresh BitSet.
- `df.filter.init(predicate)` **zeros the buffer first** before applying
  `predicate`. To intersect with the existing filter, build a fresh BitSet
  via `DG.BitSet.create(df.rowCount, i => ...)` and call
  `df.filter.and(thatBitSet)`. Don't reach for `init` to "add" a condition.

`df.selection` is also a `DG.BitSet` with the same conventions — anything
you learn here transfers to the selection skill.

## Reading filter state

| Need                                  | Use                                          |
|---------------------------------------|----------------------------------------------|
| rows currently passing                | `df.filter.trueCount`                        |
| rows currently hidden                 | `df.filter.falseCount`                       |
| total rows                            | `df.filter.length` (same as `df.rowCount`)   |
| indices of passing rows               | `df.filter.getSelectedIndexes()` (Int32Array)|
| is row `i` passing?                   | `df.filter.get(i)`                           |
| any rows passing?                     | `df.filter.anyTrue`                          |
| any rows hidden?                      | `df.filter.anyFalse`                         |

```datagrok-exec
// How many rows pass the current filter?
return {visible: t.filter.trueCount, hidden: t.filter.falseCount, total: t.rowCount};
```

```datagrok-exec
// Indices of the currently-visible rows. Fast — `getSelectedIndexes` is cached.
return t.filter.getSelectedIndexes();
```

`getSelectedIndexes()` caches its result. If you mutate via a direct buffer
write (e.g. `getBuffer()`), the cache goes stale — call `fireChanged()` to
invalidate.

## Filtering by predicate

`grokky.filterByPredicate(df, pred)` is the canonical fast path. It wraps
`df.filter.init(pred)`: buffer-direct, single notification, ~10× faster than
`df.rows.filter(row => ...)`.

**Polarity (memorize):** `pred(i)` returns `true` to **keep** row `i`. This is
the opposite of `dropRows` ("returns `true` to remove") — see "Show-only vs
drop-rows" below.

```datagrok-exec
// Keep rows where age > 40.
grokky.filterByPredicate(t, (i) => t.getCol('age').get(i) > 40);
return {visible: t.filter.trueCount};
```

```datagrok-exec
// Multi-column predicate in a single pass — strictly preferred over building
// two separate BitSets and AND-ing them.
const mw = t.getCol('MW');
const logP = t.getCol('cLogP');
grokky.filterByPredicate(t, (i) => mw.get(i) < 500 && logP.get(i) < 5);
```

`filterByPredicate` **replaces** the current filter — it does not intersect.
If you want to narrow on top of the existing filter:

```datagrok-exec
// Intersect with existing filter: build a fresh mask, then AND.
const col = t.getCol('activity');
const mask = DG.BitSet.create(t.rowCount, (i) => col.get(i) > 7);
t.filter.and(mask);
```

Avoid `df.rows.filter(row => ...)` for performance-sensitive work — it goes
through `RowList` and calls `bitset.set` per row.

## Filtering by column value (range, equals, in-set, regex)

`grokky.filterRows(target, columnName, criteria, opts?)` handles the common
single-column shapes via the UI `FilterGroup`. The filter widget appears in
the side panel and persists in saved layouts. `target` may be a `DG.TableView`
(adds a UI filter) or a `DG.DataFrame` (writes a BitSet directly — no UI).

| Criteria shape                                  | Filter type                          |
|-------------------------------------------------|--------------------------------------|
| `{min, max}` (numeric column)                   | `DG.FILTER_TYPE.HISTOGRAM`           |
| `{values: [...]}` (categorical)                 | `DG.FILTER_TYPE.CATEGORICAL`         |
| `{contains: 'substr'}` (string column)          | `DG.FILTER_TYPE.FREE_TEXT`           |
| `{regex: '/pattern/'}` (string column)          | `DG.FILTER_TYPE.FREE_TEXT` (regex)   |
| `{substructure: '...'}` (molecule column)       | `DG.FILTER_TYPE.SUBSTRUCTURE`        |

```datagrok-exec
// Range filter on a numeric column. Shows up in the side panel as a histogram.
await grokky.filterRows(view, 'MW', {min: 200, max: 500});
```

```datagrok-exec
// Equals one value (use the single-element form of `values`).
await grokky.filterRows(view, 'category', {values: ['A']});
```

```datagrok-exec
// In-set: keep rows where category is one of these.
await grokky.filterRows(view, 'category', {values: ['A', 'B', 'C']});
```

```datagrok-exec
// Free-text contains. Case-insensitive substring on a string column.
await grokky.filterRows(view, 'name', {contains: 'acid'});
```

For multi-column conditions (Lipinski-style "MW < 500 AND cLogP < 5") drop
down to `filterByPredicate` — `filterRows` does one column at a time. The UI
filters then collaborate via `onRowsFiltering`: each filter contributes its
mask, the results are AND-ed. If you call `filterByPredicate` outside that
event, your write is overwritten on the next UI-filter cycle — see
"Combining filter conditions" below.

`opts.ui` controls UI vs BitSet path. Default is the UI path when `target` is
a `TableView`. For very large tables (`df.rowCount > 50_000`) the substructure
widget is laggy; the helper auto-switches to BitSet for substructure
criteria. Pass `opts.ui: true` to force the UI path even then; `opts.ui:
false` to force BitSet.

## Substructure filter

Auto-detects the query format:

- `DG.chem.isMolBlock(s)` — pure JS, cheap (`s.includes('M  END')`).
- `DG.chem.isSmarts(s)` — async-ish, calls `Chem:isSmarts` (requires Chem package).
- Otherwise → treat as SMILES.

Converts to molblock via `DG.chem.convert(s, fromNotation, toNotation)` before
calling `Chem:substructureFilter` (which reads `state.molBlock` and silently
matches zero rows for raw SMILES/SMARTS — this is the bug `filterRows`
previously had).

Above 50k rows, the helper auto-switches to the BitSet path
(`grok.chem.searchSubstructure(col, pattern)` → write onto `df.filter`).
No UI widget appears in that case. Override with `opts.ui: true` or
`opts.ui: false`.

```datagrok-exec
// SMILES input — most common.
await grokky.filterSubstructure(view, 'smiles', 'c1ccccc1');
```

```datagrok-exec
// SMARTS input — detected automatically via DG.chem.isSmarts.
await grokky.filterSubstructure(view, 'smiles', '[#6][!#1]');
```

```datagrok-exec
// Molblock input — detected via the "M  END" sentinel.
const molblock = `\n  Mrv...\n\n  6  6  0  0  ...\nM  END\n`;
await grokky.filterSubstructure(view, 'smiles', molblock);
```

If the Chem package isn't loaded, `DG.chem.isSmarts` and `DG.chem.convert`
throw. The helper catches and logs a `console.warn`, then falls through with
the raw string as SMILES — the `Chem:substructureFilter` widget will pick up
detection itself once it loads (`molBlockFailover` field).

For the no-UI fast path explicitly:

```datagrok-exec
// BitSet path — no widget, just write the mask. Same call as
// grokky.filterSubstructure(t, 'smiles', 'c1ccccc1', {ui: false}).
const bs = await grok.chem.searchSubstructure(t.col('smiles'), 'c1ccccc1');
t.filter.copyFrom(bs);
```

## Combining filter conditions

| Intent                                                   | Code                                                                  |
|----------------------------------------------------------|-----------------------------------------------------------------------|
| Lipinski-style AND in one pass (fastest)                 | `grokky.filterByPredicate(t, (i) => mw.get(i) < 500 && logP.get(i) < 5)` |
| AND another pre-built BitSet onto current filter         | `t.filter.and(otherBitSet)`                                           |
| OR with a pre-built BitSet                               | `t.filter.or(otherBitSet)`                                            |
| Subtract (filter out rows that match `other`)            | `t.filter.andNot(otherBitSet)`                                        |
| AND / OR several predicates                              | `grokky.combineFilters(t, 'and', pred1, pred2, ...)`                  |
| Stack two single-column UI filters                       | `await grokky.filterRows(view, 'MW', {max: 500})` then `... 'cLogP'`  |
| Make a fresh mask without attaching                      | `DG.BitSet.create(t.rowCount, (i) => ...)`                            |

```datagrok-exec
// Canonical Lipinski-style combine via predicate. Single pass.
const mw = t.getCol('MW');
const logP = t.getCol('cLogP');
grokky.filterByPredicate(t, (i) => mw.get(i) < 500 && logP.get(i) < 5);
```

```datagrok-exec
// Or via combineFilters when the predicates are already-named functions.
const mwOk = (i) => t.getCol('MW').get(i) < 500;
const logPOk = (i) => t.getCol('cLogP').get(i) < 5;
grokky.combineFilters(t, 'and', mwOk, logPOk);
```

When mixing UI filters with a programmatic mask, remember: UI filters write
into `df.filter` during `onRowsFiltering`. If you write to `df.filter` outside
that event, your write is overwritten on the next filter cycle. Either:

- Subscribe to `onRowsFiltering` and AND your mask in there (collaborative
  filtering — see the lifecycle section).
- Disable the UI filter group first: `view.getFiltersGroup().setActive(false)`.
- Use `grokky.clearFilter(view)` which does the second for you.

All `BitSet` mutators return `this` so chaining works:

```datagrok-exec
t.filter.setAll(true).and(catA).and(catB);
```

## Clearing & inverting

`grokky.clearFilter(target)` is polymorphic:

- `target: DG.DataFrame` → `df.filter.setAll(true)` only.
- `target: DG.TableView` → also disables the UI filter group via
  `view.getFiltersGroup().setActive(false)`, so widget contributions don't
  immediately re-narrow your now-cleared mask. The widgets stay on screen
  (they're not removed); they're just disabled.

```datagrok-exec
// Show every row again. With a view, also disables UI filter widgets.
grokky.clearFilter(view);
```

```datagrok-exec
// DataFrame-only — leaves UI filter widgets active (they'll re-narrow next pass).
grokky.clearFilter(t);
```

`grokky.invertFilter(df)` calls `df.filter.invert()`. Same caveat as above —
if UI filters are active they'll re-clobber on the next `onRowsFiltering`.

```datagrok-exec
// Show only the rows currently hidden. Combine with clearFilter on the view
// if you want the inversion to stick past the next UI filter cycle.
grokky.clearFilter(view);
grokky.filterByPredicate(t, (i) => t.getCol('activity').get(i) > 5);
grokky.invertFilter(t);  // now shows rows where activity <= 5 (or null)
```

## Show-only vs drop-rows — the destructive split

Two distinct intents, opposite polarity, and opposite destructive-ness. Pick
the right one:

| User says...                                | Helper                                  | Polarity                       | Destructive? |
|---------------------------------------------|-----------------------------------------|--------------------------------|--------------|
| "filter to X" / "show only X"               | `grokky.filterByPredicate(t, pred)`     | `pred(i) === true` → **keep**  | no — hide-only |
| "give me a new table of just the filtered rows" | `grokky.filteredDf(t)`              | uses `df.filter`               | no — clones the DF |
| "remove rows where X" / "drop rows where X" | `grokky.dropRows(t, pred)`              | `pred(i) === true` → **drop**  | **yes** — rows gone |
| "delete the currently-hidden rows"          | `grokky.dropRows(t, (i) => !t.filter.get(i))` | invert + drop                | **yes** |

**`filterByPredicate` and `dropRows` have opposite polarity by design.**
`filterByPredicate` matches the natural reading of "filter to rows where X is
true". `dropRows` matches "drop rows where X is true". If you mix them up
you'll either hide everything or delete the wrong rows.

```datagrok-exec
// Non-destructive: hide rows. df is unchanged in size.
grokky.filterByPredicate(t, (i) => t.getCol('category').get(i) === 'X');
return t.filter.trueCount;
```

```datagrok-exec
// Non-destructive: a new DF containing only the currently-passing rows.
// Original `t` is untouched.
const subset = grokky.filteredDf(t);
return subset;
```

```datagrok-exec
// Destructive: rows where category === 'X' are gone from t.
// Note the opposite polarity from filterByPredicate above.
const removed = grokky.dropRows(t, (i) => t.getCol('category').get(i) === 'X');
return {removed};
```

`grokky.filteredDf(df, opts?)` wraps `df.clone(df.filter, opts?.cols ?? null,
opts?.withSelection ?? false)`. `opts.cols` narrows to a column subset;
`opts.withSelection` preserves the selection mask on the clone.

If the user says "remove" or "delete", **confirm intent** before reaching for
`dropRows` when "filter" might be what they meant. Filtering is reversible;
`dropRows` is not.

## Filter event lifecycle

Three events on `DG.DataFrame`:

| Event                | When                                                                   | Use for                                                   |
|----------------------|------------------------------------------------------------------------|-----------------------------------------------------------|
| `onRowsFiltering`    | Filter system is rebuilding the mask. A viewer's chance to AND in.     | Implementing a custom filter that collaborates with UI.   |
| `onRowsFiltered`     | Filter pass complete; mask is now final.                               | Observers acting on the filtered result (count, summary). |
| `onFilterChanged`    | Same as `df.filter.onChanged` — fires on any BitSet mutation.          | Most general — works for direct writes too.               |

`df.rows.requestFilter()` re-triggers the filter pass. Use it after you
change inputs **outside** `df.filter` itself (e.g. a slider that drives your
custom filter). If you wrote directly to `df.filter`, `onChanged` fires
automatically — no `requestFilter` needed.

The canonical collaborative-filtering pattern (write your contribution into
`df.filter` during `onRowsFiltering`):

```datagrok-exec
// Programmatic filter that survives UI filter cycles. The subscription
// closure dies with the exec block — to make it persistent across exec
// blocks, push the subscription onto view.subs[] (auto-disposed on detach).
const sub = t.onRowsFiltering.subscribe((_) => {
  // The mask has already been AND-folded by prior filters. AND our contribution.
  const col = t.getCol('quality_flag');
  const mask = DG.BitSet.create(t.rowCount, (i) => col.get(i) === 'pass');
  t.filter.and(mask);
});
view.subs.push(sub);
t.rows.requestFilter();  // trigger an initial pass
```

Inside a `datagrok-exec` block, a subscription not pushed onto `view.subs[]`
dies when the block closure is GC'd — no leak, but also no persistence past
the block. In a viewer / widget, push onto `viewer.subs[]` so the
subscription cancels on detach (`DG.Widget` does this automatically).

```datagrok-exec
// One-shot: count rows the filter touches each pass.
const sub = t.onRowsFiltered.subscribe((_) => console.log(`visible: ${t.filter.trueCount}`));
view.subs.push(sub);
```

## Anti-patterns

1. **`df.filter.setAll(false)` to "clear" the filter.** That hides every row.
   Use `df.filter.setAll(true)` or `grokky.clearFilter(df)`. (Polarity:
   `true` = passes.)
2. **Passing SMILES / SMARTS as raw `molBlock` to a `SUBSTRUCTURE` filter.**
   The `Chem:substructureFilter` widget reads `state.molBlock`. RDKit
   silently returns zero atoms for non-molblock input → zero matches with no
   error. Use `grokky.filterSubstructure(...)`, which auto-detects and
   converts via `DG.chem.convert`.
3. **`df.filter.init(...)` to *add* a condition on top of existing.** `init`
   zeros the buffer first. To intersect with current, build
   `DG.BitSet.create(df.rowCount, (i) => ...)` and `df.filter.and(that)`.
4. **`df.rows.removeWhere(...)` or `df.rows.removeWhereIdx(...)` when the
   user said "filter".** Destructive — rows are gone. Confirm intent before
   reaching for `grokky.dropRows`. Use `grokky.filterByPredicate` for
   non-destructive "hide" intent.
5. **Writing to `df.filter` outside `onRowsFiltering` when UI filters
   exist.** Your write is overwritten on the next filter cycle. Either
   subscribe to `onRowsFiltering`, or disable the UI filter group first
   (`fg.setActive(false)` / `grokky.clearFilter(view)`).
6. **Subscribing to `df.onFilterChanged` / `onRowsFiltering` inside a viewer
   without pushing the subscription onto `viewer.subs[]`.** Subscription
   leak. (Not an issue inside `datagrok-exec` blocks — the closure is GC'd
   with the block.)
7. **`view.filters(...)`** — deprecated. Use
   `view.getFiltersGroup({createDefaultFilters: false}).updateOrAdd(...)`
   (which `grokky.filterRows` already does).
8. **Calling `view.getFiltersGroup()` without `{createDefaultFilters:
   false}`** when you only want to add one filter. Default `true` creates a
   default histogram/categorical filter for every column. The Chem package
   consistently passes `false` (`chem/package.ts:665` and others); so do we.
9. **`df.clone(df.filter)` to "filter the DF".** Functional but wasteful —
   `clone` allocates a new DF. If the user wants rows hidden, mutate
   `df.filter`. If they want a real new DF (export, downstream pipeline),
   then yes — that's what `grokky.filteredDf(df)` is for.
10. **Forgetting `await` on `grokky.filterRows({substructure: ...})` or
    `grokky.filterSubstructure(...)`.** Both are async because conversion
    routes through Dart-backed `Chem:convertMolNotation`. The promise carries
    the actual filter write; without `await`, the rest of the exec block
    sees the old filter.
11. **Looping `df.rows` to build a filter.** Use `df.filter.init(i => ...)`
    (or `grokky.filterByPredicate`) — reads columns directly and writes the
    buffer in one pass. `RowList` is explicitly not for perf paths.
12. **Treating `BitSet.invert()` / `.and()` / `.or()` like they return a
    new BitSet.** They mutate `this` and return `this` for chaining. Clone
    first (`bs.clone().and(...)`) if you need a fresh mask.

## Out of scope

- **Selection state.** `df.selection` is the next skill (`datagrok-selection`).
  Same BitSet API, but separate semantics — "what is the user pointing at"
  vs "what is currently visible". Mixing the two confuses users.
- **Generic DataFrame cloning.** For "give me a copy of this DF with these
  columns" (no filter involved), see `datagrok-df-and-columns` (`grokky.cloneDf`).
- **Formula columns.** `addCalculatedColumn` and formula-driven columns live
  in `datagrok-calc-column`. Filters can reference any column including
  calculated ones, but the column-add path is elsewhere.
