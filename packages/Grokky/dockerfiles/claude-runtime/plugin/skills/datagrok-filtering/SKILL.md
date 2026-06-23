---
name: datagrok-filtering
description: Filter rows of a Datagrok DataFrame inside a datagrok-exec block — by predicate, range, equals/contains/in-set, regex, substructure (SMILES / SMARTS / molblock), or by combining masks. Also covers clearing, inverting, the show-only-filtered vs destructive-drop split, and the filter event lifecycle (onRowsFiltering / onFilterChanged / onRowsFiltered). Use whenever the user says "filter", "show only", "hide rows where", "narrow to subset", "find rows that", "contains", "regex", "substructure search", "categorical filter", "range filter", "invert", "clear the filter", "clear filters", "drop rows", or asks for the filtered subset as a new table. Does NOT cover selection (separate skill) or generic DataFrame cloning (datagrok-df-and-columns).
---

# datagrok-filtering

Row-shape filtering on `DG.DataFrame` via `df.filter` (a `DG.BitSet`),
`view.getFiltersGroup(...)`, and `df.rows`. Selection lives in
`datagrok-selection`; generic dataFrame and column inspection in `datagrok-df-and-columns`.

**Before writing any predicate that reads column values, open `datagrok-df-and-columns`** — it contains the correct APIs for null checks (`col.isNone(i)`), column lookup, and value access.

Globals inside every `datagrok-exec` block: `grok`, `ui`, `DG`, `view`,
`t` (the current `DG.DataFrame`, when the view is a `TableView`).

## Quick reference

| Intent                                                | Code                                                                                |
|-------------------------------------------------------|-------------------------------------------------------------------------------------|
| Predicate filter (any condition, any column count)    | `t.filter.init((i) => pred(i))`                                                     |
| Persistent predicate filter (survives UI filter pass) | Subscribe to `t.onRowsFiltering`, AND a fresh `DG.BitSet` inside; push sub to `view.subs` |
| Single-column range (UI widget)                       | `view.getFiltersGroup({createDefaultFilters:false}).updateOrAdd({type:DG.FILTER_TYPE.HISTOGRAM, column, min, max})` |
| Single-column categorical / equals / in-set (UI)      | `... .updateOrAdd({type:DG.FILTER_TYPE.CATEGORICAL, column, selected:[...]})`        |
| Free-text "contains"                                  | `... .updateOrAdd({type:DG.FILTER_TYPE.FREE_TEXT, column, value:'acid'})`            |
| Regex on a string column                              | Predicate filter — `FREE_TEXT` doesn't expose a regex shape                         |
| Substructure                                          | `... .updateOrAdd({type:DG.FILTER_TYPE.SUBSTRUCTURE, column, columnName, molBlock})` |
| Clear all filters                                     | `t.filter.setAll(true)` (and `view.getFiltersGroup(...).setActive(false)` if a UI filter group is present) |
| Invert                                                | `t.filter.invert()`                                                                 |
| New table of the currently-visible rows               | `t.clone(t.filter)`                                                                 |
| Destructively drop rows (gone for good)               | `t.rows.removeWhereIdx((i) => pred(i))`                                             |

## The mental model

`df.filter` is a `DG.BitSet` — one bit per row. **Polarity is the #1 footgun:**

| Bit value | Meaning                              |
|-----------|--------------------------------------|
| `true`    | row **passes** the filter (visible)  |
| `false`   | row is **filtered out** (hidden)     |

So "clear the filter" (show every row) is `df.filter.setAll(true)`. Reaching
for `setAll(false)` hides everything.

All bitwise mutators (`and`, `or`, `xor`, `andNot`, `invert`, `setAll`,
`init`, `copyFrom`) are **in-place** and return `this`. Clone first
(`bs.clone().and(...)`) if you need a separate mask.

`df.filter.init(predicate)` **zeros the buffer first**, then applies
`predicate`. To intersect with the existing filter, build a fresh BitSet via
`DG.BitSet.create(df.rowCount, (i) => ...)` and call `df.filter.and(thatBitSet)`.

## Reading filter state

| Need                                  | Use                                          |
|---------------------------------------|----------------------------------------------|
| rows currently passing                | `t.filter.trueCount`                         |
| rows currently hidden                 | `t.filter.falseCount`                        |
| total rows                            | `t.filter.length` (same as `t.rowCount`)     |
| indices of passing rows               | `t.filter.getSelectedIndexes()` (`Int32Array`)|
| is row `i` passing?                   | `t.filter.get(i)`                            |
| any rows passing / hidden?            | `t.filter.anyTrue` / `t.filter.anyFalse`     |

```datagrok-exec
return {visible: t.filter.trueCount, hidden: t.filter.falseCount, total: t.rowCount};
```

## Filtering by predicate

`t.filter.init(pred)` is the canonical fast path: buffer-direct, single
notification, ~10× faster than `t.rows.filter(row => ...)`.

**Polarity (memorize):** `pred(i)` returns `true` to **keep** row `i`.
Opposite of `removeWhereIdx` ("returns `true` to remove").

```datagrok-exec
// Multi-column predicate in a single pass — strictly preferred over building
// two separate BitSets and AND-ing them.
const mw = t.getCol('MW');
const logP = t.getCol('cLogP');
t.filter.init((i) => mw.get(i) < 500 && logP.get(i) < 5);
```

`t.filter.init` **replaces** the current filter; it does not intersect. To
narrow on top of the existing filter, build a fresh mask and AND:

```datagrok-exec
const col = t.getCol('activity');
const mask = DG.BitSet.create(t.rowCount, (i) => col.get(i) > 7);
t.filter.and(mask);
```

### Predicate filter that survives UI filter cycles

A bare `t.filter.init(...)` is overwritten the next time the UI filter
group re-runs (any widget add/remove, slider drag, etc.). To make a
predicate filter **collaborate** with the UI filters, subscribe to
`t.onRowsFiltering` and AND your contribution in there. Push the
subscription onto `view.subs` so it's auto-disposed on view detach.

```datagrok-exec
// Persistent predicate filter — re-applied every UI filter pass.
const sub = t.onRowsFiltering.subscribe((_) => {
  const mw = t.getCol('MW');
  const mask = DG.BitSet.create(t.rowCount, (i) => mw.get(i) < 500);
  t.filter.and(mask);
});
view.subs.push(sub);
t.rows.requestFilter();  // kick the first pass
```

Use this whenever the predicate needs to live longer than the current exec block.

## Filtering by column value (UI FilterGroup)

`view.getFiltersGroup({createDefaultFilters: false}).updateOrAdd(state)`
attaches a single-column filter widget in the side panel. Always pass
`createDefaultFilters: false` unless you actually want a histogram per column.

| Filter shape                                                                                                | Filter type                          |
|-------------------------------------------------------------------------------------------------------------|--------------------------------------|
| `{type: DG.FILTER_TYPE.HISTOGRAM, column: 'MW', min, max}` (numeric)                                        | range / numeric histogram            |
| `{type: DG.FILTER_TYPE.CATEGORICAL, column: 'category', selected: [...]}` (string / categorical)            | categorical / equals / in-set        |
| `{type: DG.FILTER_TYPE.FREE_TEXT, column: 'name', value: 'acid'}` (string substring)                       | contains / free-text                 |
| `{type: DG.FILTER_TYPE.SUBSTRUCTURE, column, columnName, molBlock}` (molecule)                              | substructure (see below)             |

```datagrok-exec
// Range filter on a numeric column — appears as a histogram in the filter panel.
view.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
  type: DG.FILTER_TYPE.HISTOGRAM, column: 'MW', min: 200, max: 500,
});
```

```datagrok-exec
// In-set: keep rows where category is one of these. Single-element `selected`
// is the equals form.
view.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
  type: DG.FILTER_TYPE.CATEGORICAL, column: 'category', selected: ['A', 'B', 'C'],
});
```

```datagrok-exec
// Free-text contains. Case-insensitive substring on a string column.
view.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
  type: DG.FILTER_TYPE.FREE_TEXT, column: 'name', value: 'acid',
});
```

For **regex** on a string column there is no first-class `FREE_TEXT` shape —
use a predicate filter:

```datagrok-exec
const col = t.getCol('name');
const re = /acid$/i;
t.filter.init((i) => {
  const v = col.get(i);
  return v != null && re.test(String(v));
});
```

For multi-column conditions (Lipinski-style "MW < 500 AND cLogP < 5"), use a
predicate filter rather than stacking single-column UI widgets — predicate
runs in one pass.

If you write to `t.filter` outside `onRowsFiltering` while UI filters exist,
your write is overwritten on the next filter cycle. Either subscribe to
`onRowsFiltering`, or disable the UI filter group first:
`view.getFiltersGroup({createDefaultFilters: false}).setActive(false)`.

## Substructure filter

The `Chem:substructureFilter` widget reads `state.molBlock` — passing raw
SMILES or SMARTS into `molBlock` silently matches **zero** rows. Convert first.

- `DG.chem.isMolBlock(s)` — pure JS, cheap (`s.includes('M  END')`).
- `DG.chem.isSmarts(s)` — calls `Chem:isSmarts` (requires the Chem package).
- Otherwise treat as SMILES.
- Convert with `DG.chem.convert(s, fromNotation, toNotation)`. Notations live
  on `DG.chem.Notation` (`Smiles`, `Smarts`, `MolBlock`).

```datagrok-exec
// SMILES or SMARTS → molblock → UI substructure filter.
const query = 'c1ccccc1';
const molBlock = DG.chem.isMolBlock(query)
  ? query
  : DG.chem.convert(query, DG.chem.isSmarts(query) ? DG.chem.Notation.Smarts : DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
view.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
  type: DG.FILTER_TYPE.SUBSTRUCTURE,
  column: 'smiles',
  columnName: 'smiles',
  molBlock,
  molBlockFailover: query,
});
```

## Combining filter conditions

| Intent                                                   | Code                                                                  |
|----------------------------------------------------------|-----------------------------------------------------------------------|
| Lipinski-style AND in one pass (fastest)                 | `t.filter.init((i) => mw.get(i) < 500 && logP.get(i) < 5)`            |
| AND another pre-built BitSet onto current filter         | `t.filter.and(otherBitSet)`                                           |
| OR with a pre-built BitSet                               | `t.filter.or(otherBitSet)`                                            |
| Subtract (filter out rows that match `other`)            | `t.filter.andNot(otherBitSet)`                                        |
| Stack two single-column UI filters                       | two `updateOrAdd` calls — they collaborate via `onRowsFiltering`      |
| Make a fresh mask without attaching                      | `DG.BitSet.create(t.rowCount, (i) => ...)`                            |

`BitSet` mutators return `this`, so chaining works:

```datagrok-exec
const catA = DG.BitSet.create(t.rowCount, (i) => i % 2 === 0);
const catB = DG.BitSet.create(t.rowCount, (i) => i % 3 === 0);
t.filter.setAll(true).and(catA).and(catB);
```

## Clearing & inverting

```datagrok-exec
// Show every row again. setAll(true) is the polarity-correct form;
// setAll(false) would hide everything.
t.filter.setAll(true);
// If UI filter widgets exist, also disable them so they don't
// immediately re-narrow the mask. The widgets stay on screen — just inactive.
view.getFiltersGroup({createDefaultFilters: false}).setActive(false);
```

```datagrok-exec
// Invert the current mask. If UI filters are active they'll re-clobber the
// inversion on the next onRowsFiltering pass — disable the group first if
// you want the inversion to persist.
t.filter.invert();
```

## Show-only vs drop-rows — the destructive split

| User says...                                | Code                                                            | Polarity                       | Destructive? |
|---------------------------------------------|-----------------------------------------------------------------|--------------------------------|--------------|
| "filter to X" / "show only X" / "hide ..."  | `t.filter.init((i) => pred(i))`                                 | `pred(i) === true` → **keep**  | no — hide-only |
| "give me a new table of just the filtered rows" | `t.clone(t.filter)`                                         | uses `t.filter`                | no — clones the DF |
| "remove rows where X" / "drop rows where X" | `t.rows.removeWhereIdx((i) => pred(i))`                         | `pred(i) === true` → **drop**  | **yes** — rows gone |
| "delete the currently-hidden rows"          | `t.rows.removeWhereIdx((i) => !t.filter.get(i))`                | inverted, drop                 | **yes** |

**`filter.init` and `removeWhereIdx` have opposite polarity.** `filter.init`
matches "filter to rows where X is true"; `removeWhereIdx` matches "drop rows
where X is true". Confusing them either hides every row or deletes the wrong ones.

```datagrok-exec
// Destructive: rows where category === 'X' are gone from t.
// Opposite polarity from filter.init.
const before = t.rowCount;
t.rows.removeWhereIdx((i) => t.getCol('category').get(i) === 'X');
return {removed: before - t.rowCount};
```

If the user says "remove" or "delete", **confirm intent** when "filter" might
be what they meant. Filtering is reversible; `removeWhereIdx` is not.

`t.clone(filter, columns?, withSelection?)` clones the visible rows. Pass an
array of column names to narrow the copy; pass `true` for `withSelection` to
carry the selection mask onto the clone.

## Filter event lifecycle

| Event                | When                                                                   | Use for                                                   |
|----------------------|------------------------------------------------------------------------|-----------------------------------------------------------|
| `onRowsFiltering`    | Filter system is rebuilding the mask. A viewer's chance to AND in.     | Implementing a custom filter that collaborates with UI.   |
| `onRowsFiltered`     | Filter pass complete; mask is now final.                               | Observers acting on the filtered result (count, summary). |
| `onFilterChanged`    | Same as `df.filter.onChanged` — fires on any BitSet mutation.          | Most general — works for direct writes too.               |

`t.rows.requestFilter()` re-triggers the filter pass. Use it after you change
inputs **outside** `t.filter` itself (e.g. a slider that drives your custom
filter). If you wrote directly to `t.filter`, `onChanged` fires automatically.

In a viewer / widget, push subscriptions onto `viewer.subs` — `DG.Widget`
auto-cleans those on detach. Inside a `datagrok-exec` block, push onto
`view.subs` to persist beyond the block.

## Anti-patterns

1. **`t.filter.setAll(false)` to "clear" the filter** — hides every row.
   Use `t.filter.setAll(true)`.
2. **Passing SMILES / SMARTS as raw `molBlock` to a `SUBSTRUCTURE` filter** —
   RDKit silently returns zero matches. Convert via
   `DG.chem.convert(query, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock)`,
   or use the BitSet path `grok.chem.searchSubstructure`.
3. **`t.filter.init(...)` to *add* a condition on top of existing** — `init`
   zeros the buffer first. Build `DG.BitSet.create(t.rowCount, (i) => ...)`
   and `t.filter.and(thatMask)` instead.
4. **`t.rows.removeWhereIdx(...)` when the user said "filter"** —
   destructive. Use `t.filter.init(...)` for non-destructive "hide" intent.
5. **Writing to `t.filter` outside `onRowsFiltering` when UI filters exist** —
   your write is overwritten next filter cycle. Either subscribe to
   `onRowsFiltering`, or `view.getFiltersGroup({createDefaultFilters: false}).setActive(false)` first.
6. **Calling `view.getFiltersGroup()` without `{createDefaultFilters: false}`**
   when you only want to add one filter — defaults to creating a filter for
   every column.
7. **`t.clone(t.filter)` to "filter the DF"** when the user wanted rows
   hidden — allocates a new DF for nothing. Reach for `t.clone(t.filter)` only
   when they actually want a new DataFrame (export, downstream pipeline).
8. **Treating `BitSet.invert()` / `.and()` / `.or()` like they return a new
   BitSet.** They mutate `this` and return `this` for chaining.
