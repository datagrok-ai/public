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

## Choosing the path: Filters panel first

Default to the UI FilterGroup (`view.getFiltersGroup(...).updateOrAdd(...)`)
for any condition expressible as per-column widgets — range, categorical /
equals / in-set, multi-value, boolean, substructure. The filter
appears in the Filters panel where the user can see, adjust, and remove it.

Programmatic `t.filter` is for conditions widgets can't express — OR across
columns, computed predicates, regex on high-cardinality columns. It doesn't
show in the Filters panel.

## Quick reference

| Intent                                                | Code                                                                                |
|-------------------------------------------------------|-------------------------------------------------------------------------------------|
| Several filters at once (Filters panel)               | `view.filters({filters: [{type: DG.FILTER_TYPE.HISTOGRAM, columnName, min, max}, ...]})` |
| Single widget — add or update (Filters panel)         | `view.getFiltersGroup({createDefaultFilters:false}).updateOrAdd(state)` — state shapes below |
| Condition widgets can't express (OR across columns, computed) | `t.filter.init((i) => pred(i))` — not shown in the Filters panel |
| Persistent predicate filter (survives UI filter pass) | Subscribe to `t.onRowsFiltering`, AND a fresh `DG.BitSet` inside; push sub to `view.subs` |
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

## Filtering by column value (UI FilterGroup)

Add several filters at once with `view.filters(...)`:

```datagrok-exec
view.filters({filters: [
  {type: DG.FILTER_TYPE.HISTOGRAM, columnName: 'height', min: 120, max: 150},
  {type: DG.FILTER_TYPE.FREE_TEXT},
  {type: DG.FILTER_TYPE.MULTI_VALUE, columnName: 'sex', mode: 'OR', include: ['F'], exclude: []},
  {type: DG.FILTER_TYPE.CATEGORICAL, columnName: 'disease'},
]});
```

For a single filter,
`view.getFiltersGroup({createDefaultFilters: false}).updateOrAdd(state)`
adds or updates one widget. Always pass `createDefaultFilters: false` unless
you actually want a histogram per column. State shapes:

| State                                                                                                        | Behavior                             |
|---------------------------------------------------------------------------------------------------------------|--------------------------------------|
| `{type: DG.FILTER_TYPE.HISTOGRAM, columnName: 'height', min: 120, max: 150}`                                | numeric range                        |
| `{type: DG.FILTER_TYPE.CATEGORICAL, columnName: 'race', selected: ['Asian', 'Black']}`                      | in-set; single element = equals; omit `selected` → widget with no constraint |
| `{type: DG.FILTER_TYPE.MULTI_VALUE, columnName: 'sex', mode: 'AND'\|'OR', include: ['F'], exclude: []}`     | cells holding several values (column's separator tag, default newline); `include`/`exclude` BOTH required (empty array ok, null not); a value in neither list is unconstrained; a `selected` key (seen in older samples) is ignored |
| `{type: DG.FILTER_TYPE.BOOL_COLUMNS, 'columnless-filter-identifier': 'control', mode: 'AND', true: [true], false: [false]}` | one widget over N bool columns (order = comma-joined name list); `true[i]` keeps rows where the column is true, `false[i]` where false, both true → no constraint |
| `{type: DG.FILTER_TYPE.FREE_TEXT, gridNames: ['age > 30', '* smith']}`                                      | each entry is a row-matcher expression (`sex = "M"`, `height > 180`) or a `* text` wildcard over all columns; `value` only prefills the box, `column` is ignored |
| `{type: DG.FILTER_TYPE.SUBSTRUCTURE, columnName, molBlock}`                                                  | substructure (see below)             |

For multi-column AND conditions ("age > 60 AND height < 170"), stack one
widget per column — they collaborate via `onRowsFiltering`, and the user can
adjust each condition in the panel.

## Filtering by predicate (when widgets can't express it)

`t.filter.init(pred)` writes the mask directly — it doesn't show in the
Filters panel. Prefer it over `t.rows.filter(row => ...)`. For value access
inside the predicate (null checks, raw typed-array views), see
`datagrok-df-and-columns`.

```datagrok-exec
// Cross-column computed condition (BMI > 30) — not expressible as
// per-column widgets.
const weight = t.getCol('weight');
const height = t.getCol('height');
t.filter.init((i) => weight.get(i) / Math.pow(height.get(i) / 100, 2) > 30);
```

A bare `t.filter.init(...)` is overwritten the next time the UI filter group
re-runs (any widget add/remove, slider drag, etc.). To make a predicate
filter collaborate with the UI filters — or live longer than the current
exec block — subscribe to `t.onRowsFiltering` and AND your contribution in
there; push the subscription onto `view.subs` so it's auto-disposed on view
detach:

```datagrok-exec
// Persistent predicate filter — re-applied every UI filter pass.
const sub = t.onRowsFiltering.subscribe((_) => {
  const age = t.getCol('age');
  const mask = DG.BitSet.create(t.rowCount, (i) => age.get(i) < 30);
  t.filter.and(mask);
});
view.subs.push(sub);
t.rows.requestFilter();  // kick the first pass
```

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
  columnName: 'smiles',
  molBlock,
});
```

## Combining masks

| Intent                                                   | Code                                                                  |
|----------------------------------------------------------|-----------------------------------------------------------------------|
| AND another pre-built BitSet onto current filter         | `t.filter.and(otherBitSet)`                                           |
| OR with a pre-built BitSet                               | `t.filter.or(otherBitSet)`                                            |
| Subtract (filter out rows that match `other`)            | `t.filter.andNot(otherBitSet)`                                        |
| Make a fresh mask without attaching                      | `DG.BitSet.create(t.rowCount, (i) => ...)`                            |

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
// Destructive: rows where race === 'Other' are gone from t.
// Opposite polarity from filter.init.
const before = t.rowCount;
t.rows.removeWhereIdx((i) => t.getCol('race').get(i) === 'Other');
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
