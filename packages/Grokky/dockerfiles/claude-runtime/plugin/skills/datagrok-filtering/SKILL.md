---
name: datagrok-filtering
description: Filter rows of a Datagrok DataFrame inside a datagrok-exec block through the Filters panel — by range, equals/contains/in-set, multi-value, boolean, free-text row expressions, or substructure (SMILES / SMARTS / molblock). Also covers clearing, inverting, the show-only-filtered vs destructive-drop split, and the filter event lifecycle (onRowsFiltering / onFilterChanged / onRowsFiltered). Use whenever the user says "filter", "show only", "hide rows where", "narrow to subset", "find rows that", "contains", "substructure search", "categorical filter", "range filter", "invert", "clear the filter", "clear filters", "drop rows", or asks for the filtered subset as a new table. Does NOT cover selection (separate skill) or generic DataFrame cloning (datagrok-df-and-columns).
---

# datagrok-filtering

Row-shape filtering on `DG.DataFrame` through the **Filters panel**
(`view.getFiltersGroup(...)` / `view.filters(...)`). Selection lives in
`datagrok-selection`; generic dataFrame and column inspection in `datagrok-df-and-columns`.

Globals inside every `datagrok-exec` block: `grok`, `ui`, `DG`, `view`,
`t` (the current `DG.DataFrame`, when the view is a `TableView`).

## Always filter through the Filters panel

Every user filter goes through the UI FilterGroup
(`view.getFiltersGroup(...).updateOrAdd(...)`), so the filter shows up in the
Filters panel where the user can see, adjust, and remove it. This covers every
condition — range, categorical / equals / in-set, multi-value, boolean,
substructure, and free-text row expressions (`age > 30`, `sex = "M"`) via the
FREE_TEXT filter. For several conditions at once, stack one widget per column.

**Never write `t.filter` directly** (`t.filter.init`, hand-built BitSet masks) to
satisfy a filter request — it bypasses the panel, leaves no visible chip, and is
overwritten the moment any UI filter re-runs. The panel is always enough.

## Quick reference

| Intent                                                | Code                                                                                |
|-------------------------------------------------------|-------------------------------------------------------------------------------------|
| Several filters at once (Filters panel)               | `view.filters({filters: [{type: DG.FILTER_TYPE.HISTOGRAM, column, min, max}, ...]})` |
| Single widget — add or update (Filters panel)         | `view.getFiltersGroup({createDefaultFilters:true}).updateOrAdd(state)` — state shapes below |
| Clear all filters                                     | `view.getFiltersGroup(...).setActive(false)` (or `t.filter.setAll(true)`)           |
| Invert                                                | `t.filter.invert()`                                                                 |
| New table of the currently-visible rows               | `t.clone(t.filter)`                                                                 |
| Destructively drop rows (gone for good)               | `t.rows.removeWhereIdx((i) => pred(i))`                                             |

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
  {type: DG.FILTER_TYPE.HISTOGRAM, column: 'height', min: 120, max: 150},
  {type: DG.FILTER_TYPE.FREE_TEXT},
  {type: DG.FILTER_TYPE.MULTI_VALUE, column: 'sex', mode: 'OR', include: ['F'], exclude: []},
  {type: DG.FILTER_TYPE.CATEGORICAL, column: 'disease'},
]});
```

For a single filter,
`view.getFiltersGroup({createDefaultFilters: true}).updateOrAdd(state)`
adds or updates one widget.
Provide **every** field a state requires — a partial state renders the widget
but silently does not filter.

State shapes:

| State                                                                                                        | Behavior                             |
|---------------------------------------------------------------------------------------------------------------|--------------------------------------|
| `{type: DG.FILTER_TYPE.HISTOGRAM, column: 'height', min: 120, max: 150}`                                | numeric range; **both `min` and `max` are required** — a min-only (or max-only) state does NOT filter. For an open-ended condition, fill the free bound from the column's extent, e.g. `"> X"` → `min: X, max: t.col(name).stats.max` |
| `{type: DG.FILTER_TYPE.CATEGORICAL, column: 'race', selected: ['Asian', 'Black']}`                      | in-set; single element = equals; omit `selected` → widget with no constraint |
| `{type: DG.FILTER_TYPE.MULTI_VALUE, column: 'sex', mode: 'AND'\|'OR', include: ['F'], exclude: []}`     | cells holding several values (column's separator tag, default newline); `include`/`exclude` BOTH required (empty array ok, null not); a value in neither list is unconstrained; a `selected` key (seen in older samples) is ignored |
| `{type: DG.FILTER_TYPE.BOOL_COLUMNS, 'columnless-filter-identifier': 'control', mode: 'AND', true: [true], false: [false]}` | one widget over N bool columns (order = comma-joined name list); `true[i]` keeps rows where the column is true, `false[i]` where false, both true → no constraint |
| `{type: DG.FILTER_TYPE.FREE_TEXT, gridNames: ['age > 30', '* smith']}`                                      | each entry is a row-matcher expression (`sex = "M"`, `height > 180`) or a `* text` wildcard over all columns; use this for computed / cross-column / OR conditions a single widget can't express; `value` only prefills the box, `column` is ignored |
| `{type: DG.FILTER_TYPE.SUBSTRUCTURE, column, molBlock}`                                                  | substructure (see below)             |

For a single categorical pick like "show only females", use one
CATEGORICAL widget: `updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: ['F']})`.

For multi-column AND conditions ("age > 60 AND height < 170"), stack one
widget per column — they collaborate via `onRowsFiltering`, and the user can
adjust each condition in the panel.

**A UI FilterGroup filter applies asynchronously** — `t.filter` updates a frame
or two after `updateOrAdd`, not synchronously, so reading `t.filter.trueCount`
right after `updateOrAdd` sees the stale count; do not fall back to writing
`t.filter` directly.

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
view.getFiltersGroup({createDefaultFilters: true}).updateOrAdd({
  type: DG.FILTER_TYPE.SUBSTRUCTURE,
  column: 'smiles',
  molBlock,
});
```

## Clearing & inverting

**Polarity is the #1 footgun:** in `t.filter`, a `true` bit = row **visible**,
`false` = **hidden**. So "show every row" is `setAll(true)`; `setAll(false)`
hides everything.

```datagrok-exec
// Show every row again. Disable the UI filter widgets so they don't
// immediately re-narrow the mask (they stay on screen, just inactive).
view.getFiltersGroup({createDefaultFilters: true}).setActive(false);
t.filter.setAll(true);
```

```datagrok-exec
// Invert the current mask. If UI filters are active they'll re-clobber the
// inversion on the next filter pass — disable the group first if you want it to persist.
t.filter.invert();
```

## Show-only vs drop-rows — the destructive split

| User says...                                | Code                                                            | Destructive? |
|---------------------------------------------|-----------------------------------------------------------------|--------------|
| "filter to X" / "show only X" / "hide ..."  | `view.getFiltersGroup({createDefaultFilters:true}).updateOrAdd({type, column, ...})` (the matching widget) | no — hide-only |
| "give me a new table of just the filtered rows" | `t.clone(t.filter)`                                         | no — clones the DF |
| "remove rows where X" / "drop rows where X" | `t.rows.removeWhereIdx((i) => pred(i))`                         | **yes** — rows gone |
| "delete the currently-hidden rows"          | `t.rows.removeWhereIdx((i) => !t.filter.get(i))`                | **yes** |

"Filter" / "show only" / "hide" is a **hide-only** (reversible) operation → the
Filters panel. `removeWhereIdx` deletes rows for good. If the user says "remove"
or "delete", **confirm intent** when "filter" might be what they meant.

`t.clone(filter, columns?, withSelection?)` clones the visible rows. Pass an
array of column names to narrow the copy; pass `true` for `withSelection` to
carry the selection mask onto the clone.

## Filter event lifecycle

| Event                | When                                                                   | Use for                                                   |
|----------------------|------------------------------------------------------------------------|-----------------------------------------------------------|
| `onRowsFiltered`     | Filter pass complete; mask is now final.                               | Observers acting on the filtered result (count, summary). |
| `onFilterChanged`    | Same as `df.filter.onChanged` — fires on any BitSet mutation.          | Most general — works for direct writes too.               |

In a viewer / widget, push subscriptions onto `viewer.subs` — `DG.Widget`
auto-cleans those on detach. Inside a `datagrok-exec` block, push onto
`view.subs` to persist beyond the block.

## Anti-patterns

1. **Writing `t.filter` directly (`t.filter.init`, BitSet masks) for a user
   filter** — bypasses the Filters panel and is overwritten next filter cycle.
   Use `view.getFiltersGroup(...).updateOrAdd(...)`.
2. **`t.filter.setAll(false)` to "clear" the filter** — hides every row. Use
   `setActive(false)` on the group, or `t.filter.setAll(true)`.
3. **Passing SMILES / SMARTS as raw `molBlock` to a `SUBSTRUCTURE` filter** —
   RDKit silently returns zero matches. Convert via
   `DG.chem.convert(query, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock)`.
4. **`t.rows.removeWhereIdx(...)` when the user said "filter"** — destructive.
   Use the Filters panel for non-destructive "hide" intent.
5. **`t.clone(t.filter)` to "filter the DF"** when the user wanted rows hidden —
   allocates a new DF for nothing. Reach for it only when they actually want a
   new DataFrame (export, downstream pipeline).
6. **Treating `BitSet.invert()` like it returns a new BitSet** — it mutates
   `this` and returns `this` for chaining.
