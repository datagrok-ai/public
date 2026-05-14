---
name: datagrok-selection
description: Manipulate row selection (`df.selection`) on a Datagrok DataFrame inside a datagrok-exec block — set, clear, invert, add to, remove from, intersect, and read the selection mask. Also covers current row (`df.currentRowIdx`), the cross-skill bridges to/from the filter, and how to materialize the selected rows as a new DataFrame. Use whenever the user says "select", "deselect", "highlight rows", "selected rows", "clear selection", "invert selection", "what's selected", "current row", or asks for selected rows as a new table. Does NOT cover row filtering (separate skill `datagrok-filtering`) or generic DataFrame cloning (`datagrok-df-and-columns`).
---

# datagrok-selection

Use the `grokky.*` selection helpers inside a `datagrok-exec` block. They
wrap `df.selection` (a `DG.BitSet`), the current-row pointer, and the bridges
to `df.filter` so Claude doesn't have to remember that **`setAll(false)` is
how you clear the selection** (the opposite polarity from `clearFilter`!),
that `df.selection.init` zeros the buffer before applying a predicate, that
selection and "current row" are distinct concepts, or that
`df.selection.copyFrom(df.filter)` is the canonical "select what's currently
visible" idiom.

## What this skill covers

Row-shape selection state on `DG.DataFrame`: reading the selection BitSet,
writing a programmatic selection by BitSet / predicate / index list /
`Int32Array`, clearing / inverting / select-all, the four selection modes
(`replace` / `add` / `remove` / `intersect`), materializing the selected rows
as a new DataFrame via `selectedDf`, the cross-skill bridges
(`filterFromSelection`, `selectionFromFilter`), the distinct concept of
current row (`df.currentRowIdx`), and subscribing to `onSelectionChanged`.

**Out of scope.** Filter state (`df.filter`) lives in the `datagrok-filtering`
skill — same BitSet API, **opposite polarity for clearing**, see the table
below. Generic DataFrame cloning (no selection involved) lives in
`datagrok-df-and-columns` (`grokky.cloneDf`). Viewer-driven interactive
selection (lasso, click, shift-click) is grid behaviour, not script-time —
covered by the viewers skill, with a footnote here for discoverability.

## Quick reference

| Helper                                            | One-liner                                                                            |
|---------------------------------------------------|--------------------------------------------------------------------------------------|
| `grokky.selectRows(df, input, opts?)`             | Polymorphic: BitSet / predicate / index list / `Int32Array` × 4 modes.               |
| `grokky.clearSelection(df)`                       | Deselect every row. `df.selection.setAll(false)` — **opposite polarity from filter**.|
| `grokky.selectAll(df)`                            | Mark every row selected. Selects filtered-out rows too.                              |
| `grokky.invertSelection(df)`                      | Flip the current mask.                                                               |
| `grokky.selectedDf(df, opts?)`                    | Clone of currently-selected rows. Symmetric to `filteredDf`.                         |
| `grokky.filterFromSelection(df)`                  | Show only the selected rows. `df.filter.copyFrom(df.selection)`.                     |
| `grokky.selectionFromFilter(df)`                  | Select every currently-visible row. `df.selection.copyFrom(df.filter)`.              |
| `grokky.setCurrentRow(df, idx)`                   | Move the current-row pointer. Not selection — distinct concept.                      |
| `grokky.describeSelection(df)`                    | JSON summary: `{count, total, indexes, currentRowIdx, sample?}`.                     |

Globals available inside every `datagrok-exec` block: `grok`, `ui`, `DG`,
`view`, `t` (the current `DG.DataFrame`, when the view is a TableView),
`grokky`.

## The polarity table — opposite of filter

**This is the #1 trap in this skill.** Selection and filter use the same
`DG.BitSet` class, but the *meaning of `true`* is different, so "clear" goes
in opposite directions:

| Action                               | Filter (`df.filter`)     | Selection (`df.selection`)  |
|--------------------------------------|--------------------------|-----------------------------|
| clear (return to default)            | `setAll(true)` (show all)| `setAll(false)` (none selected) |
| "select all" / "filter to all"       | `setAll(true)`           | `setAll(true)`              |
| meaning of `bit === true`            | row **passes** (visible) | row is **selected**         |
| meaning of `bit === false`           | row hidden               | row not selected            |

The pair to memorize, side by side:

```ts
df.filter.setAll(true);      // clear filter  → every row visible
df.selection.setAll(false);  // clear selection → none selected
```

If Claude reaches for `df.selection.setAll(true)` to "clear" the selection,
that **selects every row** instead. Use `grokky.clearSelection(df)`, or write
`df.selection.setAll(false)` explicitly — the helper exists to kill this trap.

The predicate polarity, in contrast, matches across skills: in both
`filterByPredicate` and `selectRows(..., pred)`, the predicate returns `true`
to *include* the row. Only the *clear* idiom flips.

## Selection vs current row vs mouse-over

Datagrok tracks three independent row pointers. They are distinct concepts —
none of them is "the selected row" in the singular:

| Concept       | API                                        | Cardinality            | Event                          |
|---------------|--------------------------------------------|------------------------|--------------------------------|
| Selection     | `df.selection` (BitSet)                    | many rows (set)        | `df.onSelectionChanged`        |
| Current row   | `df.currentRowIdx` / `df.currentRow`       | one row (focus)        | `df.onCurrentRowChanged`       |
| Mouse-over    | `df.mouseOverRowIdx`                       | one row (transient)    | `df.onMouseOverRowChanged`     |

Setting `df.currentRowIdx = 5` does **not** select row 5. Selecting rows does
**not** move the current row. `SYNC_TYPE.CURRENT_ROW_TO_SELECTION` exists
exactly because the link is opt-in.

When a user says "highlight row 5", the intent is ambiguous between *focus*
(current row) and *selection set membership*. Pick the right one or ask. The
SKILL examples below pick selection for "select row 5" and current-row for
"navigate to row 5".

## Reading selection state

| Need                           | Use                                          |
|--------------------------------|----------------------------------------------|
| rows currently selected        | `df.selection.trueCount`                     |
| rows not selected              | `df.selection.falseCount`                    |
| total rows                     | `df.selection.length` (same as `df.rowCount`)|
| indices of selected rows       | `df.selection.getSelectedIndexes()` (Int32Array) |
| is row `i` selected?           | `df.selection.get(i)`                        |
| any row selected?              | `df.selection.anyTrue`                       |
| any row unselected?            | `df.selection.anyFalse`                      |
| find next / prev set bit       | `df.selection.findNext(i, true)` / `findPrev(i, true)` |

```datagrok-exec
// How many rows are currently selected?
return {selected: t.selection.trueCount, unselected: t.selection.falseCount, total: t.rowCount};
```

```datagrok-exec
// Indices of the currently-selected rows. Fast — getSelectedIndexes is cached.
return Array.from(t.selection.getSelectedIndexes());
```

```datagrok-exec
// Is row 17 selected right now?
return t.selection.get(17);
```

`getSelectedIndexes()` caches its result. The cache invalidates only on the
public mutators (`set`, `init`, `setAll`, `and/or/xor/andNot`, `invert`,
`copyFrom`). If you mutate via a direct buffer write (`getBuffer()`), call
`fireChanged()` to invalidate.

## Programmatic selection

`grokky.selectRows(df, input, opts?)` is the single polymorphic entry point.
It absorbs four input shapes and four modes.

**Input shapes:**

| Shape                                    | Example                                            |
|------------------------------------------|----------------------------------------------------|
| `DG.BitSet`                              | a mask you built elsewhere                         |
| `(i: number) => boolean`                 | predicate over row index                           |
| `number[]`                               | explicit row indices                               |
| `Int32Array` / `Uint32Array` / `ArrayLike<number>` | `getSelectedIndexes()` output, typed-array indices |

**Modes** (`opts.mode`, default `'replace'`):

| Mode         | Meaning                                  | BitSet op                       |
|--------------|------------------------------------------|---------------------------------|
| `'replace'`  | New selection replaces current (default) | `init` / `copyFrom` / setAll+set |
| `'add'`      | Union with current selection             | `or`                            |
| `'remove'`   | Subtract from current selection          | `andNot`                        |
| `'intersect'`| Intersect with current selection         | `and`                           |

```datagrok-exec
// Select rows where age > 40. Default mode = 'replace'.
grokky.selectRows(t, (i) => (t.getCol('age').get(i)) > 40);
return {selected: t.selection.trueCount};
```

```datagrok-exec
// Select explicit indices. Empty array with default 'replace' clears.
grokky.selectRows(t, [0, 3, 5, 7, 11]);
```

```datagrok-exec
// Add to the current selection (don't replace it).
grokky.selectRows(t, (i) => (t.getCol('activity').get(i)) > 7, {mode: 'add'});
```

```datagrok-exec
// Remove a row range from the selection.
grokky.selectRows(t, [5, 7], {mode: 'remove'});
```

```datagrok-exec
// Intersect with another mask — keep only currently-selected rows that also
// match this predicate.
grokky.selectRows(t, (i) => t.getCol('category').get(i) === 'A', {mode: 'intersect'});
```

A predicate with `mode: 'replace'` routes through `df.selection.init(pred)` —
single buffer-direct pass, single notification. With the other modes the
helper materializes a fresh `DG.BitSet.create(df.rowCount, pred)` then ORs /
AND-NOTs / ANDs it onto `df.selection`.

`mode: 'replace'` with an empty input shape (`[]` / `Int32Array(0)`) clears the
selection. `mode: 'add'` / `'remove'` / `'intersect'` with empty input is a
no-op. Both are intentional — see anti-patterns below.

## Clearing & inverting

```datagrok-exec
// Deselect everything. SETALL(FALSE), NOT TRUE.
grokky.clearSelection(t);
```

```datagrok-exec
// Select every row in the DF (including filtered-out rows).
grokky.selectAll(t);
```

```datagrok-exec
// Flip the current selection.
grokky.invertSelection(t);
```

**Polarity warning, repeated for emphasis:**

```ts
grokky.clearFilter(view);    // shows every row    → setAll(TRUE) under the hood
grokky.clearSelection(t);    // deselects every row → setAll(FALSE) under the hood
```

The two helpers exist precisely to prevent Claude from substituting one
polarity for the other.

`grokky.selectAll(df)` is `df.selection.setAll(true)`. **It selects every row
in the DF — including filtered-out rows.** For "select all currently-visible
rows", use `grokky.selectionFromFilter(df)` (see Cross-skill bridges below).

## Selected rows as a new dataframe

`grokky.selectedDf(df, opts?)` returns a `DG.DataFrame` containing just the
currently-selected rows. Mirrors `filteredDf` from the filtering skill — same
options shape, different mask:

```ts
selectedDf(df, opts?: {cols?: string[]; withSelection?: boolean}): DG.DataFrame;
```

```datagrok-exec
// New DF with just the selected rows. Source `t` is untouched.
const subset = grokky.selectedDf(t);
return subset;
```

```datagrok-exec
// Narrow to two columns.
const subset = grokky.selectedDf(t, {cols: ['smiles', 'activity']});
return subset;
```

Contrast with `grokky.filteredDf(df)` from the filtering skill:

| Helper              | Mask used         | Source DF modified? |
|---------------------|-------------------|---------------------|
| `grokky.filteredDf` | `df.filter`       | no — clone of visible rows |
| `grokky.selectedDf` | `df.selection`    | no — clone of selected rows |

Same shape, different mask. If the selection is empty, `selectedDf` returns a
zero-row clone (no error, no toast). For user-facing flows, guard with
`df.selection.anyTrue` first.

## Cross-skill bridges

Two one-line `copyFrom` helpers that answer common verbal requests:

| Helper                          | Effect                                              | Phrasing it answers                                |
|---------------------------------|-----------------------------------------------------|----------------------------------------------------|
| `grokky.filterFromSelection(df)`| `df.filter.copyFrom(df.selection)`                  | "hide the un-selected rows", "filter (in place) to the selection" |
| `grokky.selectionFromFilter(df)`| `df.selection.copyFrom(df.filter)`                  | "select all currently-visible rows"                |

**Critical disambiguation — don't confuse with `selectedDf`:**

- "show only the selected rows **as a new table**" → `selectedDf(df)` (returns a *new* DataFrame; original unchanged).
- "copy / subset / extract the selected rows" → `selectedDf(df)`.
- "show only the selected rows" (no "new table") → ambiguous; default to `selectedDf` and surface the in-place option, OR ask. `filterFromSelection` MUTATES the current view's filter — it does NOT create a new table.

If the request mentions "new table", "subset", "copy", "extract", or asks to inspect/operate on the selected rows separately → ALWAYS `selectedDf`, never `filterFromSelection`.

These work because both masks use the same polarity in their respective
domains — `true` means "yes" in both ("visible" and "selected"). No
inversion needed.

```datagrok-exec
// "Filter down to the selected rows" — pipe selection → filter.
grokky.filterFromSelection(t);
return {visible: t.filter.trueCount, selected: t.selection.trueCount};
```

```datagrok-exec
// "Select everything I can currently see" — pipe filter → selection.
grokky.selectionFromFilter(t);
return {selected: t.selection.trueCount};
```

`selectionFromFilter` fires `onSelectionChanged` but does **not** fire a
filter event — it reads `df.filter` and writes `df.selection`. Likewise
`filterFromSelection` fires `onFilterChanged` (and the filter lifecycle) but
not `onSelectionChanged`.

## Current row

`grokky.setCurrentRow(df, idx)` moves the current-row pointer:

```datagrok-exec
// Focus row 5 (the platform usually does this on click; rarely needed in scripts).
grokky.setCurrentRow(t, 5);
```

**The platform already sets `currentRowIdx` on most user interactions** —
clicking a row in the grid, navigating with arrow keys, clicking a viewer
data point. Programmatic `setCurrentRow` is for scripted navigation
("scroll to the first row matching X").

Out-of-range index throws (`idx < 0` or `idx >= df.rowCount`). Setting to
`-1` is a valid "no current row" state, but use `df.currentRowIdx = -1`
directly for that — the helper guards against accidental negatives.

`df.currentRow` (a `Row` proxy) and `df.currentCell` (a `Cell`) are read on
demand. They follow `currentRowIdx`; setting them via `df.currentRow = row`
works because Dart normalizes through `row.idx`, but `setCurrentRow(df, idx)`
is the cleaner API from a script.

## Selection event lifecycle

One event — unlike filter, which has three (`onRowsFiltering`,
`onFilterChanged`, `onRowsFiltered`). Selection is just state; there are no
collaborators contributing into it, so there's no two-phase lifecycle:

| Event                       | When                                                        |
|-----------------------------|-------------------------------------------------------------|
| `df.onSelectionChanged`     | Any mutation on `df.selection` (`set`, `init`, `setAll`, `and/or/xor/andNot`, `invert`, `copyFrom`). Same observable as `df.selection.onChanged`. |

```datagrok-exec
// Log selection count whenever it changes. Subscription lives until view detach.
const sub = t.onSelectionChanged.subscribe((_) => {
  console.log(`selected: ${t.selection.trueCount}`);
});
view.subs.push(sub);
```

Inside a `datagrok-exec` block, a subscription not pushed onto `view.subs[]`
dies when the block closure is GC'd — no leak, but no persistence past the
block. In a viewer / widget, push onto `viewer.subs[]` so the subscription
cancels on detach (`DG.Widget` does this automatically).

**Coalescing.** Every `BitSet` mutator takes a `notify: boolean = true`
parameter. The default fires immediately on each call. The deferral pattern:

```ts
for (let i = 0; i < indices.length; i++)
  df.selection.set(indices[i], true, false);   // notify = false
df.selection.fireChanged();                    // one coalesced event
```

`grokky.selectRows` already does this for the index-array path. You only need
the manual pattern if you're writing per-bit logic outside the helper.

## Common patterns

```datagrok-exec
// Describe the current selection — useful for "what is selected?" prompts.
return grokky.describeSelection(t);
```

```datagrok-exec
// Make a new table view from the selected rows.
const sub = grokky.selectedDf(t);
grok.shell.addTableView(sub);
```

```datagrok-exec
// Select via column predicate — most common shape.
grokky.selectRows(t, (i) => (t.getCol('activity').get(i)) > 7);
```

```datagrok-exec
// Toggle selection on a single row.
const i = 5;
t.selection.set(i, !t.selection.get(i));
```

```datagrok-exec
// "Select what I see" → "filter to what I just selected" pipeline.
grokky.selectionFromFilter(t);                                   // grab visible → selection
grokky.selectRows(t, (i) => t.getCol('flag').get(i) === 'ok', {mode: 'intersect'});  // narrow
grokky.filterFromSelection(t);                                   // push back to filter
```

## Anti-patterns

1. **`df.selection.setAll(true)` to "clear" the selection.** That selects
   every row — the **opposite** of what filter's clear does. Use
   `df.selection.setAll(false)` or `grokky.clearSelection(df)`. (This is the
   #1 footgun in this skill; the helper exists to absorb it.)
2. **`df.filter.setAll(false)` to "clear" the filter.** That hides every
   row — wrong direction for filter too. Use `df.filter.setAll(true)` or
   `grokky.clearFilter(df)`. The cross-skill pair must be memorized together.
3. **"Highlight row 5" treated as selection without considering current row.**
   The user might mean focus (`df.currentRowIdx = 5`), not selection. Ask, or
   pick the one that matches the surrounding context (e.g. if they said
   "scroll to row 5" → current row; if they said "mark row 5" → selection).
4. **`df.rows.select(rowPred)` for "replace selection" in perf-critical
   code.** Functional, but iterates the `RowList` and per-bit `set` —
   materially slower than `grokky.selectRows(df, pred, {mode: 'replace'})`,
   which routes through `df.selection.init` (buffer-direct, single pass).
5. **`df.selection.init(...)` to *add* a condition to the selection.** `init`
   zeros the buffer first. To extend, build `DG.BitSet.create(df.rowCount,
   pred)` and OR onto `df.selection` — or use `selectRows(df, pred, {mode:
   'add'})`.
6. **Per-bit `set` in a loop without `notify=false`.** N writes → N events →
   laggy UI. Always `set(i, x, false)` in the loop and `fireChanged()` at
   the end. `selectRows` already does this for the array shape.
7. **Holding a stale `getSelectedIndexes()` snapshot across mutations.** The
   array is cached; the cache invalidates on the public mutators, but if you
   captured it once and the selection has since changed, your local copy is
   stale. Re-read between mutations, or rely on `df.selection.get(i)` for
   point queries. (Holding the `df.selection` BitSet wrapper itself across
   mutations is fine — it points at the same underlying Dart object.)
8. **`grokky.selectAll(df)` when the user means "select all visible".**
   `selectAll` selects every row in the DF, **including filtered-out rows**.
   For "select all currently-visible rows", use `grokky.selectionFromFilter`.
9. **Combining selections across DataFrames of different row counts.**
   `dfA.selection.and(dfB.selection)` only makes sense if rows are aligned
   1:1. The Dart side may not validate — silent garbage. Cross-table sync by
   key columns is `grok.data.linkTables(...)`, which is out of scope for
   one-shot `datagrok-exec` blocks.
10. **Mutating `df.selection` and expecting the filter to update.** They're
    independent. Use `grokky.filterFromSelection(df)` (or write
    `df.filter.copyFrom(df.selection)`) to bridge.
11. **Subscribing to `df.onSelectionChanged` inside a viewer without pushing
    onto `viewer.subs[]`.** Subscription leak. (Not an issue inside
    `datagrok-exec` blocks — the closure dies with the block, which kills
    the sub.)

## Out of scope

- **Filter operations.** `df.filter` lives in `datagrok-filtering`. Same
  BitSet API, **opposite clear polarity** — see the polarity table above.
- **Viewer-driven interactive selection.** Lasso, grid click, shift-extend,
  ctrl-toggle — that's grid/viewer behaviour, dispatched via
  `BitSet.handleClick(rowPred, mouseEvent, ...)`. Out of scope for
  `datagrok-exec` (no mouse event in hand). Covered by the viewers skill.
- **Cross-table selection sync.** `grok.data.linkTables(t1, t2, k1, k2,
  [SYNC_TYPE.SELECTION_TO_SELECTION, ...])` establishes long-lived sync
  between two tables. Out of scope for one-shot blocks; document and use the
  `grok.data` API directly when needed.
- **Generic DataFrame cloning.** For "give me a copy of this DF with these
  columns" (no selection involved), see `datagrok-df-and-columns`
  (`grokky.cloneDf`).
- **Selection groups.** The platform has no "selection group" concept like
  it has filter groups; selection is one BitSet per DataFrame. If a user
  asks for "named selections", that's a layout / bookmark feature, not a
  scripting one.
