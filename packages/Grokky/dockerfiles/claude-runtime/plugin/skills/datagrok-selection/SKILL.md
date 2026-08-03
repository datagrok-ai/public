---
name: datagrok-selection
description: Manipulate row selection (`df.selection`) on a Datagrok DataFrame inside a datagrok-exec block — set, clear, invert, add to, remove from, intersect, and read the selection mask. Also covers current row (`df.currentRowIdx`), the cross-skill bridges to/from the filter, and how to materialize the selected rows as a new DataFrame. Use whenever the user says "select rows", "deselect", "invert selection", "highlight rows", "selected rows", "clear selection", "current row", "currentRowIdx", "count selected", "list selected indexes", "selection from filter", "filter from selection", or asks for selected rows as a new table. Does NOT cover row filtering (separate skill `datagrok-filtering`) or generic DataFrame cloning (`datagrok-df-and-columns`).
---

# datagrok-selection

The active DataFrame is `t`; row selection lives at `t.selection` (a `DG.BitSet`).

The number-one trap: **`t.selection.setAll(false)` clears the selection**.
`setAll(true)` *selects* every row. That is the opposite of filter, where
`setAll(true)` is the cleared / no-filter state.

## Quick reference

| Need                                 | Code                                                                                   |
|--------------------------------------|----------------------------------------------------------------------------------------|
| select all rows                      | `t.selection.setAll(true)`                                                             |
| **clear** the selection              | `t.selection.setAll(false)`                                                            |
| invert                               | `t.selection.invert()`                                                                 |
| select by predicate (replace)        | `t.selection.init(i => pred(i))`                                                       |
| select by index list (replace)       | `t.selection.setAll(false, false); for (const i of idx) t.selection.set(i,true,false); t.selection.fireChanged();` |
| add by predicate (union)             | `t.selection.or(DG.BitSet.create(t.rowCount, i => pred(i)))`                           |
| remove indices                       | `for (const i of idx) t.selection.set(i,false,false); t.selection.fireChanged();`      |
| intersect with predicate             | `t.selection.and(DG.BitSet.create(t.rowCount, i => pred(i)))`                          |
| count selected                       | `t.selection.trueCount`                                                                |
| indices of selected rows             | `t.selection.getSelectedIndexes()` (returns `Int32Array`)                              |
| is row `i` selected?                 | `t.selection.get(i)`                                                                   |
| any row selected?                    | `t.selection.anyTrue`                                                                  |
| selected rows as a new DF            | `t.clone(t.selection)`                                                                 |
| filter → selection (select visible)  | `t.selection.copyFrom(t.filter)`                                                       |
| selection → filter (hide unselected) | `t.filter.copyFrom(t.selection)`                                                       |
| move current row                     | `t.currentRowIdx = i` (validate `0 <= i < t.rowCount`)                                 |
| current row object / cell            | `t.currentRow`, `t.currentCell`                                                        |
| mouse-over row                       | `t.mouseOverRowIdx`                                                                    |
| subscribe to selection changes       | `view.subs.push(t.onSelectionChanged.subscribe(() => …))`                              |

## The polarity table — opposite of filter

Selection and filter use the same `DG.BitSet` class, but the *meaning of `true`*
differs, so "clear" goes in opposite directions:

| Action                       | Filter (`df.filter`)         | Selection (`df.selection`)      |
|------------------------------|------------------------------|---------------------------------|
| clear (return to default)    | `setAll(true)` (show all)    | `setAll(false)` (none selected) |
| meaning of bit === `true`    | row **passes** (visible)     | row is **selected**             |
| meaning of bit === `false`   | row hidden                   | row not selected                |

The pair to memorize, side by side:

```ts
t.filter.setAll(true);      // clear filter    → every row visible
t.selection.setAll(false);  // clear selection → none selected
```

Predicate polarity matches across skills: in both `t.filter.init(pred)` and
`t.selection.init(pred)`, returning `true` *includes* the row. Only the *clear*
idiom flips.

## Selection vs current row vs mouse-over

Three independent row pointers — none is "the selected row" in the singular:

| Concept       | API                                  | Cardinality       | Event                       |
|---------------|--------------------------------------|-------------------|-----------------------------|
| Selection     | `t.selection` (BitSet)               | many rows (set)   | `t.onSelectionChanged`      |
| Current row   | `t.currentRowIdx` / `t.currentRow`   | one row (focus)   | `t.onCurrentRowChanged`     |
| Mouse-over    | `t.mouseOverRowIdx`                  | one row (transient) | `t.onMouseOverRowChanged` |

Setting `t.currentRowIdx = 5` does **not** select row 5; selecting rows does
**not** move the current row. When a user says "highlight row 5", the intent is
ambiguous — pick the path matching context ("scroll to" → current row;
"mark" → selection) or surface the ambiguity.

## Reading selection state

| Need                      | Code                                              |
|---------------------------|---------------------------------------------------|
| rows currently selected   | `t.selection.trueCount`                           |
| rows not selected         | `t.selection.falseCount`                          |
| total rows                | `t.selection.length` (same as `t.rowCount`)       |
| indices of selected rows  | `t.selection.getSelectedIndexes()` → `Int32Array` |
| is row `i` selected?      | `t.selection.get(i)`                              |
| any selected?             | `t.selection.anyTrue`                             |
| any unselected?           | `t.selection.anyFalse`                            |
| next / prev set bit       | `t.selection.findNext(i, true)` / `findPrev(i, true)` |

```datagrok-exec
// How many rows are currently selected?
return {selected: t.selection.trueCount, unselected: t.selection.falseCount, total: t.rowCount};
```

`getSelectedIndexes()` caches its result. The cache invalidates on public
mutators (`set`, `init`, `setAll`, `and/or/xor/andNot`, `invert`, `copyFrom`).
If you mutate via `getBuffer()`, call `fireChanged()` to invalidate.

## Programmatic selection — modes

The four modes map directly to BitSet operations:

| Mode        | BitSet op                              | Pattern                              |
|-------------|----------------------------------------|--------------------------------------|
| replace     | `init(pred)` / `copyFrom(other)`       | new selection replaces current       |
| add         | `or(other)`                            | union with current                   |
| remove      | `andNot(other)`                        | subtract from current                |
| intersect   | `and(other)`                           | keep only rows in both               |

### Predicate

`t.selection.init(pred)` is the fast, buffer-direct path for **replace**:

```datagrok-exec
// Select rows where age > 40 (replace).
const age = t.getCol('age');
t.selection.init(i => (age.get(i)) > 40);
return {selected: t.selection.trueCount};
```

For add / remove / intersect with a predicate, build a fresh BitSet and combine:

```datagrok-exec
// Add to the current selection (union, don't replace).
const act = t.getCol('activity');
t.selection.or(DG.BitSet.create(t.rowCount, i => (act.get(i)) > 7));
```

### Index list

`t.selection.set(i, x, notify)` writes one bit. **Pass `notify=false` in a
loop, then call `fireChanged()` once** — otherwise N writes fire N events:

```datagrok-exec
// Select explicit indices (replace). Zero the buffer, set bits with notify=false,
// fire one change at the end.
const idx = [0, 3, 5, 7, 11];
t.selection.setAll(false, false);
for (const i of idx) t.selection.set(i, true, false);
t.selection.fireChanged();
```

**Anti-pattern (laggy UI):**

```ts
// BAD — fires N events.
for (const i of indices) t.selection.set(i, true);
```

### From another BitSet

```ts
t.selection.copyFrom(other);     // replace
t.selection.or(other);           // add
t.selection.andNot(other);       // remove
t.selection.and(other);          // intersect
```

## Clearing & inverting

```datagrok-exec
// Deselect everything. setAll(FALSE), not TRUE.
t.selection.setAll(false);
```

```datagrok-exec
// Flip the current selection.
t.selection.invert();
```

`t.selection.setAll(true)` selects every row — **including rows currently
hidden by the filter**. For "select all currently-visible rows", copy from the
filter (see below).

## Selected rows as a new DataFrame

```ts
const subset = t.clone(t.selection);
```

Same shape as `t.clone(t.filter)` — different mask:

| Idiom                  | Mask used      | Source DF modified?      |
|------------------------|----------------|--------------------------|
| `t.clone(t.filter)`    | `t.filter`     | no — clone of visible rows |
| `t.clone(t.selection)` | `t.selection`  | no — clone of selected rows |

```datagrok-exec
// Open selected rows as a new table view.
const subset = t.clone(t.selection);
grok.shell.addTableView(subset);
```

If the selection is empty, `clone` returns a zero-row DataFrame (no error).
Guard with `t.selection.anyTrue` first for user-facing flows.

## Cross-skill bridges — selection ↔ filter

| Idiom                                  | Effect                                | Phrasing it answers                                |
|----------------------------------------|---------------------------------------|----------------------------------------------------|
| `t.filter.copyFrom(t.selection)`       | hide non-selected rows in place       | "show only", "filter to the selection"             |
| `t.selection.copyFrom(t.filter)`       | mark every visible row as selected    | "select all currently-visible rows"                |

**Critical disambiguation:**

- "select rows [by / where / top N]" → set `t.selection` via `t.selection.init(pred)` or index-list pattern.
- "extract" / "copy" / "open as new table" / "give me a new table" → `t.clone(t.selection)` (new DataFrame; original unchanged).

`t.selection.copyFrom(t.filter)` fires `onSelectionChanged` but does NOT fire
a filter event; the reverse holds for `t.filter.copyFrom(t.selection)`.

## Current row

```datagrok-exec
// Focus row 5 (the platform usually does this on click; rarely needed in scripts).
if (5 >= 0 && 5 < t.rowCount)
  t.currentRowIdx = 5;
```

`t.currentRow` (a `Row` proxy) and `t.currentCell` (a `Cell`) follow
`currentRowIdx` — read them on demand. Setting `currentRowIdx = -1` is the
valid "no current row" state.

## Describing the current selection

No built-in "describe selection" — hand-roll a structured report: count, capped
index list, sample row from the first selected index.

```datagrok-exec
// Describe the current selection.
const sel = t.selection;
const all = sel.getSelectedIndexes();
const MAX_INDEXES = 50;
const MAX_SAMPLE_COLS = 12;
const indexes = Array.from(all.slice(0, MAX_INDEXES));

const out = {
  count: sel.trueCount,
  total: t.rowCount,
  indexes,
  currentRowIdx: t.currentRowIdx,
  sample: undefined,
};

if (all.length > 0) {
  const firstIdx = all[0];
  const cols = t.columns;
  const sample = {};
  const colCap = Math.min(cols.length, MAX_SAMPLE_COLS);
  for (let c = 0; c < colCap; c++) {
    const col = cols.byIndex(c);
    sample[col.name] = col.get(firstIdx);
  }
  out.sample = sample;
}
return out;
```

## Selection event lifecycle

One event — unlike filter:

| Event                       | When                                                |
|-----------------------------|-----------------------------------------------------|
| `t.onSelectionChanged`      | Any mutation on `t.selection` (`set`, `init`, `setAll`, `and/or/xor/andNot`, `invert`, `copyFrom`). Same observable as `t.selection.onChanged`. |

```datagrok-exec
// Log selection count whenever it changes.
const sub = t.onSelectionChanged.subscribe(() => {
  console.log(`selected: ${t.selection.trueCount}`);
});
view.subs.push(sub);
```

In a viewer / widget, push onto `viewer.subs[]` so the subscription cancels on
detach (`DG.Widget` does this automatically).

**Coalescing.** Every `BitSet` mutator takes `notify: boolean = true`. Pass
`false` in a loop and call `fireChanged()` once at the end:

```ts
for (let i = 0; i < indices.length; i++)
  t.selection.set(indices[i], true, false);  // notify = false
t.selection.fireChanged();                   // one coalesced event
```

## Common patterns

```datagrok-exec
// Toggle selection on a single row.
const i = 5;
t.selection.set(i, !t.selection.get(i));
```

```datagrok-exec
// Top N by a column — sort, slice, set.
const idx = Array.from(t.getSortedOrder(['activity'], [false])).slice(0, 5);
t.selection.setAll(false, false);
for (const i of idx) t.selection.set(i, true, false);
t.selection.fireChanged();
```

```datagrok-exec
// "Select what I see" → narrow → "filter to what I just selected" pipeline.
t.selection.copyFrom(t.filter);                                                // visible → selection
const flag = t.getCol('flag');
t.selection.and(DG.BitSet.create(t.rowCount, i => flag.get(i) === 'ok'));      // narrow
t.filter.copyFrom(t.selection);                                                // push back to filter
```

## Anti-patterns

1. **`t.selection.setAll(true)` to "clear" the selection.** That selects every
   row — the **opposite** of what filter's clear does. Use
   `t.selection.setAll(false)`.
2. **`t.filter.setAll(false)` to "clear" the filter.** Hides every row. Use
   `t.filter.setAll(true)`. The cross-skill pair must be memorized together.
3. **`t.rows.select(rowPred)` for "replace selection" in perf-critical code.**
   Functional, but iterates the `RowList` and per-bit `set` — materially slower
   than `t.selection.init(pred)`, which is buffer-direct.
4. **`t.selection.init(...)` to *add* a condition to the selection.** `init`
   zeros the buffer first. To extend, OR in
   `DG.BitSet.create(t.rowCount, pred)`.
5. **Per-bit `set` in a loop without `notify=false`.** N writes → N events →
   laggy UI. Always `set(i, x, false)` in the loop and `fireChanged()` at the
   end.
6. **`t.selection.setAll(true)` when the user means "select all visible".**
   That includes filtered-out rows. Use `t.selection.copyFrom(t.filter)`.
7. **Combining selections across DataFrames of different row counts.**
   `dfA.selection.and(dfB.selection)` only makes sense if rows are aligned
   1:1. Cross-table sync by key columns is `grok.data.linkTables(...)`.
8. **`t.rows.removeWhereIdx(pred)` when asked to "remove from the selection".**
   That deletes rows from the DF entirely. To narrow the selection, use
   `andNot` or `set(i, false, false)` in a loop.

## Out of scope

- **Filter operations.** `t.filter` lives in `datagrok-filtering`. Same BitSet
  API, **opposite clear polarity** — see the polarity table above.
- **Viewer-driven interactive selection.** Lasso, grid click, shift-extend,
  ctrl-toggle — dispatched via `BitSet.handleClick(rowPred, mouseEvent, ...)`.
  Covered by `datagrok-viewers`.
- **Cross-table selection sync.** `grok.data.linkTables(t1, t2, k1, k2,
  [SYNC_TYPE.SELECTION_TO_SELECTION, ...])` establishes long-lived sync.
- **Generic DataFrame cloning.** For "copy of this DF with these columns" (no
  selection involved), see `datagrok-df-and-columns`.
