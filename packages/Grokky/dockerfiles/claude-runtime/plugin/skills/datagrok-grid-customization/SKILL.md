---
name: datagrok-grid-customization
description: Sort, hide, show, reorder, resize, pin, format, and color-code columns in a Datagrok TableView grid from a datagrok-exec block. Use whenever the user asks to sort by a column (any direction), multi-sort, hide / show / reorder / pin / resize columns, freeze the first N columns, change number-format display, color-code cells (defaults and grid-only tint here; full per-type reference in datagrok-df-and-columns), set row height, or reset the grid back to defaults. Distinct from datagrok-df-and-columns (which owns column-level data metadata like semType, units, friendlyName, and is also where canonical color-coding lives) and from datagrok-viewers (which owns scatter plot / histogram / etc.). Does NOT cover filtering (`datagrok-filtering`), selection (`datagrok-selection`), custom cell renderer authoring (`create-cell-renderer`), saving / restoring layouts, or grid event handlers.
---

# datagrok-grid-customization

Customize a Datagrok grid inside a `datagrok-exec` block using `view.grid`,
`grid.columns`, the per-`GridColumn` setters, and the data-side
`col.meta.colors` / `col.meta.format`.

Globals inside every `datagrok-exec` block: `grok`, `ui`, `DG`, `view`, `t`
(the active `DG.DataFrame` when `view` is a `TableView`).

Color coding is canonically data-side (`col.meta.colors.set*`) because it
propagates to every viewer that reads the column. `datagrok-df-and-columns`
owns the broader column-metadata picture.

## Quick reference

| Need                          | API                                                                |
|-------------------------------|--------------------------------------------------------------------|
| Sort (single / multi)         | `view.grid.sort(['name', 'mw'], [true, false])`                    |
| Clear sort                    | `view.grid.sort([], [])`                                           |
| Hide one column               | `view.grid.col('smiles').visible = false`                          |
| Show only these               | `view.grid.columns.setVisible(['name', 'mw', 'logp'])`             |
| Reorder columns               | `view.grid.col('name').move(1)` per column (see *Visibility and order*) |
| Width (one column, pixels)    | `view.grid.col('smiles').width = 250`                              |
| Width policy (all columns)    | `view.grid.setColumnsWidthType(DG.ColumnWidthType.Optimal)`        |
| Pin to left                   | `view.grid.col('smiles').pin()`                                    |
| Unpin                         | `view.grid.col('smiles').unpin()`                                  |
| Freeze first N visible        | `view.grid.setOptions({frozenColumns: 2})`                         |
| Row height                    | `view.grid.setOptions({rowHeight: 28})`                            |
| Color coding (all types)      | `col.meta.colors.set*` — no args = platform default; full reference in `datagrok-df-and-columns` |
| Turn off color coding         | `col.meta.colors.setDisabled()`                                    |
| Number format                 | `col.meta.format = '0.00'`                                         |

## Critical traps

1. **`grid.sortByColumns = [...]` is a getter-only assignment.** Silent no-op.
   Always call `grid.sort(columnNames, ascendingBooleans)` instead. The booleans
   are `true = asc`, `false = desc`. The pair must be parallel.

2. **`grid.columns.byIndex(0)` returns the row header,** not the first data
   column. Loops start at index 1. Prefer `grid.col('name')` whenever you know
   the name.

3. **`gridCol.name = 'New Label'` renames the underlying DataFrame column** —
   not a display rename. Every reference by the old name silently breaks.
   For display-only rename, use `col.meta.friendlyName` (df-and-columns).

4. **`view.resetLayout()` closes every viewer the user added.** Don't use it
   for "clear sort" or "reset the grid". Reset selectively (see *Resetting*).

## Data-side vs grid-side

| Need                                                         | Side  | API                                                         |
|--------------------------------------------------------------|-------|-------------------------------------------------------------|
| Color cells so every viewer respects it                      | data  | `col.meta.colors.setLinear/setCategorical/setConditional` (df-and-columns) |
| Color the grid column only (uniform tint)                    | grid  | `gridCol.backColor = 0xFFFFEEEE` (rare)                     |
| Format numbers everywhere                                    | data  | `col.meta.format = '0.00'`                                  |
| Format only in this grid                                     | grid  | `gridCol.format = '0.00'`                                   |
| Hide a column (still in DF)                                  | grid  | `gridCol.visible = false`                                   |
| Drop a column from the DF                                    | data  | `df.columns.remove(col)` (df-and-columns)                   |
| Display name across all viewers                              | data  | `col.meta.friendlyName = 'Potency'` (df-and-columns)        |
| Sort the grid                                                | grid  | `grid.sort(cols, bools)`                                    |
| Pin / freeze                                                 | grid  | `gridCol.pin()` or `grid.setOptions({frozenColumns: N})`    |

Principle: **data-side (`col.meta.*`) for what semantically belongs to the
column; grid-side (`view.grid.*`) for this view's presentation only.**

## Sorting

`grid.sort(columnNames, ascendingBooleans)`. Both arrays are parallel.
`true` is ascending, `false` is descending.

```datagrok-exec
// Single column descending.
view.grid.sort(['activity'], [false]);
```

```datagrok-exec
// Multi-column: class asc, then activity desc.
view.grid.sort(['class', 'activity'], [true, false]);
```

```datagrok-exec
// Clear the sort.
view.grid.sort([], []);
```

**Anti-pattern (silent failure):**

```ts
// WRONG — getter-only, does nothing, throws nothing.
view.grid.sortByColumns = [t.col('activity')];
```

`grid.sortByColumns` and `grid.sortTypes` read the live sort state. Assigning
to them is a no-op. The only way to apply a sort is `grid.sort(...)`.

## Visibility and order

```datagrok-exec
// Hide a single column.
view.grid.col('_internalId').visible = false;
```

```datagrok-exec
// Show only these (hides everything else).
view.grid.columns.setVisible(['SMILES', 'MW', 'LogP', 'activity', 'class']);
```

```datagrok-exec
// Re-order: SMILES first, then MW, then LogP. Use gridCol.move(position) per
// column. Index 0 is the row header, so the first data slot is 1. Apply the
// list in order: each move(i) inserts at that position.
view.grid.col('SMILES').move(1);
view.grid.col('MW').move(2);
view.grid.col('LogP').move(3);
```

**Why not `grid.columns.setOrder([...])`?** It exists but does NOT move the listed
columns to the front — the listed columns retain their existing positions while
their *relative* order among each other is enforced (e.g., ApiTests' grid setOrder
test calls `setOrder(['race', 'age'])` and asserts race lands at `byIndex(4)`, not
`byIndex(1)`). For "put these columns first, in this order," use `move()` as above.

To re-show every hidden column, loop `grid.columns` from index 1 (index 0 is
the row header):

```datagrok-exec
const cols = view.grid.columns;
for (let i = 1; i < cols.length; i++) {
  const gc = cols.byIndex(i);
  if (gc) gc.visible = true;
}
```

## Widths

```datagrok-exec
// Per-column pixel widths.
view.grid.col('SMILES').width = 250;
view.grid.col('MW').width = 80;
```

```datagrok-exec
// Width policy for the whole grid. Values: Minimal | Compact | Optimal | Maximal.
view.grid.setColumnsWidthType(DG.ColumnWidthType.Optimal);
```

`Minimal` fits just the label; `Compact` clamps narrower; `Optimal` matches
the widest visible value; `Maximal` expands to fill remaining space.

## Pinning and freezing

`gridCol.pin()` pins on the left. No positional argument — pin in
left-to-right order. `gridCol.unpin()` reverses it.

```datagrok-exec
view.grid.col('SMILES').pin();
view.grid.col('ID').pin();   // SMILES ends up leftmost.
```

`frozenColumns` freezes the first N visible columns by position (different
mechanism from `pin()`; do not mix them in one view):

```datagrok-exec
view.grid.setOptions({frozenColumns: 2});
```

## Row height

```datagrok-exec
// Tall rows — useful when rendering molecules or HTML in cells.
view.grid.setOptions({rowHeight: 100});
```

## Number / date format

`col.meta.format` writes a format tag the grid honors:

```datagrok-exec
t.col('IC50').meta.format = '0.00';
t.col('count').meta.format = '#,##0';
```

Common format strings: `'int'`, `'0.00'`, `'0.0000'`, `'#,##0.00'`,
`'compact'`, `'scientific'`, `'money'`, `'percent'`; dates: `'dd.MM.yyyy'`,
`'MMM d, yyyy'`, `'relative'`.

`gridCol.format = '0.00'` exists too — grid-only override. Use the data-side
`col.meta.format` unless the grid should display differently from other
consumers.

## Color coding

The canonical reference — all four shapes (linear / categorical /
conditional / off), options, and traps — lives in the
**`datagrok-df-and-columns`** skill; load it for anything beyond the
rules below.

- No specific colors or rules requested → call the setter with **no
  arguments** (`setLinear()` / `setCategorical()` / `setConditional()`) —
  the platform's default scheme applies. Never invent a palette the user
  didn't ask for.
- Turn off: `col.meta.colors.setDisabled()`.

Grid-only color (advanced): `view.grid.col('foo').backColor = 0xFFFFEEEE`
paints a uniform background in this grid only. Other viewers ignore it.
For value-based coloring, always use `col.meta.colors.set*` — a scatter plot
color axis on the same column will not pick up `gridCol.backColor` /
`gridCol.categoryColors`.

## Resetting the grid

There is no single "reset grid" API. Combine the relevant primitives:

```datagrok-exec
// Restore visibility (skip index 0 = row header), reset width policy, clear sort, drop color codings.
const cols = view.grid.columns;
for (let i = 1; i < cols.length; i++) {
  const gc = cols.byIndex(i);
  if (gc) gc.visible = true;
}
view.grid.setColumnsWidthType(DG.ColumnWidthType.Optimal);
view.grid.sort([], []);
for (const c of t.columns.toList())
  c.meta.colors.setDisabled();
```

This is the right pattern for "reset the grid customization but keep my
viewers". Do **not** call `view.resetLayout()` — it closes scatter plots,
histograms, filter panels the user added.

## Guarding for non-TableView

`view.grid` is only present on `DG.TableView`. Script views and bare
`ViewBase` have no grid:

```ts
if (!(view instanceof DG.TableView)) return;
```

## Out of scope

- **Custom cell renderer authoring** — see `create-cell-renderer`. From here,
  `gridCol.cellType = 'Molecule'` attaches a registered renderer.
- **Viewer configuration** — `datagrok-viewers`.
- **Filtering and selection** — `datagrok-filtering`, `datagrok-selection`.
- **Layout save / restore** — `view.saveLayout()` / `view.loadLayout()`.
- **`grid.onCellPrepare` / `grid.onCellRender`** — powerful but abuse-prone.
  One-line pointer:
  ```ts
  view.grid.onCellPrepare((gc) => {
    if (gc.tableColumn?.name === 'activity' && gc.cell.value > 100)
      gc.style.backColor = 0xFFFFCCCC;
  });
  ```
