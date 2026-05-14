---
name: datagrok-grid-customization
description: Sort, hide, show, reorder, resize, pin, format, and color-code columns in a Datagrok TableView grid from a datagrok-exec block. Use whenever the user asks to sort by a column (any direction), hide / show / reorder / pin / resize columns, freeze the first N columns, change number-format display, color-code cells (linear / categorical / conditional), set row height, rename a column for display only, or reset the grid back to defaults. Distinct from datagrok-df-and-columns (which owns column-level data metadata like semType, units, friendlyName, and is also where canonical color-coding lives) and from datagrok-viewers (which owns scatter plot / histogram / etc.). Does NOT cover filtering (`datagrok-filtering`), selection (`datagrok-selection`), custom cell renderer authoring (`create-cell-renderer`), saving / restoring layouts, or grid event handlers.
---

# datagrok-grid-customization

Use the `grokky.*` grid helpers inside a `datagrok-exec` block. They wrap
`view.grid`, `grid.columns`, and the per-`GridColumn` setters so Claude doesn't
have to remember that `view.grid.sortByColumns` is a getter-only property
(silent no-op when assigned), that `grid.columns.byIndex(0)` is the row header
(not the first data column), that `gridCol.name = '…'` renames the underlying
DataFrame column (rather than relabeling for display only), or that
`view.resetLayout()` closes every viewer the user added (rather than only
resetting the grid).

## What this skill covers

The grid is the table widget at index 0 of `view.viewers`. It exposes three
classes worth of customization: `Grid` (the view), `GridColumnList`
(`grid.columns`), and `GridColumn` (`grid.col('name')`). This skill handles
sort, visibility, reorder, width, pinning, formatting display, row height,
header rename (display-only), and column color coding. For each capability
it picks the right side of the data-vs-grid split and surfaces the right
`grokky.*` helper.

**Color coding lives on both sides.** The canonical place is data-side
(`col.meta.colors.set*`), because it propagates to *every* viewer that
consumes the column — scatter plot color axis, trellis plot facets, even
CSV exports preserve the tags. This skill provides `grokky.colorCode(col, spec)`
that calls the data-side API; the shape mirrors `setColumnMeta({colorCoding: …})`
in `datagrok-df-and-columns` so both skills teach the same vocabulary. The
df-and-columns skill remains the deeper reference for color-coding semantics
on a column; this skill is the practical "while customizing the grid, here's
how you color" path. Grid-only color (e.g. uniform `gridCol.backColor` for
one specific table view) is documented as a rare advanced override.

**Out of scope.** Row filtering / selection (separate skills). Custom
cell-renderer authoring (the `create-cell-renderer` skill). Layout save /
restore (`view.saveLayout()` / `view.loadLayout()`). Grid event handlers
(`grid.onCellPrepare`, `grid.onCellRender`) beyond a one-line pointer.
Cross-table linking. Filter / selection / current-row sync between viewers
(automatic via the shared DataFrame).

## Quick reference

| Helper                                             | One-liner                                                                                |
|----------------------------------------------------|------------------------------------------------------------------------------------------|
| `grokky.applySort(view, columns, orders?)`         | Sort the grid. Replaces the silent no-op `sortRows`. `orders[i]` matches `columns[i]`.   |
| `grokky.clearSort(view)`                           | Clear any sort applied to the grid.                                                      |
| `grokky.configureGrid(view, opts)`                 | Batch hide / show / order / widths / pin / formats / rowHeight / widthPolicy / frozenCols.|
| `grokky.colorCode(col, spec)`                      | Data-side color coding (linear / categorical / conditional / off). Throws on bad input.  |
| `grokky.resetGrid(view, opts?)`                    | Selective reset. Defaults clear visibility / widths / sort / colors. **Keeps viewers.**  |
| `grokky.pinColumn(view, name)`                     | Pin one column on the left. Convenience wrapper over `gridCol.pin()`.                    |
| `grokky.unpinColumn(view, name)`                   | Unpin one column.                                                                        |

Globals inside every `datagrok-exec` block: `grok`, `ui`, `DG`, `view`, `t`
(the current `DG.DataFrame`, when the view is a `TableView`), `grokky`.

## The data-side vs grid-side decision

The single most useful mental model for this skill: every customization
either lives on the **column** (visible to all viewers, CSV-preserved) or on
the **grid column** (display-only, this view only). Use this table before
reaching for an API:

| Need                                                            | Side       | API                                                     |
|-----------------------------------------------------------------|------------|---------------------------------------------------------|
| Color cells in a way that all viewers respect                   | data       | `grokky.colorCode(col, …)` (calls `col.meta.colors.*`)  |
| Color cells in the GRID ONLY (e.g. uniform column tint)         | grid       | `gridCol.backColor = 0xFFFFEEEE` (rare; advanced)       |
| Format numbers as 2 decimals everywhere                         | data       | `setColumnMeta(col, {format: '0.00'})` (df-and-columns) |
| Format display in the grid only                                 | grid       | `configureGrid({formats: {…}})` sugar over `gridCol.format` |
| Hide a column in the grid (still in DF)                         | grid       | `configureGrid({hide: […]})` / `gridCol.visible = false`|
| Drop a column from the DF                                       | data       | `removeColumns(df, …)` (df-and-columns)                 |
| Display name change visible to ALL viewers                      | data       | `setColumnMeta(col, {friendlyName: 'Potency'})` (df-and-columns) |
| Display name change in this grid only                           | (no API)   | use `friendlyName`. **Never `gridCol.name = …`.**       |
| Sort the grid                                                   | grid       | `grokky.applySort(view, cols, orders)`                  |
| Pin a column to the left                                        | grid       | `grokky.pinColumn(view, name)` or `gridCol.pin()`       |
| Set column width                                                | grid       | `configureGrid({widths: {name: px}})`                   |

The principle: **data-side (`col.meta.*`) for anything that semantically
belongs to the column; grid-side (`view.grid.col(name).*`) for this view's
presentation only.**

> **Critical trap.** `gridCol.name = 'Potency'` is **not** a display rename —
> it renames the underlying DataFrame column. Every existing reference by
> the old name silently breaks. To relabel for display, use
> `setColumnMeta(col, {friendlyName: 'Potency'})` from the df-and-columns
> skill. There is no grid-local rename.

## Quick wins / demo recipes

The headline use cases. Each recipe is a single `datagrok-exec` block.

### Color-code activity from green to red

```datagrok-exec
grokky.colorCode(t.col('activity'), {kind: 'linear', range: ['#00FF00', '#FFFF00', '#FF0000']});
```

Data-side. All viewers reading `activity` now see the same gradient.
Scatter-plot color axis bound to `activity` picks it up; CSV export
preserves the tags.

### Sort by potency descending

```datagrok-exec
grokky.applySort(view, ['potency'], ['desc']);
```

Replaces the long form `view.grid.sort(['potency'], [false])`. The wrong
form `view.grid.sortByColumns = ['potency']` is a silent no-op — see the
sorting section.

### Customize a whole HitTriage-style hit list in one block

> "Pin SMILES on the left, hide the index column, color activity red-to-green,
> sort by activity desc."

```datagrok-exec
grokky.colorCode(t.col('activity'), {kind: 'linear', range: ['#FF0000', '#FFFF00', '#00FF00']});
grokky.configureGrid(view, {
  hide: ['index'],
  pin:  ['SMILES'],
});
grokky.applySort(view, ['activity'], ['desc']);
```

Three logical steps, one block. The order matters only in the sense that
`pin` / `hide` mutate visible-column positions, and the sort is independent.

## Sorting the grid

`grokky.applySort` is the canonical entrypoint. `orders[i]` is parallel to
`columns[i]`; missing entries default to `'asc'`. Boolean `true` is treated
as `'asc'`.

| Need                          | Call                                                            |
|-------------------------------|-----------------------------------------------------------------|
| Sort by one column asc        | `grokky.applySort(view, ['Activity'])`                          |
| Sort by one column desc       | `grokky.applySort(view, ['Activity'], ['desc'])`                |
| Multi-column                  | `grokky.applySort(view, ['Class', 'Activity'], ['asc', 'desc'])`|
| Clear sort                    | `grokky.clearSort(view)`                                        |

```datagrok-exec
// Sort by class ascending, then by activity descending — two-level sort.
grokky.applySort(view, ['class', 'activity'], ['asc', 'desc']);
```

```datagrok-exec
// Clear any sort.
grokky.clearSort(view);
```

### Anti-patterns

**`view.grid.sortByColumns = ['x']` silently does nothing.** The `sortByColumns`
and `sortTypes` properties on `Grid` are getter-only (returning the live sort
state). Assigning to them does **not** apply a sort — neither does it throw,
which makes the bug invisible. Use `grokky.applySort` or call the underlying
`view.grid.sort(cols, orders)` directly.

**`grokky.sortRows` (legacy) returns a permutation array but never touches
the grid.** It exists for the rare case where Claude needs the sorted-index
`Int32Array` for data-side work. For visual sort, use `grokky.applySort`.

```ts
// WRONG — silent no-op:
view.grid.sortByColumns = [view.dataFrame.col('Activity')];

// RIGHT:
grokky.applySort(view, ['Activity'], ['desc']);

// Also right (underlying API; the helper is just thinner):
view.grid.sort(['Activity'], [false]);
```

## Column visibility & order

```ts
grokky.configureGrid(view, {
  hide?:  string[];     // gridCol.visible = false for each
  show?:  string[];     // restricts visible set to exactly these (others hidden)
  order?: string[];     // grid.columns.setOrder(...) — listed names go first
});
```

Three independent fields. `hide` toggles per-column; `show` is the bulk
"show only these"; `order` reorders the listed names to the front while
unlisted columns retain their relative position behind them.

```datagrok-exec
// Hide a single column.
grokky.configureGrid(view, {hide: ['_internalId']});
```

```datagrok-exec
// Show only these five columns; hide everything else.
grokky.configureGrid(view, {show: ['SMILES', 'MW', 'LogP', 'activity', 'class']});
```

```datagrok-exec
// Re-order: SMILES first, then MW, LogP — everything else keeps relative order.
grokky.configureGrid(view, {order: ['SMILES', 'MW', 'LogP']});
```

Missing column names are skipped with a `console.warn`, never thrown — Claude
will see the warning in the runtime log but the block does not abort.

### Show every column again

`grokky.resetGrid(view, {visibility: true})` walks `grid.columns` and sets
`visible = true` on each. Equivalent: a `for` loop over named columns.
**Do not** use `grid.columns.byIndex(0)` to "reset the first column" — that
index is the row header.

## Column widths

```datagrok-exec
// Width by name → pixels.
grokky.configureGrid(view, {widths: {SMILES: 250, MW: 80, LogP: 70}});
```

To apply an auto-fit policy instead of explicit pixel widths, use
`widthPolicy`:

```datagrok-exec
// Auto-size every column to its content (Optimal). Other values: 'Minimal' | 'Compact' | 'Maximal'.
grokky.configureGrid(view, {widthPolicy: 'Optimal'});
```

`widthPolicy` covers all columns at once (mirrors
`grid.setColumnsWidthType(...)`). For per-column auto-fit, set the width
explicitly with `widths`. The policies map to `DG.ColumnWidthType` —
`'Minimal'` is just enough for the label, `'Compact'` clamps narrower,
`'Optimal'` matches the widest visible value, `'Maximal'` expands to
remaining space.

## Pinning

`gridCol.pin()` / `gridCol.unpin()` work directly — no PowerGrid plugin
needed. Pin is **left-only** — there's no positional argument. To pin
multiple columns, pin them in left-to-right order; each call appends to the
pinned list.

```datagrok-exec
grokky.pinColumn(view, 'SMILES');
```

```datagrok-exec
// Pin SMILES then ID — SMILES ends up leftmost.
grokky.configureGrid(view, {pin: ['SMILES', 'ID']});
```

```datagrok-exec
// Unpin.
grokky.unpinColumn(view, 'SMILES');
```

To **freeze the first N visible columns regardless of name**, use
`frozenColumns` (a different mechanism from per-column `pin()`; freezes by
grid position):

```datagrok-exec
grokky.configureGrid(view, {frozenColumns: 2});
```

Pinning and `frozenColumns` are independent. Mixing them works but gets
visually confusing; prefer one strategy per view.

## Formatting display values

`configureGrid({formats: …})` is sugar over `col.meta.format` per column.
For a one-column format change in pipeline with column metadata work,
`setColumnMeta(col, {format: '0.00'})` (from the df-and-columns skill) is
the canonical path. The two have the same effect — they both write the
`format` tag on the column, which the grid honors.

```datagrok-exec
// Two decimals on IC50; integer with thousands separator on count.
grokky.configureGrid(view, {formats: {IC50: '0.00', count: '#,##0'}});
```

Available format strings (data-side; documented in `samples/grid/styles/data-format.js`):

- Numbers: `'int'`, `'0.00'`, `'0.0000'`, `'#,##0.00'`, `'compact'`,
  `'scientific'`, `'money'`, `'percent'`.
- Dates: `'dd.MM.yyyy'`, `'MMM d, yyyy'`, `'relative'`.

> Grid-only override: `gridCol.format = '0.00'` exists and overrides the
> column's format *for this grid only*. Use only when the grid should
> display differently from other consumers of the same column (rare).
> The wrapper sets the data-side `col.meta.format` to keep behavior
> consistent across viewers.

## Color coding

`grokky.colorCode(col, spec)` is the single entry. Spec is polymorphic on
`kind`. Data-side: writes `col.meta.colors.*`, which propagates to every
viewer reading the column.

### Linear (numeric)

```datagrok-exec
// Green → yellow → red across the value range.
grokky.colorCode(t.col('activity'), {kind: 'linear', range: ['#00FF00', '#FFFF00', '#FF0000']});
```

```datagrok-exec
// Pin the gradient endpoints to specific values, with overflow colors.
grokky.colorCode(t.col('IC50'), {
  kind: 'linear',
  range: ['#FF0000', '#FFFF00', '#00FF00'],
  min: 0, max: 100,
  belowMinColor: '#000080', aboveMaxColor: '#000000',
});
```

Linear on a non-numeric column throws — mirrors `setColumnMeta`. The error
text names the offending column and its type.

### Categorical (string)

```datagrok-exec
// Class column with explicit per-category colors. Hex strings or ARGB ints both work.
grokky.colorCode(t.col('class'), {
  kind: 'categorical',
  colors: {'A': '#FF0000', 'B': '#00FF00', 'C': '#0000FF'},
  fallbackColor: '#CCCCCC',
});
```

### Conditional (numeric or string, range-based)

```datagrok-exec
// Conditional rules. Keys are range or value patterns.
grokky.colorCode(t.col('height'), {
  kind: 'conditional',
  rules: {'20-170': '#00FF00', '170-190': '#FFFF00', '190-': '#FF0000'},
});
```

### Off

```datagrok-exec
grokky.colorCode(t.col('activity'), {kind: 'off'});
```

`{kind: 'off'}` calls `col.meta.colors.setDisabled()` — clears any previous
coding regardless of kind.

> **Grid-only color (advanced).** `gridCol.backColor = 0xFFFFEEEE` paints a
> uniform background for the whole column in this grid only. Use only when
> the customization explicitly should not leak to other viewers. For per-cell
> styling rules, see `grid.onCellPrepare` (out of scope for this skill —
> brief mention only).

## Row height & per-row styling

```datagrok-exec
// Tall rows — useful when rendering molecules or HTML in cells.
grokky.configureGrid(view, {rowHeight: 100});
```

To reset row height, `view.grid.resetRowHeight()` (or `resetGrid(view, …)`).

For per-row / per-cell custom rendering, use `grid.onCellRender` or
`grid.onCellPrepare` — both are out of scope for this skill. See the
`create-cell-renderer` skill for proper per-column renderers (e.g. a
molecule cell). Quick advanced patch:

```ts
view.grid.onCellPrepare((gc) => {
  if (gc.tableColumn?.name === 'activity' && gc.cell.value > 100)
    gc.style.backColor = 0xFFFFCCCC;
});
```

This pattern is grid-only, per-cell, and easy to abuse — for value-based
color, prefer `grokky.colorCode` (data-side).

## Resetting / clearing

`grokky.resetGrid(view, opts?)` is selective. Default opts:
`{visibility: true, widths: true, sort: true, colors: true}`. Pass `false`
to opt out of any aspect.

| Flag         | What it does                                                             |
|--------------|--------------------------------------------------------------------------|
| `visibility` | Sets `gridCol.visible = true` for every grid column (re-shows hidden).    |
| `widths`     | Applies `DG.ColumnWidthType.Optimal` globally (re-derives from content).  |
| `sort`       | Calls `clearSort(view)`.                                                  |
| `colors`     | Calls `col.meta.colors.setDisabled()` for every data column.              |

```datagrok-exec
// Reset everything in the grid back to defaults.
grokky.resetGrid(view);
```

```datagrok-exec
// Reset only the sort.
grokky.resetGrid(view, {sort: true, visibility: false, widths: false, colors: false});
```

```datagrok-exec
// Restore visibility but keep widths and colors.
grokky.resetGrid(view, {visibility: true, widths: false, sort: false, colors: false});
```

`resetGrid` does **not** close other viewers (scatter plots, histograms,
filter panels). To close every non-grid viewer, use the `datagrok-viewers`
skill's `grokky.closeAllViewers(view)`. To do both — reset the grid AND
close all other viewers — the canonical platform call is `view.resetLayout()`,
but it cannot be selective. Confirm with the user before calling
`view.resetLayout()` if they only asked for "reset the grid".

## Common pitfalls / anti-patterns

1. **`view.grid.sortByColumns = [...]` — silent no-op.**
   The setter doesn't exist (getter-only). Use `grokky.applySort` or
   `view.grid.sort(cols, orders)`.

2. **`gridCol.name = 'New Label'` for display rename — destructive.**
   This renames the underlying DataFrame column. Every reference by the
   old name breaks. Use `setColumnMeta(col, {friendlyName: 'New Label'})`
   from the df-and-columns skill.

3. **`grid.columns.byIndex(0)` returns the row header, not the first data
   column.** Use `grid.columns.byIndex(1)` for "first data column" or
   `grid.col('name')` / `grid.columns.byName('name')` whenever possible.

4. **Looping over `df.rows` (or `view.grid.cell(...)`) to apply cell
   background colors** — fragile, slow, and gets blown away on the next
   repaint. Use `grokky.colorCode(col, …)` instead.

5. **`view.resetLayout()` when the user just wanted to clear sort.** This
   closes every viewer the user added (scatter plot, histogram, filter
   panel). Use `grokky.clearSort(view)` (sort only) or
   `grokky.resetGrid(view, {…})` (grid only) instead.

6. **`gridCol.categoryColors = {…}` to color a column "by class" — grid-only.**
   Scatter plot color axis bound to the same column will *not* pick up these
   colors. Use `grokky.colorCode(col, {kind: 'categorical', colors: {…}})`
   so the coloring propagates.

7. **`gridCol.format` for "format the column".** The grid-only override is
   rarely what's intended. Use the data-side `col.meta.format` (via
   `configureGrid({formats: …})` or `setColumnMeta`).

8. **`grid.columns.setVisible([...])` followed by `gridCol.visible = true`
   on a column not in the list** — easy off-by-one. `setVisible` already
   hides everything not listed. Pick one strategy.

9. **`view.grid` on a non-`TableView`.** Script views, function views, and
   bare `ViewBase` instances have no `grid`. Guard:
   ```ts
   if (!(view instanceof DG.TableView)) throw new Error('grid customization requires a TableView');
   ```
   The wrappers do this themselves (throw with a clear message).

10. **Calling `grokky.sortRows` expecting visual sort.** That helper returns
    an `Int32Array` permutation; it does not touch the grid. Use
    `grokky.applySort` for visual sort. (The legacy helper stays for the
    rare data-side use case.)

## Out of scope

- **Custom cell renderer authoring** — subclassing `DG.GridCellRenderer`,
  `meta.cellRenderer`-driven renderers. Covered by the `create-cell-renderer`
  skill. From here, `gridCol.cellType = 'Molecule'` is the single-line way
  to attach a registered renderer.
- **Viewer configuration** — scatter plots, histograms, line charts, etc.
  All in `datagrok-viewers`.
- **Filtering and selection** — `view.getFiltersGroup()` and `df.selection`
  live in `datagrok-filtering` and `datagrok-selection`.
- **Saving and restoring layouts** — `view.saveLayout()` /
  `view.loadLayout()`, viewer-state snapshots. Future layout skill.
- **Per-cell event handlers and onCellPrepare deep dive** — powerful but
  easy to abuse. One-line pointer only.
- **Pinned rows (`grid.setOptions({pinnedRowValues, pinnedRowColumnNames})`)** —
  not commonly requested; add only if a clear use case emerges. Mention
  briefly under row height if needed.
- **Heatmap mode (`grid.setOptions({isHeatmap: true})`)** — sets the whole
  grid into a heatmap rendering mode. Niche; out of scope for the day-1
  surface.
