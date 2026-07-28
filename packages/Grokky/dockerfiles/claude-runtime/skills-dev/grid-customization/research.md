# Grid customization API research — for `datagrok-grid-customization` skill (topic 5)

Citations use `<absolute-path>:<line>`. Path abbreviations:

- `js-api/` = `/Users/aleksashka/Desktop/datagrok/reddata/public/js-api/src/`
- `samples/` = `/Users/aleksashka/Desktop/datagrok/reddata/public/packages/ApiSamples/scripts/`
- `grokky/` = `/Users/aleksashka/Desktop/datagrok/reddata/public/packages/Grokky/`

The grid surface is split across three classes (`Grid`, `GridColumn`, `GridColumnList` in
`js-api/grid.ts`) plus an inherited `IGridSettings` shape (`js-api/interfaces/d4.ts:4-296`).
Unlike the viewer surface, there is essentially no Grokky helper for it today — `grokky.sortRows`
exists but is a **silent no-op** for the grid (returns a sort order, never applies it).
The most important new conceptual content this skill teaches is the *data-side vs grid-side*
split: where to use `col.meta.*` vs `view.grid.col(name).*`.

## 1. Summary table

| # | Area | Key entrypoints | Canonical idiom | Wrapper-worthiness |
|---|------|-----------------|-----------------|--------------------|
| A | Grid access | `TableView.grid` getter (`js-api/views/view.ts:401`), `df.plot.grid()` (off-view, `js-api/dataframe/data-frame.ts:551+`) | `view.grid` when on a `TableView`; document that `view.grid` only exists on `TableView` (not `View`) | **None** — single accessor |
| B | Column visibility, width, order | `grid.col(name)` (`grid.ts:955`), `grid.columns.byName/byIndex` (`grid.ts:836,829`), `gridCol.visible/width/idx` (`grid.ts:714,735,695`), `grid.columns.setVisible(names)` (`grid.ts:850`), `grid.columns.setOrder(names)` (`grid.ts:843`), `gridCol.move(idx)` (`grid.ts:797`) | `view.grid.columns.setVisible(['a','b','c'])`, `view.grid.col('smiles').width = 200` | **High** — `configureGrid` batches show/hide/order/width |
| C | Pinning columns | `gridCol.pin()` / `gridCol.unpin()` (`grid.ts:699-700`); also `PowerGrid:addPinnedColumn` function for fancier pinning (`samples/grid/advanced/pinned-columns.js:7`) | `view.grid.col('SMILES').pin()` | **Medium** — fold into `configureGrid({pin})` |
| D | Display formatting | `col.meta.format` (`column-helpers.ts:222`) — canonical, data-side; **`gridCol.format`** also exists (`grid.ts:723`) and overrides per-grid only | `df.col('IC50').meta.format = '0.00'` (preferred) — both viewers and grid see it | **None** — cross-link to df-and-columns; trap "format the column" → `col.meta.format` |
| E | Sorting the grid | `grid.sort(cols, orders?)` (`grid.ts:1022`), `grid.sortIndexes((a,b)=>...)` (`grid.ts:1030`), `grid.setRowOrder(int32[])` (`grid.ts:1045`), `grid.getRowOrder()` (`grid.ts:1038`); read-only `grid.sortByColumns` / `grid.sortTypes` getters (`grid.ts:990-993`) | `view.grid.sort(['potency'], [false])` — descending | **High** — `applySort` (fixes the no-op), `clearSort` |
| F | Color coding — data side | `col.meta.colors.setLinear/setConditional/setCategorical/setLinearAbsolute/setDisabled` (`column-helpers.ts:75-126`) | `df.col('activity').meta.colors.setLinear(['#f00','#0f0'])` | **Already covered in df-and-columns** — keep, cross-link |
| F2 | Color coding — grid side | `gridCol.categoryColors` (`grid.ts:741`), `gridCol.backColor` (`grid.ts:718`), `gridCol.isTextColorCoded` (`grid.ts:766`); `onCellPrepare` for per-cell (`grid.ts:998`) | Use only when override must affect grid alone (e.g. PowerGrid pinned cols, custom heatmap row) | **Low** — document; data-side is default |
| G | Custom cell renderer (briefly) | `gridCol.cellType` (`grid.ts:727`), `col.meta.cellRenderer` (`column-helpers.ts:237`) | `view.grid.col('mol').cellType = 'html'`; full custom rendering = `create-cell-renderer` skill | **None** — cross-link |
| H | Row height & pinned rows | `IGridSettings.rowHeight` (`d4.ts:42`) via `grid.setOptions({rowHeight: N})`; pinned rows via `grid.setOptions({pinnedRowValues, pinnedRowColumnNames})` (`samples/grid/advanced/pinned-rows.js:5`); `grid.pinnedRows` getter (`grid.ts:1052`); `grid.resetRowHeight()` (`grid.ts:1099`) | `view.grid.setOptions({rowHeight: 100})` | **Low** — optional inside `configureGrid` |
| I | Header customization | `gridCol.name` (`grid.ts:703`) — actually renames the underlying column; `col.meta.friendlyName` (`column-helpers.ts:214`) — visual rename; `gridCol.headerCellStyle.{font,backColor,textColor,horzAlign}` (`grid.ts:710`, `GridCellStyle` `grid.ts:1173-1220`) | `df.col('x').meta.friendlyName = 'Activity (nM)'` to rename without breaking refs | **None** — trap "rename column" → friendlyName |
| J | `grokky.sortRows` bug | `grokky-api.ts:33-36` — returns `df.getSortedOrder(columns, asc)` but never mutates view | **Replace** with `applySort(view, cols, orders)` → `view.grid.sort(cols, orders)` | **Critical fix** |
| K | Clearing / resetting grid | `grid.sort([])` clears sort; `gridCol.visible = true` per-column to reshow; `grid.resetRowHeight()` (`grid.ts:1099`); `view.resetLayout()` (`view.ts:609`) **closes viewers too** | Selective: `clearSort` / `resetGrid` flags | **High** — `clearSort(view)`, `resetGrid(view, {...})` |
| L | Width policies | `DG.ColumnWidthType` enum: `Minimal / Compact / Optimal / Maximal` (`grid.ts:664-669`); `grid.setColumnsWidthType(type, cols?)` (`grid.ts:1156`); `gridCol.setWidthType(type)` (`grid.ts:807`); `grid.autoSize(maxW, maxH, ...)` (`grid.ts:1150`) | `view.grid.setColumnsWidthType(DG.ColumnWidthType.Optimal)` | **Low** — optional inside `configureGrid({widthPolicy})` |

## 2. Grid access (Area A)

`TableView.grid` getter (`js-api/views/view.ts:401`):

```ts
get grid(): Grid { return toJs(api.grok_View_Get_Grid(this.dart)); }
```

It exists **only on `TableView`** — not on `View`. The exec block's `view` global may be a
plain `DG.ViewBase`; the runtime injects it from `grok.shell.v`. Always guard:

```js
if (!(view instanceof DG.TableView)) throw new Error('grid customization requires a TableView');
const grid = view.grid;
```

`Grid` extends `Viewer<IGridSettings>` (`grid.ts:925`), so it inherits `setOptions(map)`,
`getOptions()`, `props.*`, and `props.<setting>` access from `Viewer` (`js-api/viewer.ts:146,158,118`).
That means everything in the `IGridSettings` interface (`js-api/interfaces/d4.ts:4-296`) can be
set via `grid.setOptions({...})` regardless of whether there's a dedicated getter/setter on `Grid`.

Off-view rendering of a grid is available via `df.plot.grid()` (sample: `samples/grid/advanced/auto-size.js:3`).
Not relevant to most Claude exec blocks — they should always operate on `view.grid`. The
`grid.canvas` / `grid.overlay` HTMLCanvasElement getters (`grid.ts:980,985`) exist for custom
rendering but are out of scope here.

## 3. Column visibility, width, order (Area B)

### Lookup

`GridColumnList` (`grid.ts:813-870`):

- `grid.columns` returns the list (`grid.ts:948`).
- `grid.col(name)` is a shortcut for `grid.columns.byName(name)` (`grid.ts:955` → `:836`).
- `grid.columns.byIndex(idx)` — note that index 0 is the **row header** (`grid.ts:822-824`), so
  data columns start at index 1.
- `grid.columns.rowHeader` returns `byIndex(0)`.
- `grid.columns.length` (`grid.ts:856`).

Iteration: no explicit `Symbol.iterator`. Use `for (let i = 1; i < grid.columns.length; i++) { … }`
or `grid.columns.byIndex(i)`. Some samples iterate by name via `df.columns.names()`.

### Visibility

Two paths:

1. **Per-column toggle** — `gridCol.visible = false` (`grid.ts:735-736`).
2. **Bulk "show only these"** — `grid.columns.setVisible(['a','b'])` (`grid.ts:850-852`).
   Hides all others. Sample: `samples/grid/resize/hide-columns.js:6`.

There's also a data-side trick documented in the hide-columns sample (`:9-11`): prefixing the
column name with `~` (`table.columns.byName('sex').name = '~sex'`) hides it **in all views**,
not just this grid. Document, but warn that it renames the underlying column.

### Width

- `gridCol.width = 200` (`grid.ts:714-715`) — pixels.
- `grid.columns.byName('age').width = 200` and `grid.columns.byIndex(4).width = 300` and
  `grid.columns.rowHeader.width = 100` are all canonical (sample `samples/grid/resize/resize-columns.js:5-7`).
- Read-only `gridCol.getDataWidth()` (`grid.ts:800`) gives the pixel width needed to render the
  longest value — useful for "auto-fit this one column".

### Width policies

`DG.ColumnWidthType` enum (`grid.ts:664-669`): `Minimal | Compact | Optimal | Maximal`.

- `grid.setColumnsWidthType(DG.ColumnWidthType.Maximal)` applies to all columns (or pass a
  `columns` array to scope it) (`grid.ts:1156`).
- `gridCol.setWidthType(DG.ColumnWidthType.Minimal)` per-column (`grid.ts:807`).
- Sample: `samples/grid/resize/column-width-type.js:7-10`.

### Reorder

- `grid.columns.setOrder(['sex','age','height'])` (`grid.ts:843`). Columns not listed are
  positioned after the listed ones (sample header `samples/grid/order/order-columns.js:1`).
- `gridCol.move(position)` for single-column reorder (`grid.ts:797`).
- `gridCol.idx` is read-only (`grid.ts:695`) — assigning to it does nothing; use `move`.

## 4. Pinning columns (Area C)

Two API paths:

1. **Built-in** — `gridCol.pin()` / `gridCol.unpin()` (`grid.ts:699-700`).
   Used by `PowerGrid/src/package.ts:316` (`gridCol.pin();`) and
   `PreclinicalCase/src/views/observation-timelines-view.ts:257`. **No** `pin('left'|'right')`
   parameter — pinning is left-only in the built-in API. Multiple columns can be pinned.
2. **PowerGrid plugin** — `await grok.functions.call('PowerGrid:addPinnedColumn', {gridCol})`
   (`samples/grid/advanced/pinned-columns.js:7`). Adds an enhanced pinned column with extra
   rendering. Requires the PowerGrid package installed. Prefer the built-in `pin()` unless the
   user explicitly asks for the PowerGrid features.

The `IGridSettings.frozenColumns: number` (`d4.ts:149`) is a different mechanism — it freezes
the first N columns by their grid position. Use `view.grid.setOptions({frozenColumns: 2})`
to freeze the first two data columns. Distinct from per-column `pin()` and harder to combine.

## 5. Display formatting (Area D)

There are **two `format` accessors** and they behave differently:

- `col.meta.format` (`column-helpers.ts:222`) — sets the `format` tag on the **DataFrame
  column**. Visible to grid, all viewers, CSV export. Canonical.
- `gridCol.format` (`grid.ts:723`) — overrides format **only in this grid** instance.
  Doesn't affect other viewers or exports.

For "Format IC50 with 2 decimal places", use `df.col('IC50').meta.format = '0.00'`.
Use `gridCol.format` only when the grid display should differ from the rest of the platform
(rare). The full format-string vocabulary is documented in `samples/grid/styles/data-format.js:13-74`
(`'int'`, `'two digits after comma'`, `'four digits after comma'`, `'compact'`, `'scientific'`,
`'money'`, `'percent'`; date formats like `'dd.MM.yyyy'`, `'MMM d, yyyy'`, `'relative'`).

This skill should **trap** "format the column" → recommend `col.meta.format` and cross-link to
the df-and-columns skill (see `df-and-columns/research.md:222`).

## 6. Sorting (Area E)

### The blessed path

`grid.sort(columns: string[] | Column[], orders?: boolean[])` (`grid.ts:1022-1025`):

```js
view.grid.sort(['age', 'sex'], [true, false]);  // ascending age, descending sex
view.grid.sort(['age']);                         // ascending (default)
```

`orders` is parallel to `columns`. `true` = ascending, `false` = descending. When omitted, all
ascending. This is the only mutating sort entrypoint on `Grid`.

### Read-back

- `grid.sortByColumns: Column[]` (`grid.ts:990`) — **getter only**, returns the live sort columns.
- `grid.sortTypes: boolean[]` (`grid.ts:993`) — **getter only**.

There is **no setter** for either. Setting `grid.sortByColumns = [...]` silently fails (it's a
plain `get` in TypeScript) — `grid.sort()` is the only way to mutate the sort. The original
research brief assumed setters existed; they do not.

The underlying `IGridSettings` carries `sortByColumnNames: string[]` and `sortTypes: boolean[]`
(`d4.ts:50-52`), so `grid.setOptions({sortByColumnNames: ['age'], sortTypes: [true]})` also
works in principle — but the `grid.sort()` method is shorter, type-safe (accepts `Column` too),
and is the one used by every sample.

### Custom comparer / explicit order

- `grid.sortIndexes((a, b) => bmi(a) - bmi(b))` (`grid.ts:1030`, sample `samples/grid/order/order-rows-by-comparer.js:8`)
  — custom row comparer.
- `grid.setRowOrder([1, 56, 3, 6, 4])` (`grid.ts:1045`, sample `samples/grid/order/order-rows.js:5`)
  — explicit Int32Array-like.
- `grid.getRowOrder()` returns the current `Int32Array` (`grid.ts:1038`).

### Clearing sort

`grid.sort([])` clears the sort (passing an empty array). Confirmed semantically equivalent
because the underlying setter accepts the empty list. Alternative: `grid.setOptions({sortByColumnNames: [], sortTypes: []})`.

### The `grokky.sortRows` no-op

`grokky-api.ts:33-36`:

```ts
function sortRows(df: DG.DataFrame, columns: string[], orders: ('asc' | 'desc')[] = []): Int32Array {
  const asc = columns.map((_, i) => (orders[i] ?? 'asc') === 'asc');
  return df.getSortedOrder(columns, asc);
}
```

`df.getSortedOrder` returns a permutation array — Claude can use it for data-side operations,
but it **does not touch the grid**. Calling `grokky.sortRows(t, ['age'], ['desc'])` looks like a
sort but the grid display is unchanged. Two follow-ups:

1. Add `grokky.applySort(view, cols, orders?)` that calls `view.grid.sort(cols, orders.map(o => o === 'asc'))`.
2. Either keep `sortRows` for the rare case where the caller wants the permutation array
   (rename to `getSortOrder`?), or delete it and tell Claude not to use it.

The skill should explicitly mention this trap.

## 7. Color coding — data side vs grid side (Areas F and F2)

The most important new content in this skill. Two parallel APIs, very different semantics.

### Data-side — `col.meta.colors` (`js-api/dataframe/column-helpers.ts:48-135`)

Already documented in `df-and-columns/research.md:76-90`. Surface:

- `setLinear(range?, {min, belowMinColor, max, aboveMaxColor})` (`column-helpers.ts:75`)
- `setLinearAbsolute({value: hexColor}, options)` (`column-helpers.ts:92`)
- `setCategorical(colorMap?, {fallbackColor, matchType})` (`column-helpers.ts:103`)
- `setConditional({'20-170': '#00FF00'})` (`column-helpers.ts:113`)
- `setDisabled()` (`column-helpers.ts:124`)
- `getType()` returns `DG.COLOR_CODING_TYPE.{LINEAR|CATEGORICAL|CONDITIONAL|OFF}` (`const.ts:797-802`).
- `getColor(i)`, `getColors()` (`column-helpers.ts:128,132`) — read-only.

Sets tags on the column (`.color-coding-type`, `.color-coding-linear`, etc.). **All** viewers
that read those tags pick up the new colors — grid, scatter plot color axis, etc. CSV export
preserves the tags. This is the canonical path.

Sample idiom (`samples/grid/color-coding/color-coding.js:4-13`):

```js
t.col('height').meta.colors.setConditional({'20-170': '#00FF00', '170-190': '#220505'});
t.col('age').meta.colors.setLinear(['#ff0000','#ffff00','#00ff00'], {min: 19, max: 70});
t.col('weight').meta.colors.setLinearAbsolute({58.31:'#73aff5', 137.65:'#ffa500'}, {aboveMaxColor: '#ff0000'});
t.col('race').meta.colors.setCategorical({'Asian': 4278190335, 'Black': 4286578816});
```

### Grid-side — `GridColumn` properties

- `gridCol.categoryColors: {[s: string]: number}` (`grid.ts:741-742`) — per-category colors,
  ARGB integers. Sample: `samples/grid/color-coding/category-colors.js:5-8`:
  ```js
  view.grid.col('sex').categoryColors = { 'M': 0xFF0000FF, 'F': 0xFF800080 };
  ```
  **Grid-only**. Does not affect scatter plot color axis.
- `gridCol.backColor` (`grid.ts:718`) — uniform background ARGB integer for the **whole column**.
- `gridCol.isTextColorCoded` (`grid.ts:766`) — when true, color coding paints the text instead of the cell background.
- `gridCol.headerCellStyle.backColor` (`grid.ts:710`, `GridCellStyle.backColor` `grid.ts:1203`) — column header background.
- Per-cell overrides via `grid.onCellPrepare(cell => cell.style.backColor = ...)` (`grid.ts:998`, sample `samples/grid/events/custom-cell-prepare.js`).

### When to prefer which

| Goal | Use |
|------|-----|
| Color-code a numeric column green→red across all viewers | `col.meta.colors.setLinear(...)` |
| Color-code a string column by category, picked up by scatter plot | `col.meta.colors.setCategorical(...)` |
| Color the **whole** column background uniformly (no value mapping) | `gridCol.backColor = 0xFFFFEEEE` |
| Per-cell rules that should not leak into other viewers (e.g. flag a sketchy value) | `grid.onCellPrepare(cell => …)` |
| Apply category colors only inside one specific grid view (rare) | `gridCol.categoryColors = {...}` |
| Turn coloring off | `col.meta.colors.setDisabled()` (data side) |

The answer is almost always data-side. Grid-side `categoryColors` and `backColor` are useful
when (a) the data column should stay raw for other consumers, or (b) the customization is
specific to this grid instance (e.g. a sub-view). Almost never for "color-code the activity
column" — that should be `col.meta.colors.setLinear`.

## 8. Custom cell renderers (Area G)

Two paths to register a renderer per column:

- `gridCol.cellType = 'Molecule'` (`grid.ts:727-728`) — sets the cell type string; the platform
  looks up a registered renderer by that name. Sample `samples/grid/advanced/charts-in-cells.js:6-7`:
  `tv.grid.columns.byName('ic50').cellType = 'html'`.
- `col.meta.cellRenderer = 'Molecule'` (`column-helpers.ts:237`) — alternative on the data side.

Authoring a new renderer (subclassing `GridCellRenderer` from `grid.ts:1276`) is out of scope —
the repo already has a `create-cell-renderer` skill at `public/.claude/skills/`.

Ad-hoc rendering for single cases: `grid.onCellPrepare(cell => cell.style.element = htmlEl)`
(`grid.ts:998`, sample `samples/grid/advanced/charts-in-cells.js:13`) injects an HTML element
per cell — combined with `grid.setOptions({rowHeight: 200})` and `gridCol.cellType = 'html'`.

## 9. Row height & pinned rows (Area H)

- **Row height**: `view.grid.setOptions({rowHeight: 100})` (samples: `samples/grid/advanced/charts-in-cells.js:11`, `samples/grid/html-cells/html-dynamic-cells.js:5`). Backed by `IGridSettings.rowHeight` (`d4.ts:42`).
- **Reset row height**: `grid.resetRowHeight()` (`grid.ts:1099`).
- **Pin a row to the top**:
  ```js
  view.grid.setOptions({ pinnedRowColumnNames: ['subj', 'subj'], pinnedRowValues: ['1', '0'] });
  ```
  (`samples/grid/advanced/pinned-rows.js:5-7`). Each pinned row is identified by a
  `(columnName, value)` pair — parallel arrays. Backed by `IGridSettings.pinnedRowValues`
  and `pinnedRowColumnNames` (`d4.ts:54-56`).
- Read the current pinned rows: `grid.pinnedRows` (`grid.ts:1052`, `Iterable<number>`).
- `grid.onPinnedRowsChanged` event (`grid.ts:1124`).
- Column header height: `grid.colHeaderHeight` (read-only getter, `grid.ts:1094`) or
  `setOptions({colHeaderHeight: N})` per `d4.ts:31`.

## 10. Header customization (Area I)

- **Display name**: prefer `col.meta.friendlyName = 'Activity (nM)'` (`column-helpers.ts:214`) — visible across all viewers, doesn't rename the underlying column (so existing references to `col.name` still work).
  - `gridCol.name = '...'` (`grid.ts:703-704`) actually **renames** the underlying DF column and breaks references. The user usually wants `friendlyName`.
- **Header background / text**: `gridCol.headerCellStyle.backColor = DG.Color.green; gridCol.headerCellStyle.textColor = 0xFF000000;` (sample `samples/grid/styles/column-cell-styles.js:10`).
- **Header font**: `gridCol.headerCellStyle.font = 'bold 12px Arial'`.
- **Header height**: see §9.
- **Column groups** (visual grouping in headers): `df.meta.setGroups({groupName: {color, columns}})` (sample `samples/grid/advanced/column-groups.js:2-6`). Data-side.

## 11. Clearing / resetting customization (Area K)

| Reset target | How |
|--------------|-----|
| Sort | `grid.sort([])` or `grid.setOptions({sortByColumnNames: [], sortTypes: []})` |
| Row height | `grid.resetRowHeight()` (`grid.ts:1099`) |
| Per-column visibility | `for each grid col: gridCol.visible = true` — no bulk reset |
| Per-column widths | No public reset; assign explicit values, or call `grid.setColumnsWidthType(DG.ColumnWidthType.Optimal)` to apply a default policy (`grid.ts:1156`) |
| Color coding | `col.meta.colors.setDisabled()` (data side, recommended); or remove the tags |
| Custom cell renderer | `gridCol.cellType = 'default'` (or set `col.meta.cellRenderer = null`) |
| Pinned columns | `gridCol.unpin()` (`grid.ts:700`) |
| Pinned rows | `grid.setOptions({pinnedRowColumnNames: [], pinnedRowValues: []})` |
| All viewers (closes them all) | `view.resetLayout()` (`view.ts:609`) — note: **closes added viewers AND resets the grid layout**. Per `view.ts:608`: "Resets view layout, leaving only grid visible." |

There is no single "reset all grid customization to defaults" method. `view.resetLayout()` is
the closest blunt instrument but it also closes scatter plots / histograms the user added.

The skill should expose a selective `resetGrid(view, {visibility?, widths?, sort?, colors?, pinned?})`.

## 12. Common multi-step idioms

From the RFC demo: "Pin SMILES on the left, hide the index, color activity green-red, sort by activity desc, format IC50 with 2 decimals":

```js
const grid = view.grid;
const df = view.dataFrame;
// data side
df.col('Activity').meta.colors.setLinear(['#FF0000', '#FFFF00', '#00FF00']);
df.col('IC50').meta.format = '0.00';
// grid side
grid.col('SMILES').pin();
grid.col('Index').visible = false;
grid.sort(['Activity'], [false]);
```

"Show only these 5 columns":

```js
view.grid.columns.setVisible(['name', 'smiles', 'logD', 'activity', 'ic50']);
```

"Make the molecule column wider":

```js
view.grid.col('molecule').width = 250;
```

"Color-code activity from green to red":

```js
view.dataFrame.col('activity').meta.colors.setLinear(['#00FF00', '#FFFF00', '#FF0000']);
```

## 13. Data-side vs grid-side decision table

| User intent | Use | Lives on | Rationale |
|-------------|-----|----------|-----------|
| Color a column by value | `col.meta.colors.setLinear/setCategorical/setConditional` | `Column` tags | Picked up by all viewers, persists in CSV |
| Change number format | `col.meta.format = '0.00'` | `Column` tags | Applies everywhere |
| Rename for display | `col.meta.friendlyName = '...'` | `Column` tags | Non-breaking |
| Hide a column | `view.grid.col('x').visible = false` or `view.grid.columns.setVisible([...])` | Grid | Display only — column still in df, scatter plot can still use it |
| Reorder columns | `view.grid.columns.setOrder([...])` | Grid | Display only |
| Resize column width | `view.grid.col('x').width = 200` | Grid | Display only |
| Pin column to left | `view.grid.col('x').pin()` | Grid | Display only |
| Sort | `view.grid.sort([col], [asc])` | Grid (sortByColumnNames + sortTypes) | Sort is a grid concept, not a df concept; this is correct |
| Pin a specific row to top | `view.grid.setOptions({pinnedRowValues, pinnedRowColumnNames})` | Grid | Display only |
| Set row height | `view.grid.setOptions({rowHeight: N})` | Grid | Display only |
| Set custom cell renderer for a column | `view.grid.col('x').cellType = 'Molecule'` or `col.meta.cellRenderer = 'Molecule'` | Either works; data-side travels with the column | Prefer data-side when the renderer should follow the column across views |

The principle: **data-side (`col.meta.*`) for anything that semantically belongs to the column;
grid-side (`view.grid.col(x).*`) for view-local presentation only**.

## 14. Proposed `grokky` helpers

All target `grokky-api.ts` and get exported via `export const grokky = {...}`. Type signatures:

```ts
/** Batch grid customization. All fields optional and applied in a stable order:
 *  1. visibility (show wins over hide), 2. order, 3. widths, 4. pin, 5. formats, 6. row height.
 *  Mutates `view.grid` in place; returns nothing. */
function configureGrid(view: DG.TableView, options: {
  hide?: string[];                                 // gridCol.visible = false
  show?: string[];                                 // grid.columns.setVisible([...]) — shows only these
  order?: string[];                                // grid.columns.setOrder([...])
  widths?: Record<string, number>;                 // gridCol.width = px
  widthPolicy?: 'Minimal' | 'Compact' | 'Optimal' | 'Maximal';  // global
  pin?: string[];                                  // pin these to left
  unpin?: string[];                                // unpin these
  formats?: Record<string, string>;                // col.meta.format — data-side!
  rowHeight?: number;                              // grid.setOptions({rowHeight})
  frozenColumns?: number;                          // grid.setOptions({frozenColumns})
}): void;

/** Apply column color coding. Polymorphic: kind picks the helper.
 *  - 'linear' on a non-numeric column throws (caller sees a real error, not silent miscolor).
 *  - 'off' clears all coding regardless of previous kind. */
function colorCode(col: DG.Column, spec:
  | { kind: 'linear'; range?: (string | number)[]; min?: number; max?: number;
      belowMinColor?: string; aboveMaxColor?: string; absolute?: Record<number, string> }
  | { kind: 'categorical'; colors: Record<string, string | number>;
      fallbackColor?: string | number; matchType?: 'exact' | 'regex' }
  | { kind: 'conditional'; rules: Record<string, string | number> }
  | { kind: 'off' }
): void;

/** Apply visual sort to the grid. Fixes the silent no-op in `sortRows`.
 *  `orders[i]` parallel to `columns[i]`; defaults to all 'asc'. */
function applySort(view: DG.TableView, columns: string[],
                   orders?: ('asc' | 'desc' | boolean)[]): void;

/** Clear visual sort. */
function clearSort(view: DG.TableView): void;

/** Selective reset. Each flag opts into clearing that aspect.
 *  - visibility: all columns visible again
 *  - widths: applies the 'Optimal' width policy
 *  - sort: clears sort
 *  - colors: data-side setDisabled() for every column
 *  - pinned: unpin every pinned column; clears pinned rows */
function resetGrid(view: DG.TableView, options?: {
  visibility?: boolean;
  widths?: boolean;
  sort?: boolean;
  colors?: boolean;
  pinned?: boolean;
}): void;
```

Plus existing `sortRows` should be **renamed to `getSortOrder`** and clearly documented as
"returns permutation; does not apply to the grid". Or remove and tell Claude to call
`df.getSortedOrder(...)` directly. The current name is a trap.

`colorCode` mirrors the polymorphic shape of `setColumnMeta`'s `colorCoding` field in the
df-and-columns research (`research.md:283-307`) — keep the shapes identical so both skills
teach the same vocabulary.

## 15. Anti-patterns

1. **Looping rows to color cells**:
   ```js
   for (let i = 0; i < df.rowCount; i++) view.grid.cell('activity', i).style.backColor = ...;
   ```
   Wrong. Use `col.meta.colors.setLinear(...)`. The loop is O(N) and gets blown away on rerender.

2. **`view.grid.sortByColumns = [...]` / `view.grid.sortTypes = [...]`**:
   These are read-only getters (`grid.ts:990,993`). The assignment silently does nothing in
   strict mode is a no-op; otherwise TypeScript will catch it. Use `grid.sort(cols, orders)`.

3. **Calling `grokky.sortRows` expecting visual sort**:
   It returns a `Int32Array` permutation and never touches the grid. The skill must redirect
   to `applySort(view, cols, orders)`.

4. **`gridCol.name = 'new label'` to rename for display**:
   Actually renames the underlying DataFrame column. Use `col.meta.friendlyName` instead.

5. **`gridCol.format` vs `col.meta.format`**:
   The grid-side `format` is rarely what's intended. "Format the IC50 column" means data-side.
   Use `df.col('IC50').meta.format = '0.00'`.

6. **`gridCol.categoryColors` for cross-viewer coloring**:
   Only the grid sees these colors. Scatter plot won't. Use `col.meta.colors.setCategorical`.

7. **`view.resetLayout()` to reset grid only**:
   `resetLayout` also closes added viewers (scatter plots, histograms, filters). Use the
   selective `resetGrid()` helper for grid-only resets.

8. **Mixing `setVisible([...])` with explicit `visible = true` toggles**:
   `setVisible` already hides everything not in the list. Following it with `gridCol.visible = true`
   on a column not in the list re-shows that column but leaves the others hidden — easy
   off-by-one if not careful. Pick one strategy.

9. **`grid.columns.byIndex(0)` to get the first data column**:
   Index 0 is the row header (`grid.ts:822-824`). Use `byIndex(1)` or `byName('actual_name')`.

10. **Forgetting `grid.invalidate()` after manual style mutations**:
    `gridCol.width = 200` and similar trigger a repaint automatically. Per-cell style edits
    inside `onCellPrepare` repaint on the next cycle. But ad-hoc style writes to
    `cell.style.backColor` from outside `onCellPrepare` (e.g. from an event handler) need
    `grid.invalidate()` (`grid.ts:1073`).

## 16. Open questions

1. **Should `applySort` accept `DG.Column[]` as well as `string[]`?** `grid.sort()` accepts both
   (`grid.ts:1022`). Probably yes, mirroring the underlying API.

2. **Should `configureGrid` handle `format` (and other `col.meta.*`) at all, since those are
   data-side?** Arguments for: one-stop ergonomic batch from Claude's POV ("set these grid
   options and these formats together"). Arguments against: muddies the data/grid split this
   skill is teaching. Tentative: yes, but the skill copy makes it explicit that `formats` lives
   on the column, not the grid.

3. **Pinning multiple columns on the left in a specific order**: `gridCol.pin()` doesn't take
   a position. Does the order of `pin()` calls determine left-to-right pinned column order?
   Probably yes (each call appends to the pinned list). Worth verifying once during
   implementation.

4. **`gridCol.pin()` vs `frozenColumns`**: Are they mutually exclusive or can they coexist?
   If both, what wins? Worth a short experiment; document the resolution.

5. **`resetGrid({widths: true})` semantics**: Apply `ColumnWidthType.Optimal` (re-derive from
   data) or restore to whatever was set at construction? `Optimal` is the closest equivalent
   to "the default the platform would have used"; go with that.

6. **Should `colorCode` accept an array of columns?** `grokky.colorCode([col1, col2], spec)`
   for "color these three columns the same way" is a nice ergonomic. The existing helper
   shapes elsewhere (`setColumnMeta`) work on single columns; consistency wins.

7. **Per-cell color coding via `onCellPrepare`** — should the skill teach this? It's powerful
   but easy to abuse (one cell-prepare handler per column, leaks subs). Lean: mention as an
   advanced escape hatch, don't elevate.

8. **`heatmapColors` / `globalColorScaling` flags in `IGridSettings`** (`d4.ts:203-207`) —
   the grid has a heatmap mode. Worth a one-paragraph mention or full coverage? Probably one
   paragraph: it's `grid.setOptions({isHeatmap: true})` (`d4.ts:60,256`), and once in heatmap
   mode `heatmapColors`/`globalColorScaling` apply. Not a Day-1 feature.
