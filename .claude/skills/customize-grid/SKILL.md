---
name: customize-grid
version: 0.3.1
description: |
  Programmatically alter how a Datagrok `TableView` renders its dataset
  — column format, order, visibility, width, sort/row order, and
  per-column / per-category coloring — without making the user click
  through column menus. For package functions and scripts that open
  a view and arrange it on load.
  Use when asked to "color-code rows in a table view by category",
  "format numeric column values programmatically", or "pin specific
  rows to the top of a dataframe view".
triggers:
  - format column values from code
  - auto-color rows by category
  - set column widths in code
  - hide columns programmatically
  - sort rows by a custom comparer
  - reorder dataframe columns
allowed-tools:
  - Read
  - Edit
  - Bash
harness-authored: true
---

# customize-grid

## When to use

You have a `TableView` (`grok.shell.addTableView` / `grok.shell.tv`) and
want to alter how its grid renders before the user touches it: "format
`weight` as `#.0000`", "show only these three columns", "color `IC50`
linearly red→green", "sort by a custom comparer".

## Prerequisites

- A package scaffold (`grok create <Name>`) OR an interactive script run
  from **Functions → Scripts** — both use the same JS-API.
- `datagrok-api` imports (article omits them):
  `import * as grok from 'datagrok-api/grok'; import * as DG from 'datagrok-api/dg';`
- A `TableView` — `grok.shell.addTableView(df)` or `grok.shell.tv`.

## Concepts

Three surfaces with different blast radius — pick the level first
(`DG-FACT-108`, `DG-FACT-110`, `DG-FACT-115`):

- `view.grid.col(n)` → `GridColumn`: `format`, `width`, `categoryColors`,
  `visible`, `cellType`. Mutates ONE grid.
- `df.col(n)` → `Column`: `valueComparer`, `meta.colors.set*`, tags.
  Mutates every grid / viewer over `df`.
- `df.columns.setOrder([…])` reorders the dataframe;
  `view.grid.columns.setOrder([…])` reorders one grid's layout. Plural
  is `view.grid.columns` (a `GridColumnList`) — NO `grid.cols(…)`
  (`DG-FACT-109`).

## Steps

1. **Format numbers and dates on the `GridColumn`.**
   `format` is a string on `GridColumn` (`DG-FACT-108`,
   `js-api/src/grid.ts:723-724`). Applies only to numeric and datetime
   columns; other types ignore the tag (`DG-FACT-429`).
   ```typescript
   const view = grok.shell.addTableView(grok.data.demo.demog());
   view.grid.col('height')!.format  = 'scientific';   // 1.79e+02
   view.grid.col('weight')!.format  = '#.0000';
   view.grid.col('started')!.format = 'dd.MM.yyyy';
   ```
   Expected: cell text reflows; underlying values untouched, so sorting
   and filtering produce the same result either way.

2. **Reorder columns and set visibility.**
   Both methods take a `string[]` (case-insensitive). `setOrder` columns
   NOT listed appear AFTER the listed ones — they are not hidden
   (`DG-FACT-110`). To hide, use `setVisible`. Production reference —
   PowerPack pairs the two when building an edit grid
   (`packages/PowerPack/src/dialogs/formula-lines.ts:338-339`).
   ```typescript
   view.grid.columns.setOrder(['age', 'sex', 'race']);    // pin to left
   view.grid.columns.setVisible(['age', 'sex', 'race']);  // hide rest
   ```

3. **Filter visible rows with `setRowOrder`.**
   `setRowOrder(indexes)` displays ONLY the rows whose indexes appear in
   the list — a partial list HIDES the rest (`DG-FACT-111`,
   `js-api/src/grid.ts:1045-1048`). Use when you intend to reduce the
   visible row set; pair with `column.getSortedOrder()` to reorder and
   reduce in one call.
   ```typescript
   view.grid.setRowOrder([1, 56, 3, 6, 4]);   // 5 rows visible
   ```

4. **Sort all rows with `sortIndexes` or `grid.sort`.**
   `sortIndexes(comparer)` builds a full identity array of size
   `df.rowCount`, sorts via the comparer, applies it — all rows remain
   visible (`DG-FACT-111`, `js-api/src/grid.ts:1027-1035`).
   `grid.sort(cols, orders?)` sorts by column names with optional
   `boolean[]` direction (`true`=ASC, `false`=DESC); `grid.sort([], [])`
   clears the active sort (`DG-FACT-116`, `packages/Chem/src/package.ts:2004`).
   ```typescript
   const h = df.col('height')!, w = df.col('weight')!;
   view.grid.sortIndexes((i, j) => h.get(i)/w.get(i) - h.get(j)/w.get(j));
   view.grid.sort(['disease', 'weight'], [true, false]);   // multi-sort
   ```

5. **Set a custom value comparer for non-alphabetic ordering.**
   `column.valueComparer` is a setter on the underlying `Column`, NOT on
   `GridColumn` (`DG-FACT-115`, `js-api/src/dataframe/column.ts:417-418`).
   Once set, it's used by EVERY viewer, aggregation, group-by, and chart
   over `df`. `getSortedOrder()` returns an `Int32Array` of row indexes.
   ```typescript
   const months: Record<string, number> = {Jan: 1, Feb: 2, Mar: 3, Apr: 4, May: 5};
   df.col('month')!.valueComparer = (a, b) => months[a] - months[b];
   view.grid.setRowOrder(df.col('month')!.getSortedOrder());
   ```

6. **Resize columns (and the row header) by setting `width`.**
   `width: number` on `GridColumn` (`DG-FACT-108`,
   `js-api/src/grid.ts:712-715`). `grid.columns.rowHeader` IS the same
   object as `byIndex(0)` — the row header is the first GridColumn
   (`DG-FACT-117`); the first DATAFRAME column lives at GridColumn
   index 1. For responsive layout, recompute widths on resize — see
   PowerPack's `ResizeObserver` pattern (`formula-lines.ts:354-365`).
   ```typescript
   view.grid.col('age')!.width        = 200;
   view.grid.columns.byIndex(4).width = 300;
   view.grid.columns.rowHeader.width  = 100;
   ```

7. **Color-code via `column.meta.colors` — never hand-write tags.**
   The `column.meta.colors.set*` helpers (`DG-FACT-114`,
   `js-api/src/dataframe/column-helpers.ts:48-135`) write the canonical
   `.color-coding-type` / `.color-coding-{categorical,conditional,linear}`
   tags. The enum has exactly four values — `'Off'`, `'Categorical'`,
   `'Conditional'`, `'Linear'` (`DG-FACT-113`); constraints — numerical:
   Off / Linear / Conditional; categorical: Off / Categorical; datetime:
   Off / Linear. For a GRID-only categorical override (paint cells in
   one grid without touching DataFrame tags), use
   `gridColumn.categoryColors` (`DG-FACT-108`, `js-api/src/grid.ts:741-742`).
   ```typescript
   df.col('site')!.meta.colors.setCategorical({'New York': DG.Color.orange});
   df.col('age')!.meta.colors.setLinear([DG.Color.orange, DG.Color.green]);
   df.col('weight')!.meta.colors.setConditional({
     '<100': DG.Color.green, '100-200': '#ff0000',   // hex also accepted
   });
   df.col('height')!.meta.colors.setLinear();   // no args → default palette
   df.col('age')!.meta.colors.setDisabled();    // turn coloring off
   view.grid.col('sex')!.categoryColors = {     // grid-only override
     'M': DG.Color.blue, 'F': DG.Color.fromHtml('#800080'),
   };
   ```

## Common failure modes

- **`view.grid.col('age')` returns `null` → `TypeError`.** Wrong name,
  column not in the DataFrame, or grid not yet rendered. `grid.col(name)`
  returns `GridColumn | null` (`DG-FACT-109`); verify via
  `df.columns.names()` before using `!`.
- **`setRowOrder` silently drops rows.** A partial index list HIDES
  non-listed rows (`DG-FACT-111`). If the goal was "reorder all rows",
  use `sortIndexes` or pass `column.getSortedOrder()` (length =
  `df.rowCount`).
- **Color literal renders the wrong color.** Datagrok colors are 4-byte
  ARGB with alpha in the HIGH byte (`DG-FACT-112`, `js-api/src/color.ts:30`).
  `0xFF0000FF` is `DG.Color.blue`, NOT red — prefer `DG.Color.red` /
  `DG.Color.fromHtml('#ff8800')` over raw integer literals.
- **`column.valueComparer` reorders unexpected charts.** Once set, the
  comparer is used by EVERY visualization and aggregation over that
  column (`DG-FACT-115`). To reorder ONE grid only, use `grid.sortIndexes`
  / `grid.sort` instead of mutating the column.
- **Tilde-renamed column not coming back via `setVisible`.** Renaming
  `data.columns.byName('age').name = '~age'` hides at the dataframe
  level; `grid.columns.setVisible(['~age'])` does NOT undo it — the
  column's *name* now contains the tilde (`DG-FACT-430`). Rename back
  to `'age'` first.

## Verification

- `GridColumn` setters (`width`, `format`, `categoryColors`) repaint the
  column without an explicit `grid.invalidate()`.
- `column.meta.colors.getType()` returns one of
  `'Off' | 'Categorical' | 'Conditional' | 'Linear'` (`DG-FACT-113`).
- Open a second `TableView` over the same `DataFrame`: `GridColumn`
  changes do NOT propagate; `Column` / `meta.colors` changes DO.

## See also

- Source articles: `help/develop/how-to/grid/customize-grid.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-108` … `DG-FACT-117`, plus `DG-FACT-429`, `DG-FACT-430`.
- Reference packages:
  - `packages/PowerPack/src/dialogs/formula-lines.ts:338-365` — pairs
    `setVisible` + `setOrder`, per-column `width`, responsive layout
    via `ResizeObserver`.
  - `packages/Chem/src/package.ts:2004-2005` — `grid.sort([], [])` to
    clear the active sort, then `setRowOrder(idxCol.toList())` to apply.
- Related skills: `custom-cell-renderers` (render the *contents* of
  a grid cell), `column-tooltip` (header hover preview).
