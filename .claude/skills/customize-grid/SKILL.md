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
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# customize-grid

## When to use

You have a `TableView` (`grok.shell.addTableView` / `grok.shell.tv`) and
want to alter how its grid renders before the user touches it: "format
`weight` as `#.0000`", "show only these three columns", "color `IC50`
linearly red→green", "sort by a custom comparer".

## Prerequisites

- A package scaffold (`grok create <Name>`) OR an interactive script run
  from **Functions → Scripts** — both use the same JS-API.
- A `TableView` — `grok.shell.addTableView(df)` or `grok.shell.tv`.

## Concepts

Three surfaces with different blast radius (`DG-FACT-108`, `110`,
`115`):

- `view.grid.col(n)` → `GridColumn`: `format`, `width`, `categoryColors`,
  `visible`, `cellType`. Mutates ONE grid.
- `df.col(n)` → `Column`: `valueComparer`, `meta.colors.set*`, tags.
  Mutates every grid / viewer over `df`.
- `df.columns.setOrder([…])` reorders the dataframe;
  `view.grid.columns.setOrder([…])` reorders one grid's layout
  (`DG-FACT-109`).

## Steps

1. **Format numbers and dates on the `GridColumn`** (`DG-FACT-108`,
   `429` — numeric/datetime only).
   ```typescript
   const view = grok.shell.addTableView(grok.data.demo.demog());
   view.grid.col('height')!.format  = 'scientific';   // 1.79e+02
   view.grid.col('weight')!.format  = '#.0000';
   view.grid.col('started')!.format = 'dd.MM.yyyy';
   ```

2. **Reorder columns and set visibility.** Both take case-insensitive
   `string[]`. `setOrder` puts unlisted columns AFTER listed ones —
   doesn't hide them; use `setVisible` to hide (`DG-FACT-110`).
   ```typescript
   view.grid.columns.setOrder(['age', 'sex', 'race']);    // pin to left
   view.grid.columns.setVisible(['age', 'sex', 'race']);  // hide rest
   ```

3. **Filter visible rows with `setRowOrder`.** A partial index list
   HIDES non-listed rows (`DG-FACT-111`).
   ```typescript
   view.grid.setRowOrder([1, 56, 3, 6, 4]);   // 5 rows visible
   ```

4. **Sort all rows with `sortIndexes` or `grid.sort`.** `sortIndexes`
   keeps all rows visible; `grid.sort([], [])` clears active sort
   (`DG-FACT-111`, `116`).
   ```typescript
   const h = df.col('height')!, w = df.col('weight')!;
   view.grid.sortIndexes((i, j) => h.get(i)/w.get(i) - h.get(j)/w.get(j));
   view.grid.sort(['disease', 'weight'], [true, false]);   // multi-sort
   ```

5. **Custom value comparer for non-alphabetic ordering.**
   `column.valueComparer` is on `Column` (not `GridColumn`) — used by
   EVERY viewer/aggregation over `df` once set (`DG-FACT-115`).
   ```typescript
   const months: Record<string, number> = {Jan: 1, Feb: 2, Mar: 3, Apr: 4, May: 5};
   df.col('month')!.valueComparer = (a, b) => months[a] - months[b];
   view.grid.setRowOrder(df.col('month')!.getSortedOrder());
   ```

6. **Resize columns by setting `width`.** `grid.columns.rowHeader` is
   `byIndex(0)`; the first DATAFRAME column is at GridColumn index 1
   (`DG-FACT-117`).
   ```typescript
   view.grid.col('age')!.width        = 200;
   view.grid.columns.byIndex(4).width = 300;
   view.grid.columns.rowHeader.width  = 100;
   ```

7. **Color-code via `column.meta.colors` — never hand-write tags.**
   Enum is `'Off' | 'Categorical' | 'Conditional' | 'Linear'`
   (`DG-FACT-113`, `114`). Use `gridColumn.categoryColors` for a
   GRID-only override (`DG-FACT-108`).
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

- **`grid.col('age')` returns `null`.** Wrong name or grid not yet
  rendered; verify with `df.columns.names()` (`DG-FACT-109`).
- **`setRowOrder` silently drops rows.** Partial index list HIDES
  non-listed rows (`DG-FACT-111`); use `sortIndexes` to keep all.
- **Color literal renders wrong.** Colors are ARGB with alpha in the
  HIGH byte — `0xFF0000FF` is blue, not red. Prefer `DG.Color.red` or
  `DG.Color.fromHtml('#ff8800')` (`DG-FACT-112`).
- **`column.valueComparer` reorders unexpected charts.** Affects EVERY
  viewer over `df` (`DG-FACT-115`); use `grid.sortIndexes` for
  one-grid scope.
- **Tilde-renamed column not visible via `setVisible`.** Name now
  contains the tilde — rename back first (`DG-FACT-430`).

## See also

- Source: `help/develop/how-to/grid/customize-grid.md`.
- Knowledge: `DG-FACT-108`–`117`, `429`, `430`.
- Related skills: `custom-cell-renderers`, `column-tooltip`.
