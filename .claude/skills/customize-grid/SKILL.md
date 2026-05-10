---
name: customize-grid
description: Programmatically customize a Datagrok grid — column format, order, visibility, width, sorting, and color-coding — from a package or script
---

# customize-grid

## When to use

You have a `TableView` open and want to alter how its grid renders without
asking the user to click through column menus. Triggers: "format `weight`
as `#.0000`", "show only these three columns", "color-code `IC50` linearly
red→green", "make the row header wider", "sort by a custom comparer".

## Prerequisites

- A package scaffold (`grok create <Name>`) OR an interactive script run
  from **Functions → Scripts** — both call the same JS-API.
- `datagrok-api` imports — the article omits them and won't run as written:
  ```typescript
  import * as grok from 'datagrok-api/grok';
  import * as DG   from 'datagrok-api/dg';
  ```
- A `TableView` reference — either `grok.shell.addTableView(df)` or the
  active one (`grok.shell.tv`).

## Steps

1. **Pick the right level: GridColumn vs. Column vs. DataFrame.**
   Customizations split across three surfaces with different blast radius
   (`DG-FACT-108`, `DG-FACT-110`, `DG-FACT-115`):
   - `view.grid.columns.byName(n)` / `view.grid.col(n)` → `GridColumn`.
     Hosts `format`, `width`, `categoryColors`, `visible`, `cellType`.
     Mutates ONE grid only.
   - `df.col(n)` → `Column`. Hosts `valueComparer`, `meta.colors.set*`,
     and tags. Mutates every grid/viewer over `df`.
   - `df.columns.setOrder([…])` reorders columns inside the dataframe
     itself; `view.grid.columns.setOrder([…])` reorders only the visual
     layout. There is NO `view.grid.cols(…)` — the plural is
     `view.grid.columns` (the `GridColumnList`) (`DG-FACT-109`).

2. **Format numbers and dates on the GridColumn.**
   `format` is a string written to the GridColumn (`DG-FACT-108`,
   `js-api/src/grid.ts:723-724`). Applies only to numeric and datetime
   columns; other types ignore the tag.
   ```typescript
   const view = grok.shell.addTableView(grok.data.demo.demog());
   view.grid.col('height')!.format  = 'scientific';   // 1.79e+02
   view.grid.col('weight')!.format  = '#.0000';
   view.grid.col('started')!.format = 'dd.MM.yyyy';
   ```
   Expected: cell text reflows; underlying values untouched, so sorting
   and filtering produce the same result.

3. **Reorder columns and set visibility.**
   Both methods take a `string[]` (case-insensitive). `setOrder` columns
   NOT listed appear AFTER the listed ones — they are not hidden
   (`DG-FACT-110`). To hide, use `setVisible`.
   ```typescript
   view.grid.columns.setOrder(['age', 'sex', 'race']);    // pin to left
   view.grid.columns.setVisible(['age', 'sex', 'race']);  // hide rest
   ```
   Production reference — PowerPack pairs the two when building a custom
   edit grid (`packages/PowerPack/src/dialogs/formula-lines.ts:338-339`).

4. **Reorder rows: `setRowOrder` filters, `sortIndexes` preserves.**
   Critical asymmetry (`DG-FACT-111`, `js-api/src/grid.ts:1027-1048`):
   - `grid.setRowOrder(indexes)` displays ONLY the listed rows — others
     are hidden.
   - `grid.sortIndexes((a, b) => …)` builds a full identity array of size
     `df.rowCount`, sorts via your comparer, applies the result — all
     rows remain visible.
   - `grid.sort(columns, orders?)` sorts by names with optional `boolean[]`
     direction (`true`=ASC, `false`=DESC). `grid.sort([], [])` clears
     the active sort (`DG-FACT-116`).
   ```typescript
   view.grid.setRowOrder([1, 56, 3, 6, 4]);                    // filter
   const h = df.col('height')!, w = df.col('weight')!;
   view.grid.sortIndexes((i, j) => h.get(i)/w.get(i) - h.get(j)/w.get(j));
   view.grid.sort(['disease', 'weight'], [true, false]);       // multi
   ```

5. **Set a custom value comparer for non-alphabetic ordering.**
   `column.valueComparer` is a setter on the underlying `Column`
   (`DG-FACT-115`, `js-api/src/dataframe/column.ts:417-418`), not on
   `GridColumn`. Once set, the comparer is used by EVERY viewer,
   aggregation, group-by, and chart category ordering.
   ```typescript
   const months = {Jan: 1, Feb: 2, Mar: 3, Apr: 4, May: 5};
   df.col('month')!.valueComparer = (a, b) => months[a] - months[b];
   view.grid.setRowOrder(df.col('month')!.getSortedOrder());
   ```
   `getSortedOrder()` returns an `Int32Array` of row indexes; pair with
   `setRowOrder` only when you intend to reduce/reorder the row set.

6. **Resize columns (and the row header) by setting `width`.**
   `width: number` on `GridColumn` (`DG-FACT-108`, `:712-715`).
   `grid.columns.rowHeader` IS the same object as `byIndex(0)` — the
   row header is the first GridColumn (`DG-FACT-117`); the first
   DATAFRAME column lives at GridColumn index 1.
   ```typescript
   view.grid.col('age')!.width        = 200;
   view.grid.columns.byIndex(4).width = 300;
   view.grid.columns.rowHeader.width  = 100;
   ```
   For a responsive layout, recompute on resize — see PowerPack's
   `ResizeObserver` pattern (`formula-lines.ts:354-365`).

7. **Color-code via `column.meta.colors` — never hand-write tags.**
   The `column.meta.colors.set*` helpers (`DG-FACT-114`,
   `js-api/src/dataframe/column-helpers.ts:48-135`) write the canonical
   `.color-coding-type` / `.color-coding-{categorical,conditional,linear}`
   tags for you. The article's all-caps `COLOR_CODING_TYPE` is the
   JS-API constant name, NOT the wire tag (`DG-FACT-DRIFT-046`).
   `COLOR_CODING_TYPE` has exactly four values — `'Off'`, `'Categorical'`,
   `'Conditional'`, `'Linear'` (`DG-FACT-113`). Type constraints —
   numerical: Off/Linear/Conditional; categorical: Off/Categorical;
   datetime: Off/Linear.
   ```typescript
   df.col('site')!.meta.colors.setCategorical({'New York': DG.Color.orange});
   df.col('age')!.meta.colors.setLinear([DG.Color.orange, DG.Color.green]);
   df.col('weight')!.meta.colors.setConditional({
     '<100':    DG.Color.green,
     '100-200': '#ff0000',           // hex strings also accepted
   });
   df.col('height')!.meta.colors.setLinear();   // no args → default palette
   df.col('age')!.meta.colors.setDisabled();    // turn coloring off
   ```
   For a GRID-only categorical override (paint cells in one grid without
   touching DataFrame tags), use `gridColumn.categoryColors`
   (`DG-FACT-108`, `:741-742`):
   ```typescript
   view.grid.col('sex')!.categoryColors = {
     'M': DG.Color.blue,             // prefer named constants — see step 8
     'F': DG.Color.fromHtml('#800080'),
   };
   ```

8. **Use `DG.Color` constants, not raw 32-bit literals.**
   Datagrok colors are 4-byte ARGB with alpha in the high byte
   (`DG-FACT-112`, `js-api/src/color.ts:24`). The article's `0xFF0000FF`
   is `DG.Color.blue`, NOT red — readers expecting RGB-with-alpha-suffix
   will misread it (`DG-FACT-DRIFT-045`). Prefer `DG.Color.red` /
   `.blue` / `.green` / `DG.Color.fromHtml('#ff8800')`.

## Common failure modes

- **`view.grid.col('age')` returns `null` / `TypeError`.** Wrong name,
  column not in the DataFrame, or grid not yet rendered. `grid.col(name)`
  returns `GridColumn | null` (`DG-FACT-109`); verify via
  `df.columns.names()` before using the non-null assertion (`!`).
- **`setRowOrder` silently drops rows.** A partial index list HIDES
  non-listed rows — that is the documented behavior (`DG-FACT-111`).
  If the goal was "reorder all rows", use `sortIndexes` or pass
  `column.getSortedOrder()` (size = `df.rowCount`).
- **Color literal renders the wrong color.** `0xFF0000FF` is BLUE — alpha
  is the HIGH byte (`DG-FACT-112`, `DG-FACT-DRIFT-045`). Use
  `DG.Color.red` / `.fromHtml('#ff0000')` instead of raw hex.
- **Color-coding never fires after writing tags by hand.** Constant names
  (`COLOR_CODING_TYPE`) are NOT the on-the-wire tag keys (those are
  `.color-coding-type`) (`DG-FACT-DRIFT-046`). Always go through
  `column.meta.colors.set*`.
- **`column.valueComparer` reorders unexpected charts.** Once set, the
  comparer is used by EVERY visualization and aggregation over that
  column (`DG-FACT-115`). To reorder ONE grid only, sort the grid via
  `grid.sortIndexes` / `grid.sort` instead of mutating the column.

## Verification

- After each `view.grid.col('X')!.width = N` / `.format = '…'` /
  `.categoryColors = {…}`, the grid re-renders the targeted column
  (no `grid.invalidate()` required for these setters).
- `column.meta.colors.getType()` returns one of `'Off' | 'Categorical' |
  'Conditional' | 'Linear'` — confirms the helper wrote the type tag
  (`DG-FACT-113`).
- Open a second `TableView` over the same `DataFrame`: GridColumn changes
  (`format`, `width`, `categoryColors`) do NOT propagate; Column /
  `meta.colors` changes DO.

## See also

- Source articles:
  - `help/develop/how-to/grid/customize-grid.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-108` … `DG-FACT-117` and drifts
  `DG-FACT-DRIFT-045`, `DG-FACT-DRIFT-046`.
- Reference packages:
  - `packages/PowerPack/src/dialogs/formula-lines.ts:338-365` —
    pairs `setVisible` + `setOrder`, per-column `width`, responsive
    layout via `ResizeObserver`.
  - `packages/Admetica/src/utils/admetica-utils.ts:126-156` — routes a
    model's coloring spec to `meta.colors.setLinear` / `setConditional`
    and toggles `gridColumn.isTextColorCoded` for text-vs-fill rendering.
  - `packages/Charts/src/viewers/group-analysis/group-analysis-viewer.ts:315`
    — `grid.columns.setOrder` driven by a derived list.
- Related skills:
  - `custom-cell-renderers` (sibling — render the *contents* of a grid
    cell; this skill targets the column metadata around it).
  - `column-tooltip` (sibling — register a tooltip widget for a column).
