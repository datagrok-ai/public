---
name: custom-cell-renderers
version: 0.1.0
description: |
  Register a per-cell paint routine on the Datagrok grid so the platform
  replaces the default text/number display with a domain-specific
  visualization — stars for a rating, a 2D molecule for SMILES, an SVG
  image, a colored swatch, or a sparkline summarising several numeric
  columns. For plugin authors who want the grid to bind a visualization
  automatically whenever a matching column appears, producing a
  `DG.GridCellRenderer` subclass with its decorator declaration
  codegen-wired into `package.g.ts`.
  Use when asked to "show each compound row as a 2D structure image",
  "paint colored stars in a rating column", or "draw a mini bar chart
  in every row of a summary column".
triggers:
  - per-row mini visualization
  - replace default column display
  - show molecule depiction in grid
  - paint custom shape per row
  - summary chart column
  - bind visual to column type
allowed-tools:
  - Read
  - Edit
  - Write
  - Bash
harness-authored: true
---

# custom-cell-renderers

## When to use

Your package needs the grid to draw something other than the default
text/number cell — stars for a rating, a colored swatch, an SVG molecule,
a pie/bar/radar sparkline, a hyperlink badge — and you want the platform
to bind the renderer automatically whenever a matching column appears.

## Prerequisites

- A package scaffold (`grok create <Name>`); run from the package root.
- `datagrok-api` imports (`* as DG`, `* as grok` — article omits these).
- `datagrok-tools` shipping `@grok.decorators.cellRenderer` (current 6.x).
- Familiarity with `CanvasRenderingContext2D` — `render` paints directly
  to a 2D canvas; nothing in the DOM.

## Steps

1. **Subclass `DG.GridCellRenderer`; override the three required hooks.**
   The platform discovers renderers via the function role `cellRenderer`
   (camelCase — `DG-FACT-099`). Pick a registration form before writing
   the stub — both compile to the same `package.g.ts` wrapper, and the
   article calls the class form "the recommended form" (`custom-cell-renderers.md:22-23`):
   - **Class decorator** — `@grok.decorators.cellRenderer({...})` on the
     class; no factory in `package.ts` (`DG-FACT-103`). Canonical for
     plain renderers (`StarsCellRenderer`, `SvgCellRenderer`) and for
     the article's chart-style `PieChartCellRenderer` example.
   - **Function decorator** — static factory on `PackageFunctions` with
     `@grok.decorators.func({meta:{role:'cellRenderer', cellType:'…',
     virtual:'true'}, tags:['cellRenderer'], outputs:[{type:'grid_cell_renderer',
     name:'result'}]})` (`DG-FACT-104`). Stylistic alternative when you
     keep all renderers in `PackageFunctions` — and the only form that
     accepts the `gridChart` meta flag (see Step 2 and `DG-FACT-428`).

   Do NOT paste the article's `//name:` / `//meta.role:` block into
   `package.ts` — that is the auto-emitted `package.g.ts` shape, not
   the authoring surface. `name`, `cellType`, and `render` all throw
   in the base class (`DG-FACT-102`); optional overrides —
   `defaultWidth/Height`, `renderSettings(gridColumn)`, mouse/key
   handlers (`onClick`, `onMouseEnter`, …),
   `hasContextValue`/`getContextValue`, `clip`.

   ```typescript
   // src/cell-types/stars-cell-renderer.ts — class-decorator form
   import * as DG from 'datagrok-api/dg';
   import * as grok from 'datagrok-api/grok';

   @grok.decorators.cellRenderer({
     name: 'Stars', cellType: 'Stars',
     // @ts-ignore — `tags` is on FuncOptions, not CellRendererOptions
     tags: ['cellRenderer'],
   })
   export class StarsCellRenderer extends DG.GridCellRenderer {
     get name()     { return 'Stars'; }
     get cellType() { return 'Stars'; }
     render(g: CanvasRenderingContext2D, x: number, y: number,
            w: number, h: number, gridCell: DG.GridCell,
            cellStyle: DG.GridCellStyle) {
       const value = Math.round(gridCell.cell.value);
       g.fillStyle = '#FFB400';
       g.fillText('★'.repeat(value), x + 4, y + h / 2 + 4);
     }
   }
   ```
   The `@ts-ignore tags` pattern is a PowerGrid habit
   (`stars-cell-renderer.ts:40`) — redundant under modern codegen, which
   auto-emits `//tags: cellRenderer` regardless (`DG-FACT-DRIFT-043`).

2. **Bind the renderer to columns via `cellType` (or `columnTags`).**
   The platform invokes your renderer when a column's effective cell
   type matches `meta.cellType`. Either set `column.semType = 'Stars'`
   (e.g. via a `semanticTypeDetector`) or use
   `meta.columnTags: 'units=kg,foo=bar'` for tag-based matching
   (`DG-FACT-106`; see `packages/PowerGrid/src/package.g.ts:185-201`).

   For SUMMARY/VIRTUAL renderers (sparklines, pie chart, smart form —
   columns synthesised at render time from several numeric columns),
   set `meta.virtual` so the column appears in the *Add column → New
   chart* menu (`DG-FACT-105`). The article's PieChart example uses
   `virtual: true` on the class decorator and nothing else
   (`custom-cell-renderers.md:14-21`); reach for the function form
   only if you also need the `gridChart` meta flag, which is absent
   from `CellRendererOptions` (`DG-FACT-428`):

   ```typescript
   // src/package.ts — function-decorator form when `gridChart` is needed
   import * as grok from 'datagrok-api/grok';
   import {PieChartCellRenderer} from './sparklines/piechart';
   export class PackageFunctions {
     @grok.decorators.func({
       meta: {role: 'cellRenderer', cellType: 'piechart',
              gridChart: 'true', virtual: 'true'},
       tags: ['cellRenderer'],
       outputs: [{type: 'grid_cell_renderer', name: 'result'}],
       name: 'Pie Chart',
     })
     static piechartCellRenderer() { return new PieChartCellRenderer(); }
   }
   ```
   `meta` flags use string literals on the function-decorator form
   (`virtual: 'true'`, `gridChart: 'true'`) but booleans on the class
   form (`virtual: true`) — a generic typing asymmetry across all
   Datagrok decorators (`DG-FACT-430`).

3. **Build, publish, and let the codegen wire it up.**
   ```bash
   npm install
   npm run build       # FuncGeneratorPlugin emits/updates src/package.g.ts
   grok check          # exits 0
   grok publish <host> # add --release once stable
   ```
   `package.g.ts` is meant to be committed (`DG-FACT-107`). No explicit
   `DG.GridCellRenderer.register(...)` call — the function role IS the
   registration. Expected: each renderer produces a wrapper in
   `src/package.g.ts` with `//tags: cellRenderer`,
   `//output: grid_cell_renderer <name>`, `//meta.role: cellRenderer`,
   `//meta.cellType: <Type>` (compare
   `packages/PowerGrid/src/package.g.ts:57-64` for the plain form and
   `:116-125` for the virtual/`gridChart` form).

## Common failure modes

- **Renderer never fires; cells fall back to default text.** Inspect
  `src/package.g.ts`: each entry MUST contain `//meta.role: cellRenderer`
  (camelCase, case-sensitive) and `//meta.cellType: <Type>` —
  `CellRenderer` / `cell-renderer` (the column-tag key at
  `js-api/src/const.ts:330`) won't register (`DG-FACT-099`).
- **Build emits no entry in `package.g.ts`.** `FuncGeneratorPlugin` is
  not wired into `webpack.config.js` (`DG-FACT-107`). Add
  `const FuncGeneratorPlugin = require('datagrok-tools/plugins/func-gen-plugin');`
  and `new FuncGeneratorPlugin({outputPath: './src/package.g.ts'})`.
- **`'cellType'`/`'name'`/`'Not implemented'` thrown at render time.**
  The subclass forgot to override one of `get name()`, `get cellType()`,
  or `render(...)` — all three throw in `DG.GridCellRenderer` defaults
  (`DG-FACT-102`).
- **TS error: `'tags' does not exist on type CellRendererOptions`.**
  Suppress with `// @ts-ignore` above the `tags` entry, or drop the
  line entirely — codegen auto-emits it either way (`DG-FACT-DRIFT-043`).
- **Sparkline never appears in *Add column → New chart*.** Renderer
  is registered without `meta.virtual` (`DG-FACT-105`). Add
  `virtual: true` (class form) or `virtual: 'true'` (function form);
  if you also need `gridChart: 'true'`, switch to the function form
  because the flag is absent from `CellRendererOptions`
  (`DG-FACT-428`, `DG-FACT-430`).

## Verification

- `npm run build` and `grok check` both exit `0`; `grok publish <host>`
  succeeds.
- Regenerated `src/package.g.ts` contains, per renderer, a wrapper with
  `//tags: cellRenderer`, `//output: grid_cell_renderer <name>`,
  `//meta.role: cellRenderer`, `//meta.cellType: <Type>` (compare
  `packages/PowerGrid/src/package.g.ts:13-64` plain, `116-125` virtual).
- In Datagrok, open a table whose `semType`/tags match `cellType` —
  grid paints via your `render`. For a virtual renderer, *Add column →
  New chart* must list `<Type>`. Switching to a column with a different
  `semType` falls back to the built-in renderer (binding is scoped).

## See also

- Source: `<help_develop_root>/how-to/grid/custom-cell-renderers.md`
  (mirror: `<harness_root>/docs/_internal/articles-mirror/how-to/grid/custom-cell-renderers.md`).
- Knowledge: `<harness_root>/docs/_internal/knowledge/knowledge-graph.md`
  — `DG-FACT-099`…`107`, `DG-FACT-428`, `DG-FACT-429`, `DG-FACT-430`,
  `DG-FACT-DRIFT-043`.
- Reference packages:
  - `packages/PowerGrid/src/cell-types/stars-cell-renderer.ts:37-43` —
    minimal class-decorator form.
  - `packages/PowerGrid/src/cell-types/svg-cell-renderer.ts:6-12` —
    class form with image cache + `onDoubleClick`.
  - `packages/PowerGrid/src/sparklines/piechart.ts:192-332` +
    `packages/PowerGrid/src/package.ts:85-98` — function-decorator form
    (virtual sparkline) with `renderSettings` + `getContextValue`.
  - `packages/PowerGrid/src/package.g.ts:13-201` — auto-emitted wrappers.
- Related skills: `register-identifiers` (registers the `semType` that
  `cellType` matches), `column-tooltip` (same role pattern, different role).
