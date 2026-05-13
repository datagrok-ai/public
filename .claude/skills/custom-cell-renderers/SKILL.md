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
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# custom-cell-renderers

## When to use

Your package needs the grid to draw something other than the default
text/number cell — stars for a rating, a colored swatch, an SVG molecule,
a pie/bar/radar sparkline, a hyperlink badge — and you want the platform
to bind the renderer automatically whenever a matching column appears.

## Prerequisites

- `datagrok-tools` shipping `@grok.decorators.cellRenderer` (current 6.x).
- Familiarity with `CanvasRenderingContext2D` — `render` paints directly
  to a 2D canvas; nothing in the DOM.

## Steps

1. **Subclass `DG.GridCellRenderer`; override the three required hooks.**
   Platform discovers renderers via the function role `cellRenderer`
   (camelCase — `DG-FACT-099`). Two equivalent registration forms:
   - **Class decorator** — `@grok.decorators.cellRenderer({...})` on
     the class (`DG-FACT-103`). Recommended for plain renderers.
   - **Function decorator** — static factory on `PackageFunctions`
     (`DG-FACT-104`). Required when you need the `gridChart` meta flag
     (`DG-FACT-466`).

   The three required overrides — `get name()`, `get cellType()`,
   `render(...)` — all throw in the base class (`DG-FACT-102`).
   Optional: `defaultWidth/Height`, `renderSettings(gridColumn)`,
   mouse/key handlers, `hasContextValue`/`getContextValue`, `clip`.

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
   The `@ts-ignore tags` line is redundant under modern codegen
   (`DG-FACT-DRIFT-043`).

2. **Bind the renderer to columns via `cellType` (or `columnTags`).**
   Set `column.semType = 'Stars'` (via a `semanticTypeDetector`) or
   `meta.columnTags: 'units=kg,foo=bar'` for tag-based matching
   (`DG-FACT-106`).

   For SUMMARY/VIRTUAL renderers (sparklines, pie chart, smart form),
   set `meta.virtual` so the column appears in *Add column → New
   chart* (`DG-FACT-105`). Use the function form if you also need
   `gridChart` (absent from `CellRendererOptions` — `DG-FACT-466`):

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
   `meta` flags are string literals on the function form
   (`virtual: 'true'`) but booleans on the class form (`virtual:
   true`) — typing asymmetry across all decorators (`DG-FACT-468`).

3. **Build, publish, and let the codegen wire it up.**
   ```bash
   npm install
   npm run build       # FuncGeneratorPlugin emits/updates src/package.g.ts
   grok check          # exits 0
   grok publish <host> # add --release once stable
   ```
   `package.g.ts` is meant to be committed (`DG-FACT-107`). The
   function role IS the registration — no explicit `register(...)`.

## Common failure modes

- **Renderer never fires.** `//meta.role: cellRenderer` (camelCase) or
  `//meta.cellType:` missing in `package.g.ts` (`DG-FACT-099`).
- **Build emits no entry in `package.g.ts`.** `FuncGeneratorPlugin`
  not wired into `webpack.config.js` (`DG-FACT-107`).
- **`'Not implemented'` thrown at render time.** Subclass missing one
  of `get name()`, `get cellType()`, `render(...)` (`DG-FACT-102`).
- **TS error: `'tags' does not exist on type CellRendererOptions`.**
  Suppress with `// @ts-ignore` or drop the line — codegen auto-emits
  it (`DG-FACT-DRIFT-043`).
- **Sparkline never appears in *Add column → New chart*.** Missing
  `meta.virtual` (`DG-FACT-105`); use function form if you also need
  `gridChart` (`DG-FACT-466`).

## See also

- Source: `help/develop/how-to/grid/custom-cell-renderers.md`.
- Knowledge: `DG-FACT-099`–`107`, `466`–`468`, `DRIFT-043`.
- Related skills: `register-identifiers`, `column-tooltip`.
