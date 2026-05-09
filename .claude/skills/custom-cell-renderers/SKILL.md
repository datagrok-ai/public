---
name: custom-cell-renderers
description: Register a custom Datagrok grid cell renderer so the platform paints a domain-specific visualization for matching columns
---

# custom-cell-renderers

## When to use

Your package needs the grid to draw something other than the default
text/number cell — stars for a rating, a colored swatch for a hex code,
an SVG molecule, a pie/bar/radar sparkline summarising several columns,
a custom hyperlink badge — and you want the platform to bind the
renderer automatically whenever a matching column appears.

## Prerequisites

- A package scaffold (`grok create <Name>`); run from the package root.
- `datagrok-api` imports (`* as DG`, `* as grok` — article omits these).
- `datagrok-tools` shipping `@grok.decorators.cellRenderer` (≥ 4.12.x;
  current line is 6.x — `DG-FACT-DRIFT-041`).
- Familiarity with `CanvasRenderingContext2D` — `render` paints directly
  to a 2D canvas; nothing in the DOM.

## Steps

1. **Pick a registration form.**
   The platform discovers renderers via the function role `cellRenderer`
   (camelCase — `DG-FACT-099`). Two surfaces compile to the same
   `package.g.ts` wrapper:
   - **Class decorator (canonical for plain renderers).**
     `@grok.decorators.cellRenderer({...})` on the class — no factory
     in `package.ts` (`DG-FACT-103`). Used by `StarsCellRenderer`,
     `SvgCellRenderer`, `HyperlinkCellRenderer`, `ColorCellRenderer`.
   - **Function decorator (canonical for virtual sparklines).**
     Static factory on `PackageFunctions`,
     `@grok.decorators.func({meta: {role:'cellRenderer', cellType:'…',
     virtual:'true'}, tags:['cellRenderer'],
     outputs:[{type:'grid_cell_renderer', name:'result'}]})`
     (`DG-FACT-104`). Used by `piechartCellRenderer`, `radarCellRenderer`.

   Do NOT paste the article's bare-function `//name:` / `//meta.role:`
   block into `package.ts` — that is the auto-emitted `package.g.ts`
   shape, not the authoring surface (`DG-FACT-DRIFT-042`).

2. **Subclass `DG.GridCellRenderer`; override the three required hooks.**
   `name`, `cellType`, and `render` all throw in the base class
   (`DG-FACT-102`). Optional overrides — `defaultWidth/Height`,
   `renderSettings(gridColumn)`, mouse/key handlers (`onClick`,
   `onMouseEnter`, …), `hasContextValue` / `getContextValue`, `clip`.

   ```typescript
   // src/cell-types/stars-cell-renderer.ts — class-decorator form
   import * as DG from 'datagrok-api/dg';
   import * as grok from 'datagrok-api/grok';

   @grok.decorators.cellRenderer({
     name: 'Stars',
     cellType: 'Stars',
     // @ts-ignore — `tags` is on FuncOptions, not CellRendererOptions
     tags: ['cellRenderer'],
   })
   export class StarsCellRenderer extends DG.GridCellRenderer {
     get name()     { return 'Stars'; }
     get cellType() { return 'Stars'; }

     render(
       g: CanvasRenderingContext2D,
       x: number, y: number, w: number, h: number,
       gridCell: DG.GridCell, cellStyle: DG.GridCellStyle,
     ) {
       const value = Math.round(gridCell.cell.value);
       g.fillStyle = '#FFB400';
       g.fillText('★'.repeat(value), x + 4, y + h / 2 + 4);
     }
   }
   ```
   Expected: `src/package.g.ts` (regenerated on `npm run build`) gains a
   wrapper with `//tags: cellRenderer`, `//output: grid_cell_renderer
   renderer`, `//meta.role: cellRenderer`, `//meta.cellType: Stars`
   (compare `packages/PowerGrid/src/package.g.ts:60-65`).

3. **Bind the renderer to columns via `cellType` (or `columnTags`).**
   The platform invokes your renderer when a column's effective cell
   type matches `meta.cellType`. Two binding paths:
   - **By semantic type.** Set `column.semType = 'Stars'` (e.g. via a
     `semanticTypeDetector`); the platform looks up the renderer by
     `cellType` (`DG-FACT-099`).
   - **By column tags.** Use `meta.columnTags: 'units=kg,foo=bar'` to
     match on arbitrary tag pairs instead of `semType` — see
     `packages/PowerGrid/src/package.g.ts:185-201` (`DG-FACT-106`).

   For SUMMARY/VIRTUAL renderers (sparklines, pie chart, smart form —
   columns synthesised at render time from several numeric columns),
   add `meta.virtual` so the column appears in the *Add column → New
   chart* menu (`DG-FACT-105`):

   ```typescript
   // src/package.ts — function-decorator form for a virtual sparkline
   import * as grok from 'datagrok-api/grok';
   import {PieChartCellRenderer} from './sparklines/piechart';
   export const _package = new DG.Package();
   export class PackageFunctions {
     @grok.decorators.func({
       meta: {role: 'cellRenderer', cellType: 'piechart', virtual: 'true'},
       tags: ['cellRenderer'],
       outputs: [{type: 'grid_cell_renderer', name: 'result'}],
       name: 'Pie Chart',
     })
     static piechartCellRenderer() { return new PieChartCellRenderer(); }
   }
   ```
   Note: `virtual: true` (BOOLEAN) on the class decorator,
   `virtual: 'true'` (STRING) on the function decorator — the `meta`
   field is typed `Record<string, string>` in the function form
   (`DG-FACT-DRIFT-044`).

4. **Build, publish, and let the codegen wire it up.**
   ```bash
   npm install
   npm run build       # FuncGeneratorPlugin emits/updates src/package.g.ts
   grok check          # exits 0
   grok publish <host> # add --release once stable
   ```
   `package.g.ts` is meant to be committed (`DG-FACT-107`). No explicit
   `DG.GridCellRenderer.register(...)` call — the function role IS the
   registration.

## Common failure modes

- **Renderer never fires; cells fall back to default text.** The role
  string is wrong or missing. Inspect `src/package.g.ts`: each entry
  MUST contain `//meta.role: cellRenderer` (camelCase) and
  `//meta.cellType: <Type>`. The token is case-sensitive — `CellRenderer`
  / `cell-renderer` (the lookalike COLUMN-TAG key at
  `js-api/src/const.ts:330`) won't register (`DG-FACT-099`).
- **Build emits no entry in `package.g.ts`.** `FuncGeneratorPlugin` is
  not wired into `webpack.config.js` (`DG-FACT-107`). Add
  `const FuncGeneratorPlugin = require('datagrok-tools/plugins/func-gen-plugin');`
  and `new FuncGeneratorPlugin({outputPath: './src/package.g.ts'})`.
- **`'cellType'`/`'name'`/`'Not implemented'` thrown at render time.**
  The subclass forgot to override one of `get name()`, `get cellType()`,
  or `render(...)` — all three throw in `DG.GridCellRenderer` defaults
  (`DG-FACT-102`).
- **TypeScript: `'tags' does not exist on type CellRendererOptions`.**
  The option type omits `tags`, but every canonical PowerGrid renderer
  still emits `tags: ['cellRenderer']` so the role descriptor's
  `header: 'tags'` bucket fills (`DG-FACT-DRIFT-043`). Suppress with
  `// @ts-ignore` on the line above the `tags` entry (matches
  `packages/PowerGrid/src/cell-types/stars-cell-renderer.ts:40`).
- **Sparkline never appears in *Add column → New chart*.** Renderer
  is registered without `meta.virtual`; non-virtual renderers fire
  only on existing columns (`DG-FACT-105`). Add `virtual: true` (class)
  or `virtual: 'true'` (function form).
- **Article snippet pasted verbatim, won't compile.** The article shows
  the bare-function `//name:` / `//meta.role:` header block — that is
  the auto-emitted `package.g.ts` shape, not what you author today
  (`DG-FACT-DRIFT-042`). Translate to one of the decorator forms above.

## Verification

- `npm run build` exits `0`; `grok check` exits `0`; `grok publish`
  exits `0`.
- Regenerated `src/package.g.ts` contains, for each renderer, a wrapper
  with `//tags: cellRenderer`, `//output: grid_cell_renderer <name>`,
  `//meta.role: cellRenderer`, `//meta.cellType: <Type>` (compare
  `packages/PowerGrid/src/package.g.ts:13-65` for plain renderers,
  `116-125` for virtual sparklines).
- In Datagrok, open a table with a column whose `semType`/tags match
  `cellType` — the grid paints with your `render`. For a virtual
  renderer, click *Add column → New chart* and confirm `<Type>` appears
  in the menu.
- Hover a column with a different `semType`: the platform falls back
  to the built-in renderer — confirms the binding is scoped, not
  global.

## See also

- Source articles:
  - `help/develop/how-to/grid/custom-cell-renderers.md`
  - mirror: `docs/_internal/articles-mirror/how-to/grid/custom-cell-renderers.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-099` … `DG-FACT-107` and drifts `DG-FACT-DRIFT-041` …
  `DG-FACT-DRIFT-044`.
- Reference packages:
  - `packages/PowerGrid/src/cell-types/stars-cell-renderer.ts:37-43` —
    minimal class-decorator form (`name`/`cellType`/`render`/`onClick`).
  - `packages/PowerGrid/src/cell-types/svg-cell-renderer.ts:6-12` —
    class-decorator form with image cache + `onDoubleClick`.
  - `packages/PowerGrid/src/sparklines/piechart.ts:192-326` +
    `packages/PowerGrid/src/package.ts:85-98` — function-decorator
    form with `defaultWidth/Height`, `renderSettings`, `onMouseMove`,
    `getContextValue`, virtual.
  - `packages/PowerGrid/src/package.g.ts:13-201` — auto-emitted
    header-form wrappers (what the platform reads).
- Related skills:
  - `register-identifiers` (sibling — registers the `semType` that
    `cellType` is matched against).
  - `column-tooltip` (sibling — same role-based registration pattern,
    different role token: `tooltip` vs `cellRenderer`).
