---
name: column-tooltip
version: 0.4.0
description: |
  Register a function in a Datagrok package that the platform invokes
  whenever the user hovers a column header carrying a given semantic
  type, returning a `DG.Widget` rendered in place of the default
  category/min/max summary. For plugin authors whose domain columns
  (sequences, molecules, identifiers, plate IDs) deserve a richer
  preview than the platform's built-in popup. Produces a `static`
  method on `PackageFunctions` decorated with
  `@grok.decorators.func({meta: {role: 'tooltip'}})`, plus the
  auto-emitted `package.g.ts` wrapper the platform reads at startup.
  Use when asked to "show a logo when hovering a sequence column",
  "render a molecule depiction in the column-header hover",
  or "replace the default header-hover summary for our domain columns".
triggers:
  - show a richer column header hover
  - customize the header-hover preview
  - render a domain widget on column hover
  - replace the default column header popup
  - per-semantic-type hover widget
  - molecule depiction on column hover
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# column-tooltip

## When to use

Your package owns a domain-specific column type (sequences, molecules,
plate barcodes, internal compound IDs) and you want the platform to
render a richer widget when the user hovers that column's header — not
the default category/min/max summary. Real-world phrasings: "show a
WebLogo for `Macromolecule` columns", "top-N diverse structures for
`Molecule`", "custom hover preview for `SampleId` columns".

## Prerequisites

- A semantic type already attached to the column — built-in
  (`Macromolecule`, `Molecule`, …) or one your package registers via
  the `register-identifiers` skill (`DG-FACT-096`).
- Familiarity with `DG.Widget` — wrap an `HTMLElement` via
  `DG.Widget.fromRoot(root)` or return a subclass (`DG-FACT-098`).

## Steps

1. **Declare a `tooltip`-role function on `PackageFunctions`, bound to
   a `semType`.** The platform discovers tooltips by the function role
   `tooltip` (lowercase, `DG-FACT-095`); the `semType` constraint on
   the `column` input scopes it to that type (exact-match, case-sensitive,
   `DG-FACT-096`). There is no specialized `@grok.decorators.tooltip()`
   — use the generic `@grok.decorators.func` with `meta.role: 'tooltip'`
   (`DG-FACT-097`).

   ```typescript
   // src/package.ts
   import * as grok from 'datagrok-api/grok';
   import * as DG from 'datagrok-api/dg';
   import * as ui from 'datagrok-api/ui';

   export const _package = new DG.Package();

   export class PackageFunctions {
     @grok.decorators.func({meta: {role: 'tooltip'}})
     static async sequenceTooltip(
       @grok.decorators.param({options: {semType: 'Macromolecule'}})
       col: DG.Column,
     ): Promise<DG.Widget<any> | undefined> {
       const viewer = await col.dataFrame.plot.fromType(
         'WebLogo', {sequenceColumnName: col.name});
       return DG.Widget.fromRoot(viewer.root);
     }
   }
   ```

   Three rules govern the body (`semType` only scopes the hook):

   - **Build from `col`, never from `grok.shell.tv.dataFrame`.** The
     hovered column may belong to a different table than the active
     `TableView`, so `grok.shell.tv.dataFrame` / `grok.shell.t`
     silently resolves to the wrong dataframe and any `columnName`
     option (`sequenceColumnName: col.name`) misses. The article makes
     this explicit (`column-tooltip.md:13-16, :24`, `DG-FACT-427`); Bio
     follows it at `packages/Bio/src/package.ts:146-159`.

   - **Return `undefined` to opt out per-column.** Signature is
     `tooltip(col: Column): Widget | undefined` (`DG-FACT-095`); a bare
     `return;` falls back to the default summary (`DG-FACT-098`). Chem
     uses this to skip SMARTS-pattern `Molecule` columns — see
     `packages/Chem/src/package.ts:266-280` for the early-return loop
     over `col.categories`.

   - **Wrap the result as a `DG.Widget`.** `DG.Widget.fromRoot(root)`
     for an `HTMLElement`, `DG.Widget.fromRoot(viewer.root)` for a
     viewer, or return a `DG.Widget` subclass directly (`DG-FACT-098`).

   Expected: the file type-checks (the decorators resolve, `DG.Column`
   and `DG.Widget` are in scope from the imports above).

2. **Build, and verify the auto-emitted wrapper in `package.g.ts`.**
   `npm run build` runs `FuncGeneratorPlugin`, which scans
   `PackageFunctions` for decorated methods and regenerates
   `src/package.g.ts` with one wrapper per function (`DG-FACT-428`). Do
   NOT hand-edit `package.g.ts`; it is overwritten on every build.
   ```bash
   npm install && npm run build
   ```
   Expected: `npm run build` exits 0; the regenerated
   `src/package.g.ts` gains a wrapper of the shape

   ```
   //input: column col { semType: Macromolecule }
   //output: widget result
   //meta.role: tooltip
   export async function sequenceTooltip(col: DG.Column) {
     return await PackageFunctions.sequenceTooltip(col);
   }
   ```

   The MUST-have lines are `//meta.role: tooltip` (lowercase) and
   `//input: column col { semType: <Type> }` (`DG-FACT-095`,
   `DG-FACT-096`). Chem's emitted wrapper has exactly this shape —
   compare `packages/Chem/src/package.g.ts:38-44`. A `//tags: tooltip`
   line appears ONLY if the decorator also passes `tags: ['tooltip']`
   explicitly — Bio does this, so `packages/Bio/src/package.g.ts:16-22`
   carries the extra line; the example above does not, so don't look
   for it.

3. **Publish.** The `tooltip` role IS the registration — the platform
   discovers the function from the `package.g.ts` wrapper on load
   (`DG-FACT-428`); no explicit `grok.functions.register(...)` call is
   needed.
   ```bash
   grok publish <host>   # add --release once stable
   ```
   Expected: `grok publish <host>` exits 0 and reports the package
   version it pushed.

## Common failure modes

- **Hover never fires; default summary still shows.** The role token
  is missing or misspelled. Inspect `src/package.g.ts` — it MUST
  contain `//meta.role: tooltip` (lowercase) and a
  `//input: column col { semType: <Type> }` line (`DG-FACT-095`,
  `DG-FACT-096`). `Tooltip`/`tooltips` won't match.
- **Hover fires on the wrong columns (or never).** The `semType` string
  on the `col` param doesn't match what the column actually has.
  In the UI, right-click the column → *Properties* → check `semType`;
  in code, `col.semType`. Exact-match and case-sensitive (`DG-FACT-096`).
- **Hover renders the wrong table's data when multiple `TableViews`
  are open.** The body uses `grok.shell.tv.dataFrame` (or `grok.shell.t`)
  instead of `col` / `col.dataFrame`. Switch to the column reference
  passed in (`DG-FACT-427`).
- **Article snippet copy-pasted, won't compile.** The article omits
  imports (`column-tooltip.md:18-27` shows only the class body). Add
  `import * as grok from 'datagrok-api/grok'`,
  `import * as DG from 'datagrok-api/dg'`,
  `import * as ui from 'datagrok-api/ui'` to the top of `src/package.ts`.
- **Returned a `DG.Viewer` or `HTMLElement` instead of a `DG.Widget`.**
  Wrap with `DG.Widget.fromRoot(root)`; for a viewer, return
  `DG.Widget.fromRoot(viewer.root)` (`DG-FACT-098`).

## See also

- Source article: `help/develop/how-to/grid/column-tooltip.md`.
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-095/096/097/098/427/428`.
- Reference packages:
  - `packages/Bio/src/package.ts:146-159` + `Bio/src/package.g.ts:16-22`
    — `sequenceTooltip` for `Macromolecule`.
  - `packages/Chem/src/package.ts:266-280` + `Chem/src/package.g.ts:38-44`
    — `chemTooltip` for `Molecule`, opts out on SMARTS.
- Related skills: `register-identifiers` (registers the `semType` this
  skill binds to); `custom-cell-renderers` (same `col` / `grid.dataFrame`
  discipline for per-cell rendering).
