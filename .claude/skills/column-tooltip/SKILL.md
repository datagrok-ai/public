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
   a `semType`.** Use generic `@grok.decorators.func({meta: {role:
   'tooltip'}})` — there's no specialized `@grok.decorators.tooltip()`
   (see `DG-FACT-095` role string, `DG-FACT-096` semType binding,
   `DG-FACT-097` decorator form).

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

   Three rules govern the body:

   - Build from `col` / `col.dataFrame`, never `grok.shell.tv.dataFrame`
     or `grok.shell.t` — the hovered column may belong to a different
     table (`DG-FACT-427`).
   - Return `undefined` to opt out per-column (e.g. SMARTS-pattern
     molecules); platform falls back to default summary (`DG-FACT-098`).
   - Wrap result as `DG.Widget` via `DG.Widget.fromRoot(root)` or return
     a `DG.Widget` subclass (`DG-FACT-098`).

2. **Build, and verify the auto-emitted wrapper in `package.g.ts`.**
   `npm run build` runs `FuncGeneratorPlugin`, which regenerates
   `src/package.g.ts` from decorated methods — never hand-edit it
   (`DG-FACT-428`).
   ```bash
   npm install && npm run build
   ```
   Expected: the emitted wrapper carries `//meta.role: tooltip` and
   `//input: column col { semType: <Type> }` (`DG-FACT-095`,
   `DG-FACT-096`).

3. **Publish.** The `tooltip` role IS the registration; no explicit
   `grok.functions.register(...)` call is needed (`DG-FACT-428`).
   ```bash
   grok publish <host>   # add --release once stable
   ```

## Common failure modes

- Hover never fires — check `package.g.ts` for `//meta.role: tooltip` (lowercase) and `//input: column col { semType: <Type> }` (`DG-FACT-095`, `DG-FACT-096`).
- `semType` mismatch — `col.semType` is exact-match, case-sensitive (`DG-FACT-096`).
- Wrong-table data with multiple TableViews — replace `grok.shell.tv.dataFrame`/`grok.shell.t` with `col`/`col.dataFrame` (`DG-FACT-427`).
- Article snippet won't compile — add the three imports shown above.
- Returned `DG.Viewer` or `HTMLElement` — wrap with `DG.Widget.fromRoot(...)` (`DG-FACT-098`).

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
