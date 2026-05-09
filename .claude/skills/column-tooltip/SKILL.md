---
name: column-tooltip
description: Register a custom column tooltip in a Datagrok package so the platform shows a domain-specific widget when the user hovers a column header
---

# column-tooltip

## When to use

Your package owns a domain-specific column type (sequences, molecules,
plates, identifiers, …) and you want the platform to render a richer
tooltip than the default category/min/max summary when the user hovers
that column's header. Triggers: "show a logo plot when hovering a
`Macromolecule` column header", "show top-N diverse structures for
`Molecule`", "give my `SampleId` columns a custom tooltip widget".

## Prerequisites

- A package scaffold (`grok create <Name>`); commands run from the package root.
- `datagrok-api` available (`import * as grok / DG`); the article omits
  imports and won't compile as written.
- A semantic type already attached to the column (either a built-in
  like `Macromolecule`/`Molecule`, or one your package registers — see
  `register-identifiers`). The tooltip binds to columns by `semType`.
- Familiarity with `DG.Widget` (`DG.Widget.fromRoot(...)` to wrap an
  `HTMLElement`, or a subclass producing the widget root).

## Steps

1. **Pick a registration form (decorator vs. header comments).**
   The platform discovers tooltips via the function role `tooltip`
   (camelCase, `DG-FACT-095`). There is NO specialized
   `@grok.decorators.tooltip()` (unlike `fileHandler` / `fileViewer`)
   — production packages use the generic `func` decorator with
   `meta.role: 'tooltip'` (`DG-FACT-097`). Two equivalent surfaces:
   - **Decorator (canonical, used by Bio and Chem):** a `static`
     method on `PackageFunctions` decorated with
     `@grok.decorators.func({meta: {role: 'tooltip'}})`, with the
     column input typed via
     `@grok.decorators.param({options: {semType: '<Type>'}}) col: DG.Column`.
   - **Header comments (article form):** an exported function annotated
     with `//meta.role: tooltip`, `//input: column col {semType: <Type>}`,
     `//output: widget result`. This form survives in auto-emitted
     `package.g.ts` (`DG-FACT-DRIFT-039`).

2. **Bind the tooltip to a semantic type via the `col` param.**
   The `semType` constraint on the `column` input is what scopes the
   tooltip to columns of that type — the platform calls the function
   only for columns whose `semType` matches (`DG-FACT-096`). Use the
   exact semantic-type string registered in `package.json`
   `meta.semanticTypes[]` or by `DG.SemanticValue.registerRegExpDetector`.

   ```typescript
   // src/package.ts — canonical decorator form (Bio sequenceTooltip)
   import * as grok from 'datagrok-api/grok';
   import * as DG   from 'datagrok-api/dg';
   import * as ui   from 'datagrok-api/ui';
   export const _package = new DG.Package();

   export class PackageFunctions {
     @grok.decorators.func({
       meta: {role: 'tooltip'},
       tags: ['tooltip'],          // optional; Bio sets it, Chem omits it
     })
     static sequenceTooltip(
       @grok.decorators.param({options: {semType: 'Macromolecule'}})
       col: DG.Column
     ): DG.Widget<any> {
       const root = ui.div([ui.divText(`${col.name}: ${col.length} rows`)]);
       return DG.Widget.fromRoot(root);
     }
   }
   ```
   Expected: build succeeds; `src/package.g.ts` gains a wrapper with
   `//input: column col { semType: Macromolecule }`,
   `//output: widget result`, `//meta.role: tooltip`
   (compare `packages/Bio/src/package.g.ts:16-22`,
   `packages/Chem/src/package.g.ts:38-44`).

3. **Build the widget from `col` — never from `grok.shell.tv.dataFrame`.**
   The article's example pulls the dataframe from the *active TableView*
   (`grok.shell.tv.dataFrame.plot.fromType(...)`), but the column being
   tooltipped may belong to a different dataframe than `grok.shell.tv` —
   silently rendering the wrong tooltip (`DG-FACT-DRIFT-040`).
   Production packages take `col` (or `col.dataFrame`) directly:
   - Bio: `new MacromoleculeColumnWidget(col, _package.seqHelper)`
     (`packages/Bio/src/package.ts:152`).
   - Chem: iterates `col.categories` directly
     (`packages/Chem/src/package.ts:273-277`).

   ```typescript
   // good — bound to the column being hovered
   static sequenceTooltip(
     @grok.decorators.param({options: {semType: 'Macromolecule'}})
     col: DG.Column
   ): Promise<DG.Widget<any>> {
     const viewer = await col.dataFrame.plot.fromType(
       'WebLogo', {sequenceColumnName: col.name});
     return DG.Widget.fromRoot(viewer.root);
   }
   ```

4. **Return `undefined` to opt out per-column.**
   The signature is `tooltip(col: Column): Widget` (`DG-FACT-095`), but
   returning `undefined` (or omitting the return) signals "no tooltip
   applies here" — the platform falls back to the default tooltip
   (`DG-FACT-098`). Chem uses this to skip SMARTS-only `Molecule`
   columns (`packages/Chem/src/package.ts:267-277`):
   ```typescript
   static async chemTooltip(
     @grok.decorators.param({options: {semType: 'Molecule'}})
     col: DG.Column
   ): Promise<DG.Widget | undefined> {
     for (let i = 0; i < Math.min(col.categories.length, 100); ++i)
       if (col.categories[i] && _isSmarts(col.categories[i])) return; // → default tooltip
     // … build and return widget
   }
   ```

5. **Build, publish, and let the platform auto-register.**
   ```bash
   npm install
   grok check                     # exits 0
   grok publish <host>            # add --release once stable
   ```
   No explicit `register(...)` call is needed — the function role
   `tooltip` IS the registration. After deploy, hover any column
   whose `semType` matches the constraint to see the widget
   (`DG-FACT-095`, `DG-FACT-096`).

## Common failure modes

- **Tooltip never fires; the platform shows the default tooltip.**
  The role token is wrong or missing. Inspect `src/package.g.ts`: it
  MUST contain `//meta.role: tooltip` (camelCase) and a
  `//input: column col { semType: <Type> }` line (`DG-FACT-095`,
  `DG-FACT-096`). The token is case-sensitive — `Tooltip` / `tooltips`
  won't register.
- **Tooltip fires for the wrong columns (or never).** The `semType`
  string on the `col` param doesn't match what the column actually has.
  Verify in the UI: right-click the column → *Properties* → check
  `semType`; or in code, `col.semType`. The constraint is exact-match,
  case-sensitive (`DG-FACT-096`).
- **Tooltip renders wrong data when multiple TableViews are open.**
  The body builds the widget from `grok.shell.tv.dataFrame` instead of
  `col` / `col.dataFrame`. Switch to the column reference passed in
  (`DG-FACT-DRIFT-040`).
- **Article snippet copy-pasted, won't compile.** The article shows a
  bare `export async function` with header comments and no imports
  (`DG-FACT-DRIFT-039`). Translate to the decorator form on
  `PackageFunctions` with explicit
  `import * as grok from 'datagrok-api/grok';` and
  `import * as DG from 'datagrok-api/dg';`, and either `return` a
  `DG.Widget` (sync) or `Promise<DG.Widget | undefined>` (async).
- **Function returned a viewer / `HTMLElement` instead of a `Widget`.**
  The signature is `tooltip(col: Column): Widget` (`DG-FACT-095`,
  `DG-FACT-098`). Wrap an element with `DG.Widget.fromRoot(root)`;
  for a `DG.Viewer`, return `DG.Widget.fromRoot(viewer.root)`.

## Verification

- `grok check` exits `0`; `grok publish <host>` exits `0`.
- The regenerated `src/package.g.ts` contains, for each tooltip
  function, a wrapper with `//meta.role: tooltip`,
  `//input: column col { semType: <Type> }`, and
  `//output: widget result` (compare
  `packages/Bio/src/package.g.ts:16-22`,
  `packages/Chem/src/package.g.ts:38-44`).
- In Datagrok, open a table with a column whose `semType` matches
  the constraint and hover the column header: the platform invokes
  your function and renders the returned `DG.Widget` in the tooltip
  popup.
- Hover a column with a *different* `semType`: the platform falls
  back to the default tooltip — confirms the binding is scoped, not
  global.

## See also

- Source articles:
  - `help/develop/how-to/grid/column-tooltip.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-095` … `DG-FACT-098` and drifts
  `DG-FACT-DRIFT-039`, `DG-FACT-DRIFT-040`.
- Reference packages:
  - `packages/Bio/src/package.ts:146-159` — decorator-form
    `sequenceTooltip` for `semType: Macromolecule`, returns
    `DG.Widget<any>` (sync) backed by `MacromoleculeColumnWidget`.
  - `packages/Chem/src/package.ts:263-277` — decorator-form
    `chemTooltip` for `semType: Molecule`, returns
    `Promise<DG.Widget | undefined>` and opts out on SMARTS columns.
  - `packages/Bio/src/package.g.ts:16-22`,
    `packages/Chem/src/package.g.ts:38-44` — auto-emitted header-form
    wrappers (the surface the platform actually reads).
- Related skills:
  - `register-identifiers` (sibling — registers the `semType` that a
    column tooltip binds to).
