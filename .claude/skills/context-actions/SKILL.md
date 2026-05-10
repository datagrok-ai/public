---
name: context-actions
description: Register a context-specific action in a Datagrok package so the platform offers it on right-click and in the Actions accordion when the user clicks an item of the matching semantic type
---

# context-actions

## When to use

Your package owns a domain-specific item (a molecule cell, a sequence
column, a multi-column selection, …) and you want the platform to offer
a custom verb on it: "Use as filter" on a `Molecule` cell, "To Atomic
Level..." on a `Macromolecule` column header, "Sparklines..." on
selected numerical columns. Triggers: "add a right-click action for
`<SemType>`", "add an entry to the Actions pane for `<SemType>`".

## Prerequisites

- A package scaffold (`grok create <Name>`); commands run from the package root.
- `datagrok-api` available — the article snippet omits imports and won't
  compile as written.
- The dispatch key already exists: a semantic type registered for the
  target columns (built-in like `Molecule` / `Macromolecule`, or
  package-registered via `package.json` `meta.semanticTypes[]` /
  `DG.SemanticValue.registerRegExpDetector` — see `register-identifiers`).

## Steps

1. **Pick the input type — it is the invocation surface.**
   Four type tokens are dispatched (`DG-FACT-182`, `DG-FACT-DRIFT-075`);
   the article documents only `string`, which loses cell context.
   - `semantic_value <name> { semType: <T> }` → bound to `DG.SemanticValue`;
     fires on **cell** right-click and per-row Actions pane. Production
     canonical for cell actions (`packages/Chem/src/package.g.ts:843-849`,
     `packages/Bio/src/package.g.ts (Edit Helm…)`). Use this by default.
   - `column <name> { semType: <T> }` → bound to `DG.Column`; fires on
     **column-header** right-click (`packages/Bio/src/package.g.ts:391-396`).
   - `list <name> { type: numerical }` → bound to `DG.Column[]`; fires on
     **multi-column selection** (`packages/PowerGrid/src/package.g.ts:171-176`).
   - `string <name> { semType: <T> }` → article form; works but gives
     only the raw value, no cell/column reference (`DG-FACT-DRIFT-072`).

2. **Declare the action with `meta.action: <Display Text>`.**
   The presence of `meta.action` in the function's annotations is the
   registration — no explicit `register(...)` call needed
   (`DG-FACT-180`). The text is the menu label. Exactly one input,
   carrying a `semType` annotation, is required (`DG-FACT-181`); without
   `semType` there is no dispatch surface and the action is silently
   never offered.

   ```typescript
   // src/package.ts — canonical decorator form (Chem useAsSubstructureFilter)
   import * as grok from 'datagrok-api/grok';
   import * as DG   from 'datagrok-api/dg';
   export const _package = new DG.Package();

   export class PackageFunctions {
     @grok.decorators.func({
       name: 'Use as filter',
       description: 'Adds this structure as a substructure filter',
       meta: {action: 'Use as filter'},
     })
     static useAsSubstructureFilter(
       @grok.decorators.param({options: {semType: 'Molecule'}})
       value: DG.SemanticValue
     ): void {
       const tv = grok.shell.tv;
       if (tv == null) throw new Error('Requires an open table view.');
       const molCol = value.cell?.column;
       if (molCol == null) throw new Error('Molecule column not found.');
       // … build molblock from value.value, branching on
       // value.cell.column.meta.units (see step 4)
       tv.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
         type: DG.FILTER_TYPE.SUBSTRUCTURE,
         column: molCol.name,
         columnName: molCol.name,
         molBlock: /* … */ '',
       }, false);
     }
   }
   ```
   Expected: `src/package.g.ts` (auto-generated) gains a wrapper with
   `//name: Use as filter`, `//description: …`,
   `//input: semantic_value value { semType: Molecule }`, and
   `//meta.action: Use as filter` (compare
   `packages/Chem/src/package.g.ts:843-849`). The header-comment form
   (article `//meta.action: …`) is equivalent (`DG-FACT-184`).

3. **Scope WHERE the action surfaces with the two `meta.exclude-*` toggles.**
   Both surfaces (cell right-click menu and Actions accordion pane) light
   up by default. Suppress one with a sibling `meta.*: true` flag — the
   article does not mention either (`DG-FACT-183`, `DG-FACT-DRIFT-074`):
   - `meta.exclude-actions-panel: true` → keeps right-click, hides from
     Actions pane. Chem applies this to every `Copy as <format>` sub-action
     so the panel only shows the parent `Copy as...` chooser
     (`packages/Chem/src/package.ts:2066,2079`,
     `packages/Chem/src/package.g.ts:864,873,882`).
   - `meta.exclude-current-value-menu: true` → keeps Actions pane, hides
     from right-click. Chem uses it on the parent `Copy as...` chooser
     (`packages/Chem/src/package.ts:2045`).
   ```typescript
   @grok.decorators.func({
     name: 'Copy as SMILES',
     meta: {action: 'Copy as SMILES', 'exclude-actions-panel': 'true'},
   })
   static copyAsSmiles(
     @grok.decorators.param({options: {semType: 'Molecule'}})
     value: DG.SemanticValue): void { /* … */ }
   ```

4. **For "Use as filter" on `Molecule` — branch on units, then `updateOrAdd`.**
   The article uses `.add(...)` and a bare `FILTER_TYPE.SUBSTRUCTURE`.
   Production uses `.updateOrAdd(..., false)` (so a second invocation
   replaces the existing filter row instead of stacking duplicates) and
   the fully-qualified `DG.FILTER_TYPE.SUBSTRUCTURE` (`DG-FACT-185`,
   `DG-FACT-DRIFT-073`). When the cell's `meta.units == DG.chem.Notation.Smiles`,
   convert via `DG.chem.convertMolNotation` to preserve orientation;
   otherwise call `molToMolblock(molecule, getRdKitModule())`.

5. **Build and publish — registration is automatic.**
   ```bash
   npm install
   grok check                     # exits 0
   grok publish <host>            # add --release once stable
   ```
   No explicit `register(...)` — the `meta.action` annotation IS the
   registration. After deploy, right-click a cell / column whose
   `semType` matches the input constraint to see the action; expand the
   Actions pane on that item to confirm panel placement.

## Common failure modes

- **Action never appears on the right-click menu or Actions pane.** The
  `meta.action` token is missing or the input lacks a `semType`.
  Inspect the auto-generated `src/package.g.ts`: there MUST be both a
  `//meta.action: <Text>` line and an `//input: ...{ semType: <T> }`
  annotation (`DG-FACT-180`, `DG-FACT-181`). Without `semType` the
  function registers but has no dispatch surface.
- **Action fires on the wrong column when the table has two `Molecule`
  columns.** The body uses `string` input + `tv.dataFrame.columns.bySemType(...)`
  to look up the column — that picks the first match, not the column
  the user clicked (`DG-FACT-DRIFT-072`). Switch the input to
  `semantic_value` and read `value.cell?.column` directly.
- **Action surfaces on the wrong invocation surface** (e.g. appears on
  the column header but not on cells). Wrong input type: `column` fires
  on the header, `semantic_value` on the cell, `list` on multi-column
  selection (`DG-FACT-182`).
- **Article snippet won't compile / `FILTER_TYPE` is undefined.** The
  article omits imports and uses bare `FILTER_TYPE.SUBSTRUCTURE`. Add
  `import * as DG from 'datagrok-api/dg';` and use
  `DG.FILTER_TYPE.SUBSTRUCTURE` (`DG-FACT-DRIFT-073`).
- **Each invocation of "Use as filter" stacks a new filter row.** The
  article uses `.add(...)`; switch to `.updateOrAdd(filterDescriptor, false)`
  so the existing row for that column is replaced
  (`DG-FACT-185`, `DG-FACT-DRIFT-073`).
- **`Copy as <format>` sub-actions clutter the Actions pane.** Each
  format variant lacks `meta.exclude-actions-panel: true`. Add the flag
  to every variant; keep `meta.exclude-current-value-menu: true` on the
  parent `Copy as...` chooser (`DG-FACT-183`).

## Verification

- `grok check` exits `0`; `grok publish <host>` exits `0`.
- The regenerated `src/package.g.ts` contains, for each registered
  action, a wrapper carrying both `//meta.action: <Text>` and
  `//input: <type> <name> { semType: <T> }` (compare
  `packages/Chem/src/package.g.ts:843-849`).
- In Datagrok, open a table with a column whose `semType` matches the
  input constraint. Right-click a cell (for `semantic_value` input), the
  column header (for `column`), or a multi-column selection (for `list`)
  — the action label appears with the text from `meta.action`.
- Expand the Actions accordion in the property panel for the same item
  — the action appears unless `meta.exclude-actions-panel: true` was set.
- Right-click an item with a different `semType` — the action does NOT
  appear, confirming the dispatch is scoped, not global.

## See also

- Source articles:
  - `help/develop/how-to/ui/context-actions.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-180` … `DG-FACT-185` and drifts
  `DG-FACT-DRIFT-072` … `DG-FACT-DRIFT-075`.
- Reference packages:
  - `packages/Chem/src/package.ts:2010-2087` — decorator-form
    `useAsSubstructureFilter` (`semType: Molecule`,
    `meta: {action: 'Use as filter'}`) and the `Copy as...` family
    showing `exclude-actions-panel` / `exclude-current-value-menu`.
  - `packages/Chem/src/package.g.ts:843-905` — auto-emitted header-form
    wrappers (the surface the platform actually reads).
  - `packages/Bio/src/package.g.ts:391-396` — `column`-input action
    (`To Atomic Level...`, `semType: Macromolecule`) firing on the
    column-header right-click.
  - `packages/PowerGrid/src/package.g.ts:171-183` — `list`-input
    actions (`Sparklines...`, `Smart form...`) firing on multi-column
    selection.
- Related skills:
  - `register-identifiers` (sibling — registers the `semType` that a
    context action dispatches on).
