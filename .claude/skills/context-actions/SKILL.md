---
name: context-actions
version: 0.1.0
description: |
  Register a Datagrok package function so the platform offers it as a
  custom verb when the user right-clicks (or opens the Actions
  accordion on) an item whose semantic type matches the declared input.
  For plugin authors who want a domain-specific entry — "Use as filter"
  on a `Molecule` cell, "To Atomic Level..." on a `Macromolecule` column
  header, "Sparklines..." on a multi-column selection — without writing
  any explicit `register(...)` call. Produces an annotated
  `PackageFunctions` method plus the auto-emitted `package.g.ts`
  wrapper the platform reads at startup.
  Use when asked to "add a right-click action for compound cells", "put
  an entry in the Actions pane for a sequence column", or "extend the
  right-click menu conditionally on a column type".
triggers:
  - extend the right-click menu
  - add a menu item for a column type
  - entry in the Actions accordion
  - custom verb on a cell
  - hook a context-menu action to a semantic type
  - column-header right-click handler
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# context-actions

## When to use

Your package owns a domain-specific item (a molecule cell, a sequence
column, a multi-column selection, …) and wants the platform to offer a
custom verb on it: "Use as filter" on a `Molecule` cell, "To Atomic
Level..." on a `Macromolecule` column header, "Sparklines..." on
selected numerical columns.

## Prerequisites

- A semantic type already attached to the target columns — built-in
  (`Molecule`, `Macromolecule`, …) or package-registered (see
  `register-identifiers`).

## Steps

1. **Pick the input type — it IS the invocation surface.** Four tokens dispatch on different surfaces (see DG-FACT-182):
   - `semantic_value <name> { semType: <T> }` → cell right-click + per-row Actions pane.
   - `column <name> { semType: <T> }` → column-header right-click.
   - `list <name> { type: numerical }` → multi-column selection.
   - `string <name> { semType: <T> }` → article form; loses cell/column context — avoid.

2. **Declare the action with `meta.action: <Display Text>`.** The annotation IS the registration; exactly one input carrying `semType` is required (see DG-FACT-180, DG-FACT-181).

   ```typescript
   // src/package.ts — canonical decorator form
   export class PackageFunctions {
     @grok.decorators.func({
       name: 'Use as filter',
       description: 'Adds this structure as a substructure filter',
       meta: {action: 'Use as filter'},
     })
     static useAsSubstructureFilter(
       @grok.decorators.param({options: {semType: 'Molecule'}})
       value: DG.SemanticValue,
     ): void {
       const tv = grok.shell.tv;
       if (tv == null) throw new Error('Requires an open table view.');
       const molCol = value.cell?.column;
       if (molCol == null) throw new Error('Molecule column not found.');
       // build molblock from value.value, branching on units (step 4)
       tv.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
         type: DG.FILTER_TYPE.SUBSTRUCTURE,
         column: molCol.name,
         columnName: molCol.name,
         molBlock: /* … */ '',
       }, false);
     }
   }
   ```

   Header-comment form is equivalent (see DG-FACT-184).

3. **Scope WHERE the action surfaces with the `meta.exclude-*` toggles.** Both surfaces light up by default; the `exclude-actions-panel` / `exclude-current-value-menu` sibling toggles selectively hide one (see DG-FACT-183).

   ```typescript
   @grok.decorators.func({
     name: 'Copy as SMILES',
     meta: {action: 'Copy as SMILES', 'exclude-actions-panel': 'true'},
   })
   static copyAsSmiles(/* @param ... value: DG.SemanticValue */) { /* … */ }
   ```

4. **For "Use as filter" on `Molecule` — branch on units, then `updateOrAdd`.** Use `.updateOrAdd(filterDescriptor, false)` (not `.add(...)`) so a second invocation REPLACES the existing filter row; branch on `value.cell.column.meta.units` to convert SMILES via `DG.chem.convertMolNotation` and other notations via `molToMolblock` (see DG-FACT-185).

5. **Build and publish — registration is automatic.**
   ```bash
   npm install && npm run build   # regenerates src/package.g.ts
   grok publish <host>            # add --release once stable
   ```

## Common failure modes

- **Action never appears.** `meta.action` missing OR input lacks `semType` — check `src/package.g.ts` for both (see DG-FACT-180, DG-FACT-181).
- **Wrong invocation surface / fires on the wrong column.** Wrong input type. Switch to `semantic_value` and read `value.cell?.column` so you get the column the user clicked (see DG-FACT-182).
- **Article snippet won't compile / `FILTER_TYPE` is undefined.** Add `import * as DG from 'datagrok-api/dg';` and use `DG.FILTER_TYPE.SUBSTRUCTURE`.
- **Each invocation stacks a new filter row.** Use `.updateOrAdd(...)` not `.add(...)` (see DG-FACT-185).
- **`Copy as <format>` sub-actions clutter the Actions pane.** Add `meta.exclude-actions-panel: 'true'` to each variant; `exclude-current-value-menu: 'true'` on the parent chooser (see DG-FACT-183).

## See also

- Source articles: `help/develop/how-to/ui/context-actions.md`.
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-180`/`181`/`182`/`183`/`184`/`185`.
- Reference packages:
  - `packages/Chem/src/package.ts:2014-2087` + `package.g.ts:843-905` —
    `useAsSubstructureFilter` and the `Copy as...` family (both
    `exclude-*` toggles).
  - `packages/Bio/src/package.g.ts:391-396` — `column`-input action
    (`To Atomic Level...`, `semType: Macromolecule`).
  - `packages/PowerGrid/src/package.g.ts:171-183` — `list`-input
    actions on multi-column selection.
- Related skills: `register-identifiers` (registers the `semType` this
  skill dispatches on); `add-info-panel` (sibling — same `semType`
  dispatch, renders a widget instead of a menu verb).
