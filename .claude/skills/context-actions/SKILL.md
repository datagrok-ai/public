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

1. **Pick the input type — it IS the invocation surface (`DG-FACT-182`).**
   Four tokens dispatch on different surfaces; the article shows only
   `string`, which loses cell context.
   - `semantic_value <name> { semType: <T> }` → `DG.SemanticValue`;
     **cell** right-click + per-row Actions pane. Default
     (`packages/Chem/src/package.g.ts:843-849`).
   - `column <name> { semType: <T> }` → `DG.Column`; **column-header**
     right-click (`packages/Bio/src/package.g.ts:391-396`).
   - `list <name> { type: numerical }` → `DG.Column[]`;
     **multi-column selection** (`packages/PowerGrid/src/package.g.ts:171-176`).
   - `string <name> { semType: <T> }` → article form; raw value only,
     no cell/column reference. Production never uses it here.

2. **Declare the action with `meta.action: <Display Text>`.**
   Adding `meta.action` IS the registration — no explicit `register(...)`
   call (`DG-FACT-180`). The text becomes the menu label. Exactly one
   input, carrying `semType`, is required (`DG-FACT-181`); without
   `semType` the function registers but has no dispatch surface.

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

   Expected: after `npm run build`, the regenerated `src/package.g.ts`
   carries a wrapper with `//name: Use as filter`,
   `//input: semantic_value value { semType: Molecule }`, and
   `//meta.action: Use as filter` (compare
   `packages/Chem/src/package.g.ts:843-849`). The header-comment form
   shown in the article is equivalent (`DG-FACT-184`).

3. **Scope WHERE the action surfaces with the `meta.exclude-*` toggles.**
   Both surfaces light up by default. Suppress one with a sibling
   `meta.*: 'true'` flag — the article mentions neither (`DG-FACT-183`):
   - `meta.exclude-actions-panel: 'true'` → keeps right-click, hides
     from Actions pane. Chem uses it on each `Copy as <format>`
     sub-action so the panel only shows the parent chooser
     (`packages/Chem/src/package.ts:2066,2079`).
   - `meta.exclude-current-value-menu: 'true'` → keeps Actions pane,
     hides from right-click. Chem uses it on the parent `Copy as...`
     chooser (`packages/Chem/src/package.ts:2049`).

   ```typescript
   @grok.decorators.func({
     name: 'Copy as SMILES',
     meta: {action: 'Copy as SMILES', 'exclude-actions-panel': 'true'},
   })
   static copyAsSmiles(/* @param ... value: DG.SemanticValue */) { /* … */ }
   ```

4. **For "Use as filter" on `Molecule` — branch on units, then `updateOrAdd`.**
   The article uses `.add(...)` and bare `FILTER_TYPE.SUBSTRUCTURE`;
   production uses `.updateOrAdd(filterDescriptor, false)` so a second
   invocation REPLACES the existing filter row instead of stacking, plus
   fully-qualified `DG.FILTER_TYPE.SUBSTRUCTURE` (`DG-FACT-185`). When
   `value.cell.column.meta.units == DG.chem.Notation.Smiles`, convert
   via `DG.chem.convertMolNotation` to preserve orientation; otherwise
   `molToMolblock(molecule, getRdKitModule())`
   (`packages/Chem/src/package.ts:2032-2036`).

5. **Build and publish — registration is automatic.**
   ```bash
   npm install && npm run build   # regenerates src/package.g.ts
   grok publish <host>            # add --release once stable
   ```
   The `meta.action` annotation IS the registration. After deploy,
   right-click an item whose `semType` matches the declared input to
   see the action.

## Common failure modes

- **Action never appears on right-click OR Actions pane.** The
  `meta.action` token is missing or the input lacks a `semType`.
  Inspect the auto-generated `src/package.g.ts`: it MUST carry both
  `//meta.action: <Text>` AND `//input: ...{ semType: <T> }`
  (`DG-FACT-180`, `DG-FACT-181`).
- **Action surfaces on the wrong invocation surface** (e.g. column
  header but not cells) **or fires on the wrong column when two columns
  share a `semType`.** Wrong input type: `column` fires on the header,
  `semantic_value` on the cell, `list` on multi-column selection. A
  `string` body that calls `tv.dataFrame.columns.bySemType(...)` picks
  the first match, not the column the user clicked — switch the input
  to `semantic_value` and read `value.cell?.column` (`DG-FACT-182`).
- **Article snippet won't compile / `FILTER_TYPE` is undefined.** The
  article omits imports and uses bare `FILTER_TYPE.SUBSTRUCTURE`. Add
  `import * as DG from 'datagrok-api/dg';` and use
  `DG.FILTER_TYPE.SUBSTRUCTURE`.
- **Each invocation of "Use as filter" stacks a new filter row.** The
  article uses `.add(...)`; switch to
  `.updateOrAdd(filterDescriptor, false)` so the existing row for that
  column is replaced (`DG-FACT-185`).
- **`Copy as <format>` sub-actions clutter the Actions pane.** Each
  format variant lacks `meta.exclude-actions-panel: 'true'`. Add the
  flag to every variant; keep `meta.exclude-current-value-menu: 'true'`
  on the parent `Copy as...` chooser (`DG-FACT-183`).

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
