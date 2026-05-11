---
name: add-info-panel
description: Register a Datagrok info panel that appears in the context panel when a visibility condition matches a column, value, row, or table
harness-authored: true
---

# add-info-panel

## When to use

Your package needs to push extra context into the right-hand
**Context Panel** whenever the user lands on a relevant object — a
column with a specific `semType`, a `semantic_value` cell, a row from
a particular table, or a dataframe matching some shape. Triggers:
"show ADMET predictions for `Molecule` cells", "render a composition
widget for `Macromolecule` columns", "surface a 'Flag as suspicious'
action for transaction rows".

## Prerequisites

- A package scaffold (`grok create <Name>`); commands run from the package root.
- `datagrok-api` available (`import * as grok / DG / ui`) — the article's
  snippets omit imports and won't compile as written.
- For `semType`-gated panels: the semantic type already attached to the
  column (built-in like `Molecule`/`Macromolecule`/`Text`, or one your
  package registers — see `register-identifiers`,
  `define-semantic-type-detectors`).
- Familiarity with `DG.Widget` for the JS path (`new DG.Widget(htmlEl)`
  / `DG.Widget.fromRoot(htmlEl)`), or with Grok-script header parameters
  (`# name:`, `# input:`, `# output:`, `# condition:`) for scripts.

## Steps

1. **Pick the implementation path** (`DG-FACT-276`).
   - **Client-side JS/TS panel function** — a `static` method on
     `PackageFunctions` decorated with `@grok.decorators.panel({...})`,
     OR an exported function annotated with `//meta.role: panel`.
     Body runs in the browser. Default for most plugins.
   - **Server-side panel script** — `.py` / `.R` / `.grok` under
     `scripts/` with a header block `# meta.role: panel`. Body runs on
     the Datagrok server; use for heavy compute or server-only libs.

2. **Declare EXACTLY ONE primary input** — its type narrows the
   binding (`DG-FACT-280`). Pick the narrowest match:
   - `column <name> { semType: <X> }` — active column (Bio
     `getRegionPanel` for `Macromolecule`).
   - `semantic_value <name> { semType: <X> }` — JS-only; active cell
     value, gives `value.value`, `value.cell`, `value.cell.column`,
     `value.gridCell`, `value.semType` (`DG-FACT-282`). Used by
     Admetica, Bio `monomerInfoPanel`, Bio `compositionAnalysisWidget`.
   - `dataframe <name>` — active table (`condition: table.name == "demog"`).
   - `row <name>` — active row (`condition: activity.table.name == "..."`).
   - `cell`, `file`, `string`, `int`, `object` — less common.

   Extra user/dataset filters MUST go in `condition`, NOT as additional
   inputs (`DG-FACT-280`).

3. **Write the JS panel function (decorator form, canonical).**
   ```typescript
   // src/package.ts
   import * as grok from 'datagrok-api/grok';
   import * as DG   from 'datagrok-api/dg';
   import * as ui   from 'datagrok-api/ui';
   export const _package = new DG.Package();

   export class PackageFunctions {
     @grok.decorators.panel({
       name: 'Chem | Molecule Properties',
       meta: {role: 'widgets', domain: 'chem'},  // DO NOT add 'panel' — codegen appends it (DG-FACT-281)
     })
     static moleculePanel(
       @grok.decorators.param({options: {semType: 'Molecule'}})
       sv: DG.SemanticValue
     ): DG.Widget<any> {
       if (!sv) return new DG.Widget(ui.divText('value is empty')); // DG-FACT-282: null-check
       return new DG.Widget(ui.divText(`column: ${sv.cell.column.name}`));
     }
   }
   ```
   Expected: build succeeds; `src/package.g.ts` gains a wrapper with
   `//input: semantic_value sv { semType: Molecule }`,
   `//output: widget result`, `//meta.role: widgets,panel` (no space,
   `panel` last — `DG-FACT-DRIFT-AIP-001`). Compare
   `packages/Admetica/src/package.g.ts:7-16`,
   `packages/Bio/src/package.g.ts:153-161`.

4. **OR write a server-side panel script.**
   ```python
   #name: Spectrogram info panel
   #language: grok
   #meta.role: panel
   #input: column signal {type:numerical}
   #output: graphics pic
   #condition: signal.name == "F3"
   pic = Spectrogram("eeg", signal, 256.0, 1024, 0.1, true)
   ```
   Allowed script output types (`DG-FACT-278`): `widget`, `viewer`,
   `graphics`, `dataframe` (typically `{action: join(<input-table>)}`
   to merge into the source frame), or `string {action: markup}` for
   inline action buttons (`#{button("Flag", "http.Post(...)")}`).

5. **Author the `condition` (Grok-script, regardless of language)** —
   re-evaluated on every context change (`DG-FACT-279`). Reference the
   input by name, plus platform-supplied `user` and `table`:
   ```text
   condition: smiles.semType == "Molecule"                       // DG-FACT-283 (a) — semType gate
   condition: x.isnumerical && x.stats.missingvaluecount > 0
   condition: table.name == "demog" && table.columns.containsAll(["height","weight"])
   condition: table.gettag("database") == "northwind"            // dataset filter
   condition: user.hasrole("chemist") || user.inteam("HTS")      // user filter
   condition: Boltz1:isApplicableBoltz(molecule)                 // DG-FACT-283 (b) — cross-pkg predicate
   ```
   Method names are case-insensitive but emit lowercase
   (`isnumerical`, `endswith`, `gettag` — `DG-FACT-DRIFT-AIP-002`).
   For `semType`, prefer camelCase to mirror the JS-API
   (`DG-FACT-DRIFT-AIP-003`).

6. **Build and publish.**
   ```bash
   npm install
   grok check                     # exits 0
   grok publish <host>            # add --release once stable
   ```
   No explicit `register(...)` — `meta.role: panel` IS the registration
   (`DG-FACT-276`). After deploy, put the platform into the panel's
   context: the panel appears in the **Context Panel** on the right.

## Common failure modes

- **Panel never appears.** Role token missing/misspelled. Inspect
  `src/package.g.ts`: MUST contain `//meta.role: panel` (or
  `<other>,panel` — codegen appends `panel`, `DG-FACT-281`). `Panel` /
  `panels` won't register.
- **Panel appears for wrong objects.** `condition` references a name
  that doesn't match the input, or `semType` is mistyped. Verify the
  input name matches the variable in `condition`, and `col.semType`
  (UI: column → *Properties*) matches the constraint exactly —
  case-sensitive (`DG-FACT-279`, `DG-FACT-280`).
- **TypeError on `value.cell.column`.** The `semantic_value` input was
  empty (no cell selected). Null-check before dereferencing
  (`DG-FACT-282`): `if (!value) return new DG.Widget(ui.divText('empty'));`.
- **Article snippet won't compile.** It shows bare
  `export function ...` with header comments and no imports. Translate
  to the decorator form on `PackageFunctions` with explicit imports;
  return `DG.Widget` (sync) or `Promise<DG.Widget>` (async)
  (`DG-FACT-277`).
- **Decorator declares `panel` twice.** Writing
  `meta: {role: 'panel,widgets'}` yields a doubled emit like
  `panel,widgets,panel`. Declare only the non-panel role
  (`DG-FACT-281`).
- **Script returns `widget` with no rendering backend, or JS returns
  `graphics`.** Scripts may return
  `widget`/`viewer`/`graphics`/`dataframe`/`string{action: markup}`
  (`DG-FACT-278`); JS panels MUST return `widget` (`DG-FACT-277`).

## Verification

- `grok check` exits `0`; `grok publish <host>` exits `0`.
- The regenerated `src/package.g.ts` contains, for each panel function,
  a wrapper with `//input: <type> <name> [{ semType: <X> }]`,
  `//output: widget result`, and `//meta.role: panel` (or
  `<role>,panel`). Compare `packages/Admetica/src/package.g.ts:7-16`,
  `packages/Bio/src/package.g.ts:50-58,153-161`.
- In Datagrok, put the platform into the panel's context (click a
  matching cell, select a matching column, etc.): the panel appears in
  the **Context Panel** on the right.
- Switch context to an object the `condition` excludes: the panel
  disappears — confirms the condition is re-evaluated on every change.

## See also

- Source articles:
  - `help/develop/how-to/ui/add-info-panel.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-276` … `DG-FACT-283`; drifts
  `DG-FACT-DRIFT-AIP-001` (codegen emits `<role>,panel`, no space),
  `DG-FACT-DRIFT-AIP-002` (Grok-script lowercase preferred),
  `DG-FACT-DRIFT-AIP-003` (camelCase `semType` preferred).
- Reference packages:
  - `packages/Admetica/src/package.ts:30-40` — decorator-form
    `admeticaWidget` with `semantic_value {semType: Molecule}`;
    codegen at `packages/Admetica/src/package.g.ts:7-16`.
  - `packages/Bio/src/package.g.ts:50-58` — header-form
    `getRegionPanel` with `column {semType: Macromolecule}`.
  - `packages/Bio/src/package.g.ts:153-161` —
    `compositionAnalysisWidget` with multi-role
    `//meta.role: widgets,panel`.
- Related skills:
  - `column-tooltip` (sibling — same role-based pattern, role token
    `tooltip`).
  - `register-identifiers` / `define-semantic-type-detectors`
    (prerequisite — registers the `semType` a panel binds to).
  - `home-page-widgets` (sibling — `meta.role: dashboard`).
