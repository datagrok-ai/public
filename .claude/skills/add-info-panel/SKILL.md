---
name: add-info-panel
version: 0.1.0
description: |
  Register a Datagrok package function whose output renders in the
  right-hand context panel whenever the active selection matches a
  declared semantic type or visibility condition. For plugin authors
  who want extra computed properties, mini-viewers, or action buttons
  to appear automatically beside the current cell, column, row, or
  table without the user invoking anything. Produces a `panel`-role
  function (JS for client-side, script for server-side) plus the
  auto-emitted `package.g.ts` wrapper the platform reads at startup.
  Use when asked to "show computed properties beside the selected
  molecule", "attach a widget to compound cells in the context panel",
  or "react to selection changes with a server-computed result".
triggers:
  - show widgets in the right side panel
  - react to selection in the context panel
  - attach computed properties to a semantic type
  - panel that reacts when a cell changes
  - show enrichments beside the current row
  - server-side widget bound to a column type
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# add-info-panel

## When to use

Your package needs a widget in the right-hand context panel whenever
the user's selection (cell, column, row, table, file) matches a
declared semantic type or condition — without an explicit invocation.
Real-world phrasings: "show ADMET predictions next to a `Molecule`
cell", "render a sequence logo for `Macromolecule` columns", "add a
`Flag as suspicious` button beside each row".

## Prerequisites

- A semantic type already attached to the column — built-in
  (`Molecule`, `Macromolecule`, `Text`, …) or one your package
  registers via `register-identifiers` /
  `define-semantic-type-detectors` (`DG-FACT-276`).
- Familiarity with `DG.Widget` — wrap an `HTMLElement` via
  `new DG.Widget(el)` or `DG.Widget.fromRoot(el)` (`DG-FACT-277`).

## Steps

1. **Pick the implementation path.** Client-side JS (decorator or header form) or server-side script — both register identically via `meta.role: panel` (see DG-FACT-276).

2. **Declare the panel (JS path) on `PackageFunctions`.** Use `@grok.decorators.panel({...})` with exactly one primary input (see DG-FACT-280) and return `Promise<DG.Widget>` (see DG-FACT-277). Null-check the input — `semantic_value` arrives empty when no cell is selected (see DG-FACT-282).

   ```typescript
   // src/package.ts
   export class PackageFunctions {
     @grok.decorators.panel({
       name: 'Compound | Properties',
       meta: {role: 'widgets', domain: 'chem'},
     })
     static async compoundPropertiesPanel(
       @grok.decorators.param({options: {semType: 'Molecule'}})
       value: DG.SemanticValue,
     ): Promise<DG.Widget> {
       if (!value)
         return new DG.Widget(ui.divText('no selection'));
       return new DG.Widget(ui.divText(`column: ${value.cell.column.name}`));
     }
   }
   ```

   Header-comment form (no decorator) works identically — top-level `export function` prefixed with `//meta.role: panel`, `//input: <type> <name>`, `//output: widget result`.

3. **Add a visibility condition (optional).** `//condition:` is a Grok-script boolean re-evaluated on every context change (see DG-FACT-279). Common idioms:

   - **User role.** `condition: user.hasrole("chemist")`.
   - **Column stats.** `condition: x.isnumerical && x.stats.missingvaluecount > 0`.
   - **Cross-package predicate.** `condition: Boltz1:isApplicableBoltz(molecule)` (see DG-FACT-283).

4. **Or declare a panel script (server-side).** A script under `scripts/<name>.py` / `.r` / `.grok` uses the same role/condition keys with `#` prefixes; supported output types are `widget`, `viewer`, `graphics`, `dataframe`, `string` (see DG-FACT-278 for action: directives).

   ```python
   # scripts/spectrogram.grok
   #language: grok
   #meta.role: panel
   #input: column signal { type: numerical }
   #output: graphics pic
   #condition: signal.name == "F3"
   pic = Spectrogram("eeg", signal, 256.0, 1024, 0.1, true)
   ```

5. **Build, and verify the auto-emitted wrapper in `package.g.ts`.** `npm run build` runs the function-generator that rewrites `src/package.g.ts`. Do NOT hand-edit it. The generator auto-appends `panel` to non-`panel` role strings (see DG-FACT-281).

   ```bash
   npm install && npm run build
   ```

6. **Publish.** The `panel` role IS the registration — the platform discovers the function from the `package.g.ts` wrapper on load.

   ```bash
   grok publish <host>   # add --release once stable
   ```

## Common failure modes

- **Panel never appears.** `meta.role` token missing/misspelled — check `src/package.g.ts` for lowercase `//meta.role: panel` (see DG-FACT-276).
- **Panel fires for everything (or nothing).** The `semType` string on the input doesn't match what the column carries (exact, case-sensitive — see DG-FACT-280).
- **`Cannot read property 'cell' of null`.** Missing null guard on `semantic_value` input (see DG-FACT-282).
- **Output type rejected.** JS panels MUST return `DG.Widget` (see DG-FACT-277); scripts have a wider output set (see DG-FACT-278).

## See also

- Source article: `help/develop/how-to/ui/add-info-panel.md`.
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-276/277/278/279/280/281/282/283`.
- Reference packages: `packages/Bio/src/package.g.ts:50-58`
  (header-form JS panel, `Macromolecule`);
  `packages/Admetica/src/package.ts:30-40` +
  `packages/Admetica/src/package.g.ts:9-17` (decorator form, source
  vs. wrapper); `packages/Boltz1/src/package.g.ts:44-58`
  (`semantic_value` + cross-package predicate `condition`).
- Related skills: `register-identifiers` and
  `define-semantic-type-detectors` (define the `semType` this skill
  binds to); `home-page-widgets` (similar `DG.Widget` return contract,
  unbound to selection).
