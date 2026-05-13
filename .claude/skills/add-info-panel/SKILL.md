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

1. **Pick the implementation path.** Client-side JS panel function
   (runs in the browser, direct DOM/`DG.Widget` access) or server-side
   script in Python/R/Grok (runs on the server, useful for heavy compute
   or non-JS libraries). Both register identically via `meta.role: panel`
   and re-evaluate visibility on every context change (`DG-FACT-276`).
   The article shows both at `help/develop/how-to/ui/add-info-panel.md:13-22`.

2. **Declare the panel (JS path) on `PackageFunctions`.** Use
   `@grok.decorators.panel({...})` with exactly one input — typically
   `semantic_value` (or `column` / `dataframe` / `row` / `string` /
   `file`) carrying the `semType` constraint — and return
   `Promise<DG.Widget>`. A panel function accepts EXACTLY ONE primary
   input (`DG-FACT-280`); extra user/dataset filters go in `condition`.

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

   `semantic_value` gives `.value`, `.cell`, `.gridCell`, `.column`,
   `.semType`, `.units` (`DG-FACT-282`) — use it when the widget needs
   the cell's context. For column-wide panels declare `column col`; for
   table-wide, `dataframe table`. Null-check the input —
   `semantic_value` arrives empty when no cell is selected.

   Header-comment form (no decorator) works identically — top-level
   `export function` prefixed with `//meta.role: panel`,
   `//input: <type> <name>`, `//output: widget result`. The article
   shows it at `help/develop/how-to/ui/add-info-panel.md:96-105`;
   `packages/Bio/src/package.g.ts:50-58` is a production example.

3. **Add a visibility condition (optional).** `//condition:` is a
   Grok-script boolean re-evaluated on every context change
   (`DG-FACT-279`), in Grok-script regardless of the panel's
   implementation language (`help/develop/how-to/ui/add-info-panel.md:88`).
   Common idioms beyond the `semType:` input constraint:

   - **User role.** `condition: user.hasrole("chemist")`.
   - **Column stats.** `condition: x.isnumerical && x.stats.missingvaluecount > 0`.
   - **Cross-package predicate.** `condition: Boltz1:isApplicableBoltz(molecule)` —
     declare the predicate as a sibling `//output: bool result` function;
     see `packages/Boltz1/src/package.g.ts:44-58` (`DG-FACT-283`).

4. **Or declare a panel script (server-side).** A script under
   `scripts/<name>.py` / `.r` / `.grok` uses the same role/condition
   keys with `#` prefixes; supported output types are `widget`,
   `viewer`, `graphics`, `dataframe` (`{action: join(<table>)}` merges
   results back as virtual columns), and `string` (`{action: markup}`
   for inline action buttons) (`DG-FACT-278`).

   ```python
   # scripts/spectrogram.grok
   #language: grok
   #meta.role: panel
   #input: column signal { type: numerical }
   #output: graphics pic
   #condition: signal.name == "F3"
   pic = Spectrogram("eeg", signal, 256.0, 1024, 0.1, true)
   ```

5. **Build, and verify the auto-emitted wrapper in `package.g.ts`.**
   `npm run build` runs the function-generator that scans
   `PackageFunctions` and rewrites `src/package.g.ts`. Do NOT hand-edit
   `package.g.ts`; it is overwritten on every build.

   ```bash
   npm install && npm run build
   ```

   Expected: build exits 0; a decorated panel declared with
   `meta: {role: 'widgets', ...}` emits `//meta.role: widgets,panel`
   — the generator auto-appends `panel` (`DG-FACT-281`). Compare
   `packages/Admetica/src/package.ts:30-40` against
   `packages/Admetica/src/package.g.ts:9-17`.

6. **Publish.** The `panel` role IS the registration — the platform
   discovers the function from the `package.g.ts` wrapper on load.

   ```bash
   grok publish <host>   # add --release once stable
   ```

   Expected: `grok publish <host>` exits 0 and reports the published
   version.

## Common failure modes

- **Panel never appears.** `meta.role` token is missing or misspelled.
  Inspect `src/package.g.ts` — it MUST contain `//meta.role: panel` or
  `//meta.role: <other>,panel` (lowercase). `Panel`/`panels` won't
  match (`DG-FACT-276`).
- **Panel fires for everything (or nothing).** The `semType` string on
  the input doesn't match what the column carries. Right-click the
  column header → *Properties* → check `semType`; in code, `col.semType`.
  Exact-match, case-sensitive (`DG-FACT-280`).
- **`Cannot read property 'cell' of null` at runtime.** The
  `semantic_value` input arrived empty (no current cell). Add an early
  `if (!value) return new DG.Widget(ui.divText('no selection'));`
  guard (`DG-FACT-282`).
- **Output type rejected.** A JS panel function MUST declare
  `//output: widget result` and return a `DG.Widget` (`DG-FACT-277`);
  scripts may instead emit `viewer` / `graphics` / `dataframe` (with
  `{action: join(table)}`) / `string` (with `{action: markup}`)
  (`DG-FACT-278`). Wrap raw `HTMLElement` or `DG.Viewer` as
  `new DG.Widget(el)` or `DG.Widget.fromRoot(viewer.root)`.

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
