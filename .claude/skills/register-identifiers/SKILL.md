---
name: register-identifiers
version: 0.2.0
description: |
  Teach Datagrok to recognise a class of domain identifier strings
  (compound codes, accession numbers, ticket keys) anywhere on the
  platform — grids, free text, **Search Everywhere** — and decorate
  them with hover popups and context-panel cards. Targets plugin
  authors whose tables already render identifiers as plain text and
  who want rendering driven by their own code, not a DB schema.
  Use when asked to "make our compound codes clickable", "show a hover
  card for accession numbers", or "tag our ticket strings as a known
  entity".
triggers:
  - make domain codes clickable
  - hover popup for compound ids
  - tag strings as a known entity
  - search everywhere for our ids
  - detect accession numbers in cells
  - render card for our identifier
allowed-tools:
  - Read
  - Edit
  - Write
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# register-identifiers

## When to use

Your package owns a class of identifier strings (CHEMBL IDs, supplier
SKUs, UniProt accessions, ticket keys) rendered as plain text, and you
want Datagrok to detect them anywhere — grids, free text,
**Search Everywhere** — then dispatch hover / card / context-panel
rendering through code you control. If the identifier maps to a
**database row** and you want declarative drill-down into joined
tables, prefer the sibling skill `db-explorer-identifiers` (uses
`@datagrok-libraries/db-explorer`, the article calls it the
**recommended method**, `DG-FACT-064`). For one-off config from the
admin UI, see the "No-code alternative" pointer at the bottom.

## Prerequisites

- `package.json` editable (`DG-FACT-058`).
- Familiarity with `DG.SemanticValue` and `DG.ObjectHandler`
  (`DG-FACT-059`, `DG-FACT-061`).

## Steps

1. **Declare the pattern in `package.json`.**
   `meta.semanticTypes` detects identifiers in any text without code
   (`DG-FACT-058`). Each entry needs `semType`, `description`, and a
   `parsers` ARRAY (multiple regexes per type).
   ```json
   {
     "name": "@datagrok/chembl",
     "meta": {"semanticTypes": [{
       "semType": "CHEMBL_ID",
       "description": "Compound id in the CHEMBL database",
       "parsers": [{"regexp": "CHEMBL\\d+"}]
     }]}
   }
   ```
   Expected: after `grok publish`, the platform highlights `CHEMBL1234`
   literals as tagged tokens. Shipped example:
   `packages/Chemspace/package.json:67-77`.

2. **Write an `ObjectHandler` for hover/card/properties rendering.**
   Override `get type()`, `isApplicable(x)`, and at least one of
   `renderTooltip` / `renderCard` / `renderProperties` (`DG-FACT-059`).
   Canonical guard: `x instanceof DG.SemanticValue && x.semType ===
   '<TYPE>'`; raw id is `(x as DG.SemanticValue).value` (`DG-FACT-061`).
   ```typescript
   // src/handlers.ts — imports: grok, ui, DG from 'datagrok-api/{grok,ui,dg}'
   export class ChemblIdHandler extends DG.ObjectHandler {
     get type(): string { return 'CHEMBL_ID'; }
     isApplicable(x: any): boolean {
       return x instanceof DG.SemanticValue && x.semType === 'CHEMBL_ID';
     }
     renderTooltip(x: any): HTMLElement {
       const id = (x as DG.SemanticValue).value;
       return ui.divV([ui.h3(id), ui.wait(async () => ui.bind(x,
         grok.chem.drawMolecule(await grok.functions.call(
           'Chembl:chemblIdToSmiles', {id}))))]);
     }
   }
   ```
   Expected: `src/handlers.ts` compiles. Shipped reference:
   `packages/JiraConnect/src/jira-grid-cell-handler.ts:24-80`.

3. **Register the handler from an autostart init.**
   `DG.ObjectHandler.register(...)` is the static entry point
   (`DG-FACT-060`); must run from autostart — decorator form canonical
   (`DG-FACT-063`).
   ```typescript
   // src/package.ts — imports: DG, grok, ChemblIdHandler from './handlers'
   export const _package = new DG.Package();
   export class PackageFunctions {
     @grok.decorators.autostart()
     static init() { DG.ObjectHandler.register(new ChemblIdHandler()); }
   }
   ```
   Expected: after `grok publish --release`, hovering `CHEMBL1234`
   triggers your `renderTooltip`. Reference:
   `packages/JiraConnect/src/package.ts:24-32`.

## Common failure modes

- **Handler never fires after publish.** Init isn't autostart — use
  `@grok.decorators.autostart()` or `//meta.role: autostart`
  (`DG-FACT-063`).
- **`x.semType` undefined inside `isApplicable`.** Detector regex
  never fired — recheck `meta.semanticTypes[].parsers[].regexp`
  (`DG-FACT-058`) and confirm the package republished.
- **Hover fires but `renderTooltip` never runs.** `isApplicable`
  returned false — confirm `x instanceof DG.SemanticValue` guard and
  exact-match `semType` string (`DG-FACT-061`); a bare string with
  the right shape is NOT a `SemanticValue`.
- **Identifiers detected in grid cells but not in free text /
  Search Everywhere.** Regex too narrow or anchored — `parsers[].regexp`
  is matched against arbitrary text fragments, so avoid `^…$`
  anchors (`DG-FACT-058`).

## No-code alternative

For one-off identifier configs against an existing database connection
(no package code), use **Browse → Databases → \<connection\> →
right-click → Configure Identifiers…** and fill the accordion
(Identifiers / Joins / Explicit References / Header Names / Unique
Columns / Custom Selected Columns / Renderers), then **Save**
(`DG-FACT-071`). Import/Export JSON at the bottom moves config between
environments.

## See also

- Source: `help/develop/how-to/db/register-identifiers.md`
- Sibling skill (DB-backed IDs with declarative drill-down):
  `db-explorer-identifiers`.
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` —
  `DG-FACT-058`..`DG-FACT-063`, `DG-FACT-071`.
- Related skills: `access-data`, `db-in-plugin`.
