---
name: register-identifiers
description: Register an identifier pattern + handler so Datagrok detects and explores domain IDs platform-wide
---

# register-identifiers

## When to use

Your package owns a class of domain identifiers (CHEMBL IDs, supplier
SKUs, drug codes) and you want Datagrok to detect them anywhere — grids,
free text, **Search Everywhere** — and render hover/context-panel content
for them. Triggers: "make my IDs clickable", "show a structure when I
hover a CHEMBL", "drill into related rows in our DB". For DB-backed IDs
prefer the DB-explorer path (Step 4) — the article calls it the
**recommended method** (knowledge `DG-FACT-064`).

## Prerequisites

- A package scaffold (`grok create <Name>`); commands run from the package root.
- `package.json` editable (basic path declares patterns there — knowledge
  `DG-FACT-058`). DB-explorer adds `@datagrok-libraries/db-explorer`.
- Familiarity with `DG.SemanticValue` / `DG.ObjectHandler` (`DG-FACT-059`,
  `DG-FACT-061`).
- For DB-explorer: a working data connection in Datagrok plus the schema
  name. See the `access-data` skill for connection JSON grammar (`DG-FACT-033`).

## Steps

1. **Declare the pattern in `package.json` (basic path).**
   `meta.semanticTypes` lets the platform detect identifiers in any text
   without code (`DG-FACT-058`). Each entry needs `semType`, `description`,
   and a `parsers` ARRAY (multiple regexes per type).
   ```json
   {
     "name": "@datagrok/chembl",
     "meta": {
       "semanticTypes": [{
         "semType": "CHEMBL_ID",
         "description": "Compound id in the CHEMBL database",
         "parsers": [{"regexp": "CHEMBL\\d+"}]
       }]
     }
   }
   ```
   Expected: after `grok publish`, the platform highlights `CHEMBL1234`
   literals as tagged tokens (compare `packages/Chemspace/package.json:67-77`).

2. **Write an `ObjectHandler` for hover/card/properties rendering.**
   Override `get type()`, `isApplicable(x)`, and at least one of
   `renderTooltip` / `renderCard` / `renderProperties` (`DG-FACT-059`).
   Canonical `isApplicable`: `x instanceof DG.SemanticValue && x.semType
   === '<TYPE>'`; the raw id is `(x as DG.SemanticValue).value`
   (`DG-FACT-061`; article omits the `DG.` prefix — drift
   `DG-FACT-DRIFT-023`).
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
   Expected: `src/handlers.ts` compiles; class extends `DG.ObjectHandler`.

3. **Register the handler from an autostart init.**
   `DG.ObjectHandler.register(...)` is the static entry point
   (`DG-FACT-060`). The init MUST be autostart — decorator form is
   canonical (`DG-FACT-063`). Drop the article's `init,` token in
   `//meta.role: init, autostart`; only `autostart` is documented (drift
   `DG-FACT-DRIFT-021`).
   ```typescript
   // src/package.ts — imports: DG, grok, ChemblIdHandler from './handlers'
   export const _package = new DG.Package();
   export class PackageFunctions {
     @grok.decorators.autostart()
     static init() { DG.ObjectHandler.register(new ChemblIdHandler()); }
   }
   ```
   Expected: after `grok publish --release`, hovering `CHEMBL1234` on any
   view triggers your `renderTooltip`.

4. **(Recommended for DB-backed IDs) Declare a DB-explorer config.**
   Replaces Steps 2–3 with declarative config. The library auto-registers
   the regex detector via `SemanticValue.registerRegExpDetector`
   (`DG-FACT-062`), builds handlers, and renders drill-downs from DB
   foreign keys (`DG-FACT-064`, `DG-FACT-065`). Use the colon-form
   `nqName: '<Package>:<ShortName>'` — the bare-name form in the article
   doesn't trigger namespace-qualified lookup (drift `DG-FACT-DRIFT-020`).
   ```bash
   npm install @datagrok-libraries/db-explorer
   ```
   ```typescript
   // src/explorer-config.ts
   import {DBExplorerConfig} from '@datagrok-libraries/db-explorer/src/types';
   export const explorerConfig: DBExplorerConfig = {
     connectionName: 'CHEMBL', schemaName: 'public',
     dataSourceName: 'postgres', nqName: 'Chembl:Chembl',
     entryPoints: {
       CHEMBL_ID: {
         table: 'molecule_dictionary', column: 'chembl_id',
         regexpExample: {example: 'CHEMBL1234',
           nonVariablePart: 'CHEMBL', regexpMarkup: 'CHEMBL[0-9]+'},
         matchRegexp: 'CHEMBL\\d+',
       },
     },
     joinOptions: [{
       fromTable: 'molecule_dictionary', columnName: 'molregno',
       tableName: 'compound_structures', onColumn: 'molregno',
       select: ['canonical_smiles', 'standard_inchi'],
     }],
     customRenderers: [{table: 'compound_structures',
       column: 'canonical_smiles', renderer: 'molecule'}],
   };
   ```
   Expected: `entryPoints` keys are semantic-type names; each value carries
   `table`+`column` (`DG-FACT-066`). `customRenderers[].renderer` accepts
   only `molecule` | `helm` | `imageURL` | `rawImage` — exact match
   (`DG-FACT-069`).

5. **Initialize DB-explorer from autostart.**
   `DBExplorer.initFromConfig(...)` returns the instance OR `null` —
   null-check is mandatory (`DG-FACT-064`). For content-based renderer
   overrides (e.g. SMILES validation) use `addCustomRenderer`
   (`DG-FACT-070`). Mirror `packages/Chembl/src/handlers.ts`.
   ```typescript
   // src/handlers.ts — imports: grok; DBExplorer, moleculeRenderer from
   // '@datagrok-libraries/db-explorer/src/{db-explorer,renderer}'; explorerConfig
   export function registerChemblIdHandler() {
     const exp = DBExplorer.initFromConfig(explorerConfig);
     if (!exp) { grok.shell.error('Failed to load db-explorer config'); return; }
     exp.addCustomRenderer((_, colName, value) => {
       const lc = colName?.toLowerCase() || '';
       return (lc === 'structure' || lc.includes('smiles')) &&
         typeof value === 'string' && grok.chem.checkSmiles(value);
     }, (value) => moleculeRenderer(value as string));
   }
   ```
   Wire `registerChemblIdHandler()` from the same `@autostart static init()`
   used in Step 3. Expected: after publish, clicking a `CHEMBL1234` cell
   opens the **Context Panel** (F4) with the primary row plus joined
   structures rendered as molecules.

6. **(No-code alternative) Configure from the connection UI.**
   For one-off configs without a plugin: **Browse → Databases →
   \<connection\> → right-click → Configure Identifiers…**, pick the
   schema, fill the accordion (Identifiers / Joins / Explicit References /
   Header Names / Unique Columns / Custom Selected Columns / Renderers),
   Save (`DG-FACT-071`). Use Import/Export JSON at the bottom to promote
   between environments. Expected: **Save** validates required fields and
   persists the config.

## Common failure modes

- **Handler never fires after publish.** The init function isn't
  autostart. Use `@grok.decorators.autostart()` on `static init()`, OR
  header `//meta.role: autostart` — drop the `init,` prefix the article
  shows (drift `DG-FACT-DRIFT-021`).
- **`x.semType` undefined inside `isApplicable`.** Detector regex never
  fired — recheck `meta.semanticTypes[].parsers[].regexp` (`DG-FACT-058`)
  and confirm the package republished after the `package.json` change.
- **`DBExplorer.initFromConfig(...)` returns `null`.** Schema couldn't be
  loaded or connection wasn't found (`DG-FACT-064`). Verify
  `connectionName` (case-insensitive — `DG-FACT-033`), use the colon-form
  `nqName` (drift `DG-FACT-DRIFT-020`), confirm the connection is deployed.
- **Custom renderer doesn't render the structure.** `customRenderers[].renderer`
  is an exact-match enum (`DG-FACT-069`). Anything else (`smiles`, `Mol`)
  silently falls through to text. For content-based detection use
  `addCustomRenderer` (Step 5).
- **Cross-schema join returns nothing.** `fromSchema`/`onSchema` were
  omitted on a join that crosses schemas (`DG-FACT-067`). Both keys are
  optional only when the join stays inside `schemaName`.

## Verification

- `grok publish <host> --release` exits `0`.
- Open a view with a known identifier — Datagrok highlights it as a tagged
  token; hover triggers your tooltip.
- Click the value: **Context Panel** (F4) shows the primary row + joined
  data (DB-explorer path).
- Type the identifier into **Search Everywhere** — a card with the
  rendered structure appears in results.

## See also

- Source articles: `help/develop/how-to/db/register-identifiers.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-058` through `DG-FACT-072` and drifts `DG-FACT-DRIFT-020..023`.
- Related skills:
  - `access-data` (sibling — connection JSON, query namespacing,
    `grok.data.query` used inside handlers).
  - `db-in-plugin` (sibling — when the IDs live in a Postgres DB you
    own; produces the connection that DB-explorer points at).
