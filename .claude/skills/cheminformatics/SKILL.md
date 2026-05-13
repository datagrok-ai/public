---
name: cheminformatics
description: Use the grok.chem JS API for substructure search, similarity, MCS, R-group, descriptors, and molecule rendering
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# cheminformatics

## When to use

You're writing a Datagrok package or script that needs to operate on
molecules — substructure or similarity search, picking diverse subsets,
finding the maximum common substructure, R-group analysis, computing
descriptors, or rendering structures into a panel/dialog. Triggers:
"filter rows by SMARTS", "rank by Tanimoto", "compute MolWt", "draw a
molecule into a `div`/`canvas`".

## Prerequisites

- A package scaffold (`grok create <Name>`) or an Apps-style script that
 runs against a Datagrok host. Imports come from `datagrok-api`:
 `grok.chem.*`, `ui.*`, `DG.*` (knowledge `DG-FACT-235`).
- A column of String type holding molecules in one of the five RDKit JS
 notations: `smiles`, `cxsmiles`, `molblock`, `v3Kmolblock`, `inchi`
 (knowledge `DG-FACT-236`). Format detection is automatic; molfiles are
 recognised by the `M END` marker.
- Most `grok.chem.*` calls are async — every step below uses `await`.

## Steps

1. **Substructure search → BitSet.**
 ```typescript
 import * as grok from 'datagrok-api/grok';

 const t = await grok.data.loadTable(
 'https://public.datagrok.ai/demo/sar_small.csv');
 const bs = await grok.chem.searchSubstructure(
 t.col('smiles')!,
 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1');
 t.filter.copyFrom(bs);
 grok.shell.addTableView(t);
 ```
 Returns a `BitSet` aligned 1:1 with the column; unparseable molecules
 leave bit `0` (no throw). Settings is `{ molBlockFailover?: string }`
 only (see `DG-FACT-238`).

2. **Find most similar molecules.**
 ```typescript
 const mols = await grok.data.loadTable(
 'https://public.datagrok.ai/demo/sar_small.csv');
 const top = await grok.chem.findSimilar(
 mols.col('smiles')!,
 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1',
 { limit: 20, cutoff: 0.4 });
 grok.shell.addTableView(top!);
 ```
 Returns a `DataFrame` of `molecule` / `score` / `index` columns,
 sorted descending by Tanimoto score; Morgan fingerprints (see
 `DG-FACT-237`, `DG-FACT-239`).

3. **Score every row, or prime the fingerprint cache.**
 ```typescript
 const scores = await grok.chem.getSimilarities(
 t.col('smiles')!, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1');
 t.columns.add(scores!);
 // Cache-priming form (returns null, no scoring):
 await grok.chem.getSimilarities(t.col('smiles')!, '');
 ```
 Returns a single `Column` (not a DataFrame). Pass `molecule = ''`
 (not `null`) to prime the cache; returns `null` in that branch (see
 `DG-FACT-240`).

4. **Pick a diverse subset.**
 ```typescript
 const diverse = await grok.chem.diversitySearch(
 t.col('smiles')!, { limit: 20 });
 grok.shell.addTableView(diverse);
 ```
 Expected: a `DataFrame` of up to `limit` representative molecules.
 The actual signature is `(column, settings = {limit: MAX_VALUE})` —
 there is NO `metric` parameter and NO `METRIC_TANIMOTO` constant
 despite what the article shows.

5. **Compute the Maximum Common Substructure (MCS).**
 ```typescript
 const smarts = await grok.chem.mcs(
 t, 'smiles', /*returnSmarts*/ true);
 grok.shell.info(smarts);
 ```
 Signature is `(DataFrame, columnName, returnSmarts?, exactAtomSearch?,
 exactBondSearch?)`. Passing a Column will fail at runtime.

6. **R-group analysis with a scaffold SMARTS.**
 ```typescript
 const rg = await grok.chem.rGroup(
 t, 'smiles', 'O=C1Nc2ccccc2C(C2CCCCC2)=NC1');
 grok.shell.addTableView(rg);
 ```
 Returns a new DataFrame, one column per R-group position (`R1`, `R2`, …).

7. **Compute descriptors.**
 ```typescript
 const desc = await grok.chem.descriptors(
 grok.data.testData('molecules', 100), 'smiles',
 ['MolWt', 'NumHAcceptors', 'NumHDonors', 'TPSA']);
 grok.shell.addTableView(desc);
 ```
 Lipinski / Crippen / MolSurf / Fragments group names also accepted;
 catalogue in `packages/ApiSamples/scripts/domains/chem/descriptors.js`.

8. **Render a molecule into a `div` (sync, OpenChemLib).**
 ```typescript
 import * as ui from 'datagrok-api/ui';
 ui.dialog
.add(grok.chem.svgMol(
 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1', 300, 200))
.show;
 ```
 `svgMol` is sync — returns the div immediately, SVG fills in on the
 next microtask (see `DG-FACT-241`).

9. **Render to a `canvas` with scaffold highlighting (async, RDKit).**
 ```typescript
 const canvas = document.createElement('canvas');
 canvas.width = 300; canvas.height = 200;
 await grok.chem.canvasMol(
 0, 0, canvas.width, canvas.height, canvas,
 'COc1ccc2cc(ccc2c1)C(C)C(=O)OCCCc3cccnc3',
 'c1ccccc1');
 ui.dialog({ title: 'Molecule' }).add(canvas).show;
 ```
 `canvasMol` is async; scaffold is SMARTS (see `DG-FACT-241`).

## Common failure modes

- `searchSubstructure` missing rows — unparseable molecules silently
 stay `0` (`DG-FACT-238`). Use `molBlockFailover` SMARTS.
- `mcs(column)` — wrong signature; pass `(DataFrame, columnName, …)`.
- `diversitySearch(col, METRIC_TANIMOTO, 10)` — constant doesn't exist;
 second arg is `{ limit }`.
- `getSimilarities(col, null)` — cache-priming wants `''`, not `null`.
- `canvasMol` shows nothing — must `await` (`DG-FACT-241`).
- `svgMol` empty — append div first, don't read `innerHTML` sync.

## See also

- Source articles:
 - `help/develop/domains/chem/cheminformatics.md`
- Knowledge:
 - `docs/_internal/knowledge/knowledge-graph.md` — facts `DG-FACT-235`
 through `DG-FACT-241`..094`.
- Reference samples:
 - `packages/ApiSamples/scripts/domains/chem/{substructure-search,
 similarity-scoring-sorted,descriptors,r-group,mcs,diversity-search,
 mol-rendering,mol-rendering-canvas}.js`
- Related skills:
 - `custom-cell-renderers` (rendering molecules as table cells uses the
 same `svgMol`/`canvasMol` engines).
 - `access-data` (load a `DataFrame` of molecules from SQL/CSV before
 feeding it to `grok.chem.*`).
