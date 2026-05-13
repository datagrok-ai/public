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
 Expected: a `BitSet` whose length equals `column.length`; bit `i` is
 `1` iff the `i`-th molecule contains the pattern, `0` for non-matches
 AND for unparseable molecules — invalid SMILES are silently skipped,
 not raised (knowledge `DG-FACT-238`). `settings` is
 `{ molBlockFailover?: string }` only — the `substructLibrary` flag
 shown in the article is not read by the wrapper.

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
 Expected: a `DataFrame` with three columns — `molecule` (original
 strings), `score` (Tanimoto, sorted descending, range `0.0..1.0`),
 `index` (positions in the input column). Defaults are
 `{ limit: Number.MAX_VALUE, cutoff: 0.0 }`, i.e. ALL molecules ranked
 (knowledge `DG-FACT-239`). Engine: Morgan fingerprints + Tanimoto
 (knowledge `DG-FACT-237`).

3. **Score every row, or prime the fingerprint cache.**
 ```typescript
 const scores = await grok.chem.getSimilarities(
 t.col('smiles')!, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1');
 t.columns.add(scores!);
 // Cache-priming form (returns null, no scoring):
 await grok.chem.getSimilarities(t.col('smiles')!, '');
 ```
 Expected: a single `Column` (NOT a DataFrame — drift
 ) of length equal to the input column; null
 scores for unparseable rows. With `molecule = ''` (not `null`) the
 per-column fingerprint cache is populated and the call returns
 `null` (knowledge `DG-FACT-240`, `DG-FACT-237`).

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
 Expected: a `Promise<string>` resolving to SMILES (or SMARTS when
 `returnSmarts=true`). The signature is
 `(table: DataFrame, column: string, returnSmarts?, exactAtomSearch?,
 exactBondSearch?)` — the article's `mcs(column)` form passes a
 Column where a DataFrame is expected and will fail at runtime.

6. **R-group analysis with a scaffold SMARTS.**
 ```typescript
 const rg = await grok.chem.rGroup(
 t, 'smiles', 'O=C1Nc2ccccc2C(C2CCCCC2)=NC1');
 grok.shell.addTableView(rg);
 ```
 Expected: a new `DataFrame` with one column per R-group position
 (`R1`, `R2`, …) holding the matched substituents per row.

7. **Compute descriptors.**
 ```typescript
 const desc = await grok.chem.descriptors(
 grok.data.testData('molecules', 100), 'smiles',
 ['MolWt', 'NumHAcceptors', 'NumHDonors', 'TPSA']);
 grok.shell.addTableView(desc);
 ```
 Expected: a `DataFrame` with the source column plus one numeric
 column per requested descriptor. Lipinski / Crippen / MolSurf /
 Fragments group names are also accepted; see
 `packages/ApiSamples/scripts/domains/chem/descriptors.js` for the
 full catalogue.

8. **Render a molecule into a `div` (sync, OpenChemLib).**
 ```typescript
 import * as ui from 'datagrok-api/ui';
 ui.dialog
.add(grok.chem.svgMol(
 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1', 300, 200))
.show;
 ```
 Expected: a dialog with the molecule rendered as inline SVG.
 `svgMol` is synchronous and returns the `HTMLDivElement`
 immediately; the SVG is filled in after a dynamic
 `import('openchemlib/full.js')` resolves (knowledge `DG-FACT-241`).

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
 Expected: the molecule drawn into the canvas with the benzene
 scaffold highlighted. Pass the scaffold as SMARTS; `canvasMol` is
 async (knowledge `DG-FACT-241`).

## Common failure modes

- **`searchSubstructure` is missing rows you expected.** Unparseable
 molecule strings are silently skipped — bit stays `0`, no throw
 (knowledge `DG-FACT-238`). Fix: pre-validate or supply
 `molBlockFailover` SMARTS. The article's `substructLibrary: true` is
 dropped by the wrapper.
- **`mcs(column)` throws / returns garbage.** Signature is
 `(DataFrame, columnName,...)`, not `(Column)`. Fix: pass the table and column name string.
- **`diversitySearch(col, METRIC_TANIMOTO, 10)` — `METRIC_TANIMOTO is
 undefined`.** The constant doesn't exist; second arg is a settings
 object. Fix: `(col, { limit: 10 })`.
- **`getSimilarities(col, null)` returns nothing useful.** Cache-priming
 expects `''`, not `null`.
- **`canvasMol` shows nothing.** It's async — `await` it before the
 dialog opens (knowledge `DG-FACT-241`).
- **`svgMol` returns an empty `div`.** OpenChemLib loads dynamically;
 SVG is injected on the next microtask. Append the div first; do not
 read `innerHTML` synchronously (knowledge `DG-FACT-241`).

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
 - `create-cell-renderer` (rendering molecules as table cells uses the
 same `svgMol`/`canvasMol` engines).
 - `access-data` (load a `DataFrame` of molecules from SQL/CSV before
 feeding it to `grok.chem.*`).
