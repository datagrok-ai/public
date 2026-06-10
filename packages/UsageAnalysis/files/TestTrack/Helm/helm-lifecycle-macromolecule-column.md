# Helm — lifecycle chain for the Macromolecule HELM column

## Setup

1. Authenticate to Datagrok as the test user.
2. Dataset: `System:AppData/Helm/samples/HELM.csv` (platform-shipped
   540-row HELM peptide sample). Opening the file triggers
   `grok.data.detectSemanticTypes`; the Bio Macromolecule detector
   classifies the `HELM` column with `semType=Macromolecule`,
   `units=helm`, `cell.renderer=helm` so the grid auto-applies
   `HelmGridCellRenderer`. The Helm package `@init` (`initHelm`) is
   expected to have run at platform startup; on first cold open the
   init chain (RDKit module load, Dojo 1.10.10 + JSDraw2 + Pistoia
   HELM Web Editor, monomer dictionary rewrite via
   `_package.completeInit`) takes a few seconds. Subsequent opens
   are warm.
3. Datagrok monomer library: default platform monomer library loaded
   by Bio is sufficient. Do NOT swap the monomer lib during the run
   — the `helm-input-bio-monomer-lib` cross-feature interaction
   covers that flow and lives outside this scenario.

## Scenarios

### Scenario 1: render_helm_cell — renderer attaches on the Macromolecule HELM column and invalidates on tag change

Steps:
1. Open `System:AppData/Helm/samples/HELM.csv`.
2. Wait for the table view to render. Confirm the `HELM` column
   shows monomer hexagons via `HelmGridCellRenderer` and the
   `Activity` column shows numeric values.
3. Scroll the grid by ~50 rows and back to row 0 to force a
   re-render of the same cells.
4. Re-apply the column tag via the API path
   (`column.setTag(DG.TAGS.CELL_RENDERER, 'helm')`) to force renderer
   cache invalidation, then scroll again.

Expected:
- The `HELM` column initially renders monomer hexagons in every
  visible row; no raw `PEPTIDE1{...}$$$$` text leaks through.
- After scrolling, re-rendered cells still show monomer hexagons —
  the `HelmGridCellRendererBack` LRU cache rehydrates without
  console errors.
- After the tag re-apply forces cache invalidation, cells re-render
  cleanly (the renderer reattaches; no orphaned bitmap state).

### Scenario 2: edit_helm_cell — cellEditor opens, edits round-trip through HELM string and commit to the cell

Steps:
1. With `HELM.csv` open from Scenario 1, double-click a non-empty
   HELM cell.
2. Confirm the full-screen HELM Web Editor dialog opens with the
   sequence drawn as connected monomer hexagons. The bottom tabs
   read `Sequence`, `HELM`, `Properties`, `Structure View`.
3. Click the bottom `HELM` tab and note the raw HELM string for the
   loaded cell.
4. Remove the trailing monomer from the raw HELM text and click
   `Apply` (this re-draws the structure from raw text without
   committing to the grid).
5. Click footer `OK`.
6. Verify the grid cell now renders the edited HELM string (one
   fewer monomer hexagon than before).
7. Cancel out of a re-opened dialog to confirm CANCEL discards
   pending edits.

Expected:
- Double-click activates `editMoleculeCell` (`role: cellEditor`,
  tag-gated on `quality=Macromolecule, units=helm`).
- `Apply` replaces the drawing from raw text but does NOT mutate
  the grid (Common Pitfall #5).
- Footer `OK` calls `cell.setValue(helm)`; the grid cell renders
  the edited sequence on dialog close.
- `CANCEL` on a re-opened dialog leaves the cell unchanged.

### Scenario 3: compute_properties — Properties panel formula / MW / extinction coefficient with the >1000-char guard

Steps:
1. Single-click a non-empty HELM cell so it becomes the current
   cell.
2. Open the Datagrok Context Pane (right side); locate the
   `Properties` panel and expand it.
3. Read the formula, molecular weight, and extinction coefficient
   from the panel's key/value table.
4. Click a different HELM cell; confirm the Properties panel
   updates to the new sequence's values.
5. Synthesize a >1000-char HELM string (e.g. by chaining
   `r(A)p.r(C)p.r(G)p.r(U)p` ~300 times into an `RNA1{…}$$$$V2.0`
   wrapper) and assign it to a single cell of a scratch
   Macromolecule HELM column via the JS API; reopen the Properties
   panel against that cell.

Expected:
- Properties panel shows three key/value rows: `formula` (e.g.
  `C101H140N23O31P`), `molecular weight`, and `extinction
  coefficient`, sourced from the off-screen JSDraw2 editor's
  `getFormula(true)` / `getMolWeight` / `getExtinctionCoefficient`.
- The panel updates whenever the current cell changes.
- For a >1000-char sequence, the panel raises a
  `SeqPropertiesError` and renders a `Too long sequence` warning
  instead of the key/value table — no UI freeze.

### Scenario 4 (bundled — convert_helm_to_molfile): batch HELM → molfile conversion across the column

Steps:
1. With `HELM.csv` open, call `Helm:getMolfiles` against the `HELM`
   column via the JS API:
   `await grok.functions.call('Helm:getMolfiles', { col: helmColumn })`.
2. Wait for the call to resolve and inspect the returned column.
3. Wait > 1 second after the call returns, then re-issue the same
   call to confirm the cached off-screen JSDraw2 editor evicts and
   re-instantiates without error.

Expected:
- The returned object is a string column with one molfile per row;
  row count equals the source column's row count.
- Each non-null entry is a valid molfile (`V2000` or `V3000`
  header; atom/bond counts > 0 for non-trivial peptides). Empty
  source rows yield empty / null entries.
- Re-issuing the call after the 1-second eviction window succeeds;
  no `JSDraw2 editor is null` console error.

### Scenario 5 (bundled — parse_helm): DOM-free HELM parser handles peptide and RNA fixtures

Steps:
1. Acquire the HELM helper:
   `const hh = await grok.functions.call('Helm:getHelmHelper')`.
2. Call `hh.parse('PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et}$$$$')`
   and inspect the returned `HelmMol` (`atoms`, `bonds`).
3. Call the lighter pure-string utility on the same input:
   `parseHelm('PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et}$$$$')`
   (imported from `public/packages/Helm/src/utils/index.ts`).
4. Repeat with an RNA fixture
   `'RNA1{r(A)p.r(C)p.r(G)p.r(U)p}$$$$'`.

Expected:
- `hh.parse(...)` returns a non-null `HelmMol`; atom count matches
  the monomer count of the source HELM string; bond count matches
  the connecting-bond count (n-1 for a linear chain of n
  monomers).
- The pure-string `parseHelm(...)` returns a flat array of monomer
  symbols whose length equals the monomer count above; the
  multi-char square-bracketed names (`[L-hArg(Et,Et)]`) are
  preserved as single tokens.
- Both fixtures parse without throwing; no DOM access (the parser
  uses `JSDraw2.MolHandler`, DOM-free).

### Scenario 6 (bundled — remove_gaps): gap-stripping transform preserves bonds across linear and cyclic positions

Steps:
1. With the HELM helper from Scenario 5, call
   `hh.removeGaps('PEPTIDE1{A.*.G.*.K}$$$$')` (gap monomer `*` at
   positions 2 and 4 of a 5-monomer linear chain).
2. Inspect the returned `{ helm, mapping }`.
3. Call the pure-string utility
   `removeGapsFromHelm('PEPTIDE1{A.*.G.*.K}$$$$')` and compare.
4. Edge boundary: call `hh.removeGaps('PEPTIDE1{*.A.G.K}$$$$')`
   (gap at start) and `hh.removeGaps('PEPTIDE1{A.G.K.*}$$$$')`
   (gap at end).

Expected:
- The first call returns a HELM string with the two `*` monomers
  removed and surrounding bonds re-linked across the gaps
  (`PEPTIDE1{A.G.K}$$$$` shape); the `mapping` object maps
  original positions {1,3,5} → {1,2,3} and omits the removed
  positions.
- The pure-string utility returns the gap-stripped HELM but does
  NOT return a position mapping (it is a pure regex strip).
- Start / end edge cases strip cleanly without throwing.
- A gap monomer with > 2 bonds (test a synthetic 3-bonded `*`)
  throws `HelmNotSupportedError`.
