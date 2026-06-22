# Helm — Bio menu cross-feature integration on a HELM column

> Cross-feature smoke: walk the **Bio** top menu on a Macromolecule `units=helm`
> column and execute every leaf that is applicable to a standalone HELM column,
> asserting that Bio's consumption of the column does NOT break the Helm
> renderer / service / editor (the 2026-06 SVG-editor rewrite). The Bio-menu
> leaves are owned by sibling packages (Bio, Dendrogram, SequenceTranslator,
> BiostructureViewer) — see grok-browser `references/bio.md` § "Menu ownership
> map"; this scenario only validates that they coexist cleanly with Helm.

## Setup

1. Authenticate to Datagrok as the test user.
2. Dataset: `System.AppData/Helm/samples/helm-showcase.csv` (55-row showcase).
   **Open it directly via its instance-derived file URL** —
   `<instance>/file/System.AppData/Helm/samples/helm-showcase.csv?browse=files`
   (e.g. on dev: `https://dev.datagrok.ai/file/...`) — single navigation, no
   open-platform-then-readCsv. The Bio Macromolecule detector tags the `HELM`
   column (`semType=Macromolecule`, `units=helm`, `cell.renderer=helm`,
   `quality=Macromolecule`).
3. Make the `HELM` column the current column so the `Bio` top menu appears.
4. Drive each Bio leaf via the documented dispatchEvent vector (per
   `references/bio.md` § "Top-menu Bio entries"): click `[name="div-Bio"]` →
   hover the group `[name="div-Bio---<Group>"]` → click the leaf
   `[name="div-Bio---<Group>---<Leaf>"]`. For dialogs, accept the defaults and
   click OK (only when OK is enabled).

## Scenarios

### Scenario 1: execute every HELM-applicable Bio leaf; Helm integration survives

For each applicable leaf below: open it; if a dialog appears and OK is enabled,
click OK with default inputs; otherwise observe the docked viewer / view; then
clean up (cancel dialogs, close added viewers, remove added columns).

**Applicable leaves (executed):**

| Bio leaf | Owner | Expected on a HELM column |
|---|---|---|
| `Analyze \| Composition` | Bio | Docks a **WebLogo** viewer over the HELM column |
| `Analyze \| Sequence Space...` | Bio | Runs (UMAP); adds embedding columns. *(DR convergence warning on the small showcase is tolerated — not a Helm error.)* |
| `Analyze \| Activity Cliffs...` | Bio | Runs; adds embedding columns. *(DR warning tolerated.)* |
| `Analyze \| MSA...` | Bio | Dialog opens; OK opens an **MSA result table** (with an `msa(HELM)` aligned column) |
| `Analyze \| Hierarchical Clustering...` | Dendrogram | Dialog opens; OK runs clustering on HELM-derived distances |
| `Transform \| Convert Sequence Notation...` | Bio | **HELM consumer** — adds a converted-notation column |
| `Transform \| Split to Monomers...` | Bio | **HELM consumer** — adds per-position monomer columns |
| `Calculate \| Extract Region...` | Bio | Dialog opens; OK adds a region column |
| `Annotate \| Apply Numbering Scheme...` | Bio | Dialog opens; OK runs |
| `Annotate \| Scan Liabilities...` | Bio | Dialog opens; OK adds a liabilities column |
| `Annotate \| Manage Annotations...` | Bio | Dialog opens; OK runs |
| `Search \| Similarity Search` | Bio | Docks a **Sequence Similarity Search** viewer |
| `Search \| Diversity Search` | Bio | Docks a **Sequence Diversity Search** viewer |
| `Search \| Subsequence Search ...` | Bio | Docks the **Filters** panel (HELM substructure filter) |
| `PolyTool \| Convert...` | SequenceTranslator | Opens the conversion dialog (routes to To-Atomic-Level on this column) |
| `PolyTool \| Enumerate HELM...` | SequenceTranslator | Dialog opens; OK runs HELM enumeration |

Expected:
- The test **waits for and asserts each leaf's concrete result** (does not merely
  fire the function): viewer-docking leaves (Composition, Similarity / Diversity /
  Subsequence Search) MUST dock their viewer; column-adding leaves (Sequence Space,
  Activity Cliffs, Convert Notation, Split to Monomers, Extract Region, Scan
  Liabilities) MUST grow the column count; MSA MUST open a result table. The
  remaining leaves (Hierarchical Clustering, Apply Numbering Scheme, Manage
  Annotations, PolyTool Convert / Enumerate) run cleanly but produce no
  deterministic column/viewer artifact with defaults on this dataset, so only the
  no-error guarantee is asserted for them.
- Every applicable leaf opens / runs **without a Helm-related error balloon**
  (no error mentioning `helm` / `monomer render` / `JSDraw` / `pseudo-molfile` /
  `getMolfiles` / cell renderer). Generic Bio compute artifacts (e.g. UMAP
  "Dimensionality reduction failed" on the tiny showcase) are tolerated and
  logged, not failed.
- After the whole sweep, the `HELM` column is **still intact**:
  `semType=Macromolecule`, `cell.renderer=helm`, and the grid still renders the
  HELM structures — proving Bio's operations did not corrupt the Helm
  renderer / service.

### Skipped leaves (NOT applicable to a standalone HELM column)

Documented so the sweep's omissions are intentional, not gaps:

| Bio leaf | Why skipped |
|---|---|
| `Folding \| EsmFold...`, `Folding \| Boltz...` | Protein-structure prediction (Docker), not a HELM-column operation. |
| `Analyze \| SAR...` | Peptides SAR requires a numeric **activity** column; the showcase has none (no dialog opens). |
| `Analyze \| Compare sequences...` | Requires **≥2** Macromolecule columns; the showcase has one. Guards: "needs at least two Macromolecule columns". |
| `Transform \| To Atomic Level...` | A HELM operation, but covered separately (verified independently); excluded to avoid duplication. |
| `Transform \| Molecules to HELM...` | Requires a small-**molecule** column as input, not a HELM column. |
| `Transform \| Fetch PDB Sequences...` | Requires a PDB id (BiostructureViewer); not a HELM-column operation. |
| `Manage \| Match with Monomer Library...` | Requires a small-molecule column; on a HELM-only table OK throws an input NPE. |
| `Calculate \| Identity...`, `Calculate \| Similarity...` | Applicable to macromolecules but require a **reference sequence**; OK stays disabled with defaults, so nothing executes. |
| `PolyTool \| Combine Sequences...` | Requires ≥2 sequence inputs / placeholders; guards "Please fill all the fields" with defaults. |
| `Manage \| Monomer Libraries`, `Manage \| Monomers` | Open library-management **Views** (navigate away), not HELM-column operations. |
