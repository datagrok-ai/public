---
feature: chem
target_layer: playwright
coverage_type: regression
priority: p0
realizes: [chem.cp.chemical-space]
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/chemical-space.md
migration_date: 2026-05-11
source_text_fixes:
  - dataset-list-extended-from-1-entry-to-3-entries
candidate_helpers:
  - helpers.playwright.chem.openChemicalSpaceDialog
  - helpers.playwright.chem.changeChemicalSpaceMethod
  - helpers.playwright.chem.toggleClusterMCS
unresolved_ambiguities:
  - change-the-parameters-arbitrarily-which-specific-parameter-to-change
  - method-selector-default-value
  - fresh-scatter-plot-overwrite-vs-append-semantics
  - cluster-mcs-toggle-visibility-on-the-editor
scope_reductions:
  - id: SR-01
    check: A-STRUCT-02
    rationale: "A-STRUCT-02 carryforward (chain-level edge/perf coverage)"
    verdict_status: SCOPE_REDUCTION
related_bugs: []
---

# Chem | Analyze | Chemical Space walk

Multi-format Chemical Space dialog walk: open a dataset, run **Chem | Analyze | Chemical Space...**
with default parameters, then re-run with custom parameters (method, similarity metric, Cluster MCS
toggle). Repeated across SMILES, molV2000, and molV3000 dataset variants — `smiles.csv` is the
representative case, and the molV2000 / molV3000 variants exercise notation-handling branches in the
Chemical Space editor and the dimensionality-reduction transform.

## Setup

1. **Provision linked datasets.** All three files are bundled in the platform's System file shares
   (no external provisioning required):
   - `System:AppData/Chem/tests/smiles-50.csv` — SMILES-notation molecule column (Variant A — canonical
     SMILES).
   - `System:AppData/Chem/mol1K.sdf` — molV2000 SDF (Variant B — V2000 molblock notation).
   - `System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf` — molV3000 SDF (Variant C — V3000 molblock
     notation).
2. **Confirm Chem package is loaded** so that the **Chem | Analyze | Chemical Space...** top-menu entry is
   registered (package source `Chem/src/package.ts#L823`). The custom editor (`ChemSpaceEditor`) must
   surface the molecule column picker + method selector (UMAP / t-SNE) + similarity metric + Cluster MCS
   toggle.
3. **No fixture consumed.** Per chain YAML `depends_on: []`. Each dataset is opened fresh inside the
   scenario.

## Scenarios

### Chemical Space dialog walk per dataset

Data-driven walk: for each dataset variant (D1–D3), open the dataset, run **Chem | Analyze | Chemical
Space...** with default parameters, verify the embedding scatter plot is produced, then re-run with edited
parameters.

Dataset matrix:

| id | dataset path                                       | format    |
|----|----------------------------------------------------|-----------|
| D1 | `System:AppData/Chem/tests/smiles-50.csv`                 | smiles    |
| D2 | `System:AppData/Chem/mol1K.sdf`                    | molV2000  |
| D3 | `System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf`  | molV3000  |

For each dataset `D`:

1. Open the dataset `D` (close any previously open views first to start from a clean state). Wait for the
   table view to render with the molecule column populated and the RDKit cell renderer applied.
2. From the top menu, run **Chem | Analyze | Chemical Space...**. The Chemical Space editor dialog opens
   (`ChemSpaceEditor` surfaces the molecule column picker, method selector, similarity metric, and Cluster
   MCS toggle).
3. Click **OK** to run with default parameters. Verify a 2D scatter plot is created showing the embedding,
   with two new columns tagged `chem-space-embedding-col` (x / y axes) added to the active table
   (`chemSpaceTransform`); no console errors.
4. Run **Chem | Analyze | Chemical Space...** a second time on the same active table.
5. In the editor dialog, change at least one parameter from its default value (e.g. switch the method
   selector between UMAP and t-SNE, change the similarity metric, or toggle Cluster MCS on / off — any
   change that makes the run materially different from the default-parameters pass at step 3).
6. Click **OK** to run with edited parameters. Verify a fresh embedding scatter plot is produced that
   reflects the edited parameters (point distribution differs from step 3, or — if Cluster MCS toggled
   on — an additional `chem-space-cluster-col` column appears on the table), with no console errors.
7. Close the active view before moving to the next dataset variant — viewer state and per-table tags do
   not need to leak between cells.

Implicit cross-cell invariants:

- Each dataset variant is exercised independently — no inter-cell state coupling. The Chemical Space
  dialog + embedding scatter all complete without console errors regardless of the input notation format.
- The dialog shape (column picker + method selector + Cluster MCS toggle) is consistent across the three
  notation formats; only the auto-detected molecule column source varies.

## Notes

- **Sibling spec.** A Playwright spec already exists at `chemical-space-spec.ts` ("Chem: Chemical
  Space"); it is extended to cover the full 3-format dataset variant matrix described here.
- **Parallel bug-focused coverage.** GROK-18407 (Chemical Space dialog silently closing on empty
  column input) has no dedicated bug-focused spec yet — this scenario walks
  the happy path only, across the three notation formats.
