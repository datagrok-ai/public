---
feature: chem
sub_features_covered: [chem.analyze.chemical-space, chem.analyze.chemical-space.top-menu, chem.analyze.chemical-space.transform, chem.analyze.chemical-space.editor, chem.analyze.chemical-space.embeddings]
target_layer: playwright
coverage_type: regression
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

Multi-format Chemical Space dialog walk: open dataset, run **Chem | Analyze | Chemical Space...** with default
parameters, then re-run with custom parameters (method, similarity metric, Cluster MCS toggle). Iterated across
SMILES, molV2000, and molV3000 dataset variants per chain rev 2 directive (Step 1 prose: "should be tested on
smiles, molV2000, molV3000 formats" — the original JSON footer listed `smiles.csv` only; extended during
migration to match body prose, sibling activity-cliffs.md / calculate.md / elemental-analysis.md precedent).

Per chain YAML (`scenario-chains/chem.yaml` rev 2): independent scenario (`depends_on: []`),
`classification: medium`, `pyramid_layer: integration`, `target_layer: playwright`, strategy `simple`.
Multi-format walk — Variant C representative source: smiles.csv; the molV2000 and molV3000 variants surface
notation-handling branches in `ChemSpaceEditor` and the dim-reduction transform. UI coverage owned
(`ui_coverage_delegated_to: null`) over the Chemical Space editor dialog, method selector, and Cluster MCS
toggle. Overlaps atlas critical path `chem.cp.chemical-space` (p1).

## Setup

1. **Provision linked datasets.** All three files are bundled in the platform's System file shares
   (no external provisioning required):
   - `System:AppData/Chem/tests/smiles-50.csv` — SMILES-notation molecule column (Variant A — canonical
     SMILES).
   - `System:AppData/Chem/mol1K.sdf` — molV2000 SDF (Variant B — V2000 molblock notation).
   - `System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf` — molV3000 SDF (Variant C — V3000 molblock
     notation).
2. **Confirm Chem package is loaded** so that the **Chem | Analyze | Chemical Space...** top-menu entry is
   registered (per atlas `chem.analyze.chemical-space.top-menu`; package source `Chem/src/package.ts#L823`).
   The custom editor (`ChemSpaceEditor`, atlas `chem.analyze.chemical-space.editor`) must surface the
   molecule column picker + method selector (UMAP / t-SNE) + similarity metric + Cluster MCS toggle.
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
   (atlas `chem.analyze.chemical-space.editor` — `ChemSpaceEditor` surfaces the molecule column picker,
   method selector, similarity metric, and Cluster MCS toggle).
3. Click **OK** to run with default parameters. Verify a 2D scatter plot is created showing the embedding,
   with two new columns tagged `chem-space-embedding-col` (x / y axes) added to the active table
   (`chemSpaceTransform` per atlas `chem.analyze.chemical-space.transform`); no console errors.
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

- **`coverage_type: regression`** — multi-format dialog + transform walk across the Chemical Space
  surface. Not `smoke` (the section's smoke is `Advanced/scaffold-tree-functions.md` per chain
  `ui_coverage_plan.smoke_scenario`); not `edge` / `perf` (no specific failure-mode invariant being
  asserted here — the empty-column-input silent-close failure mode flagged in atlas
  `chem.cp.chemical-space` is owned by bug-focused candidate `chem-grok-18407-spec.ts` per chain
  `bug_focused_candidates[]`, parallel coverage on the same `chem.analyze.chemical-space.*`
  sub_features). This scenario walks the happy path only across three notation formats.
- **No JS API substitution.** Every entry in chain `ui_coverage_responsibility`
  (`chem-add-chemical-space`, `chem-chemical-space-editor-dialog`,
  `chem-chemical-space-method-selector`, `chem-chemical-space-cluster-mcs-toggle`) is exercised via UI
  driving — top-menu walk + editor dialog interaction (column picker, method selector, Cluster MCS
  toggle) + OK click. Direct invocation of `chemSpaceTransform` or `getChemSpaceEmbeddings` is NOT a
  substitute for the dialog-driven flow that this scenario asserts.
- **Existing sibling spec.** A test file already exists at
  `public/packages/UsageAnalysis/files/TestTrack/Chem/chemical-space-spec.ts` (per
  `existing-test-index.yaml`: `test_name: "Chem: Chemical Space"`, `category: Chem`, `layer: playwright`).
  The Automator will extend it to cover the full enumerated 3-format dataset variant matrix per this
  migrated scenario.
- **Parallel bug-focused coverage.** `GROK-18407` (Chemical Space silent dialog close on empty column
  input, atlas `chem.cp.chemical-space` p1 critical-path edge case) is owned by chain
  `bug_focused_candidates[].chem-grok-18407-spec.ts` (matched_steps = [3, 5] in the original body —
  default-parameters pass and edited-parameters pass, both routed through `OK` click). That spec
  exercises the empty-column-input failure-mode invariant; this scenario exercises happy-path
  multi-format coverage. The two are parallel-coverage on the same `chem.analyze.chemical-space.*`
  sub_features. `related_bugs: []` in frontmatter — the bug is owned elsewhere (see migration report
  Decisions § Cross-cutting bug citations).
- **Helpers usage.** Standard `loginToDatagrok` + `softStep` + `closeAllViews` from
  `helpers-registry.yaml` cover the section harness. No registered helper currently abstracts the
  Chemical Space dialog walk pattern; candidate helpers surfaced in migration report Decisions.
- **Order in chain.** `order: 8` per source JSON.
- **Source-text fixes silently applied during migration.** The original JSON footer listed only
  `System:AppData/Chem/tests/smiles-50.csv`, while body Step 1 prose said "should be tested on smiles,
  molV2000, molV3000 formats". The dataset list is extended here to match body prose — three format
  variants (D1 / D2 / D3) per chain rev 2 directive footer note (c) "dataset list extension" and
  sibling activity-cliffs.md / calculate.md / elemental-analysis.md precedent. See migration report
  Decisions § Source-text fixes.
