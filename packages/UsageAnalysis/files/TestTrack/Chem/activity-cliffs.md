---
feature: chem
sub_features_covered: [chem.analyze.activity-cliffs, chem.analyze.activity-cliffs.top-menu, chem.analyze.activity-cliffs.transform, chem.analyze.activity-cliffs.init, chem.analyze.activity-cliffs.editor]
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/activity-cliffs.md
migration_date: 2026-05-11
source_text_fixes:
  - step-numbering-gap-closed
  - smiles-2-columns-csv-todo-resolved
  - spgi-dataset-reconciled
candidate_helpers:
  - helpers.playwright.chem.addActivityCliffsViewer
  - helpers.playwright.chem.toggleActivityCliffsShowOnlyCliffs
  - helpers.playwright.chem.clickActivityCliffsCliffCountLink
  - helpers.playwright.chem.clickActivityCliffsGridRow
unresolved_ambiguities:
  - click-on-the-link-with-number-of-cliffs-exact-dom-selector-and-label-format
  - first-row-of-the-grid-first-by-dom-render-order-or-first-by-cliff-magnitude-pair
  - hovering-over-line-line-locator-in-the-zoomed-scatter-plot
  - change-the-parameters-arbitrarily-which-specific-parameter-to-change
  - double-click-to-unzoom-exact-zoom-reset-target
  - cliffs-grid-docking-position-on-the-bottom
scope_reductions:
  - id: SR-01
    check: A-STRUCT-02
    rationale: "A-STRUCT-02 carryforward (chain-level edge/perf coverage)"
    verdict_status: SCOPE_REDUCTION
related_bugs: []
---

# Chem | Analyze | Activity Cliffs walk

Multi-format Activity Cliffs viewer walk: open dataset, run **Chem | Analyze | Activity Cliffs** with default
parameters, exercise the Show-only-cliffs toggle, the cliff-count link, the cliffs-grid row click, the
scatter-plot zoom / line-click interaction, then re-run with custom parameters. Iterated across SMILES,
molV2000, molV3000, and 2-column SMILES dataset variants per chain rev 2 directive.

Per chain YAML (`scenario-chains/chem.yaml` rev 2): independent scenario (`depends_on: []`),
`classification: medium`, `pyramid_layer: integration`, `target_layer: playwright`, strategy `simple`.
Multi-format walk — Variant C representative source: smiles.csv; the remaining variants surface dialog +
column-detection branches and notation handling. UI coverage owned (`ui_coverage_delegated_to: null`) over
the Activity Cliffs editor dialog, Show-only-cliffs toggle, cliff-count link, cliffs-grid row click,
scatter-plot zoom + line-click, and property-panel pair-of-molecules display.

## Setup

1. **Provision linked datasets.** All five files are bundled in the platform's System file shares
   (no external provisioning required):
   - `System:AppData/Chem/tests/smiles-50.csv` — SMILES-notation molecule column (Variant A — canonical
     SMILES).
   - `System:AppData/Chem/mol1K.sdf` — molV2000 SDF (Variant B — V2000 molblock notation).
   - `System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf` — molV3000 SDF (Variant C — V3000 molblock
     notation).
   - `System:AppData/Chem/tests/smiles_2_columns.csv` — 2-column SMILES dataset (Variant D — two
     molecule columns; exercises column-selection branches in the Activity Cliffs editor dialog).
   - `System:AppData/Chem/tests/spgi-100.csv` — SMILES + numeric-activity reference dataset
     (Variant E — JSON footer in the original includes SPGI alongside the four format-variant files;
     kept as completeness reference for the editor dialog's molecule + activity column auto-detect).
2. **Confirm Chem package is loaded** so that the **Chem | Analyze | Activity Cliffs...** top-menu entry
   is registered (per atlas `chem.analyze.activity-cliffs.top-menu`; package source
   `Chem/src/package.ts`). The custom editor (`ActivityCliffsEditor`,
   atlas `chem.analyze.activity-cliffs.editor`) must surface similarity threshold + method + metric +
   preprocessing inputs.
3. **No fixture consumed.** Per chain YAML `depends_on: []`. Each dataset is opened fresh inside the
   scenario.

## Scenarios

### Activity Cliffs viewer walk per dataset

Data-driven walk: for each dataset variant (D1–D5), open the dataset, run **Chem | Analyze | Activity
Cliffs** with default parameters, exercise the scatter + cliff-grid + property-panel interactions, then
re-run with custom parameters.

Dataset matrix:

| id | dataset path | format |
|----|--------------|--------|
| D1 | `System:AppData/Chem/tests/smiles-50.csv` | smiles |
| D2 | `System:AppData/Chem/mol1K.sdf` | molV2000 |
| D3 | `System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf` | molV3000 |
| D4 | `System:AppData/Chem/tests/smiles_2_columns.csv` | smiles, 2-column |
| D5 | `System:AppData/Chem/tests/spgi-100.csv` | smiles + numeric activity (reference) |

For each dataset `D`:

1. Open the dataset `D` (close any previously open views first to start from a clean state). Wait for
   the table view to render with the molecule column populated and the RDKit cell renderer applied.
2. From the top menu, run **Chem | Analyze | Activity Cliffs...**. The Activity Cliffs editor dialog
   opens (atlas `chem.analyze.activity-cliffs.editor`). For Variant D (two molecule columns) select the
   first molecule column in the column picker; defaults are acceptable for the remaining inputs.
3. Click **OK** to run with default parameters. Verify the scatter plot is created with the cliffs
   overlay drawn (`activityCliffsInitFunction` runs without console errors; the
   `activityCliffsParams` tag is set on the table).
4. Toggle **Show only cliffs** in the scatter-plot viewer's property panel / context-menu. Verify only
   points that participate in a cliff pair remain visible on the scatter plot (non-cliff points are
   hidden).
5. Click the **cliffs link** in the viewer (the link showing the cliff count, e.g. *N cliffs*). Verify a
   cliffs grid appears docked below the scatter plot and is populated with the cliff pair rows.
6. Click the **first row** of the cliffs grid. Verify the scatter plot zooms to the corresponding cliff
   region, the property panel updates to show the pair of molecules (with uncommon parts highlighted)
   and the activity difference for the selected pair.
7. Hover over the cliff line in the zoomed scatter plot. Verify a tooltip appears containing the same
   pair-of-molecules + activity-difference information as the property panel.
8. Double-click the scatter plot to unzoom. Zoom into any cluster until at least one cliff line is
   visible, then click the line. Verify the scatter plot zooms to that cliff and the corresponding row
   becomes current (highlighted) in the cliffs grid below.
9. Run **Chem | Analyze | Activity Cliffs...** a second time on the same active table.
10. In the editor dialog, change at least one parameter from its default value (e.g. similarity
    threshold, method, metric, or preprocessing option — any change that makes the run materially
    different from the default-parameters pass at step 3).
11. Click **OK** to run with edited parameters. Verify a fresh scatter plot + cliffs overlay is produced
    that reflects the edited parameters (cliff count and / or cliff distribution differs from step 3),
    with no console errors.
12. Close the active view before moving to the next dataset variant — viewer state and per-table tags
    do not need to leak between cells.

Implicit cross-cell invariants:

- Each dataset variant is exercised independently — no inter-cell state coupling. The cliff scatter +
  grid + property panel interactions all complete without console errors regardless of the input
  notation format.
- The cliffs link, cliffs grid, scatter-plot zoom interactions, and property-panel pair display all
  remain consistent across the four format variants (D1–D4) and the SPGI reference (D5).

## Notes

- **`coverage_type: regression`** — multi-format × interaction-walk across the Activity Cliffs viewer
  surface. Not `smoke` (the section's smoke is `Advanced/scaffold-tree-functions.md` per chain
  `ui_coverage_plan.smoke_scenario`); not `edge` / `perf` (no specific failure-mode invariant or
  threshold being asserted — the row-count cap (10,000 rows, atlas
  `chem.analyze.activity-cliffs.top-menu`) and the long-running-confirm dialog above 5000 / 2000 rows
  per method are not exercised by these datasets — they sit below the threshold). The dialog-level
  empty-input / invalid-input failure mode is owned by sibling `chemical-space.md` + bug-focused spec
  `chem-grok-18407-spec.ts` (parallel-coverage, different sub_feature scope — `chem.analyze.chemical-
  space.*`, not `chem.analyze.activity-cliffs.*`).
- **No JS API substitution.** Every entry in chain `ui_coverage_responsibility`
  (`chem-add-activity-cliffs`, `chem-activity-cliffs-editor-dialog`,
  `chem-activity-cliffs-show-only-cliffs-toggle`, `chem-activity-cliffs-cliff-count-link`,
  `chem-activity-cliffs-grid-row-click`, `chem-activity-cliffs-scatter-zoom-line-click`,
  `chem-activity-cliffs-property-panel-pair-molecules`) is exercised via UI driving — top-menu walk +
  editor dialog interaction + scatter-plot click / double-click / zoom + cliffs-grid row click + line
  hover. The bounded-row dialog confirm (above 5000 / 2000 rows depending on method) is NOT exercised
  because the linked datasets sit below the threshold; the cap invariant is observable only on a
  large-volume dataset and is intentionally out of scope for this regression walk.
- **Existing sibling spec.** A test file already exists at
  `public/packages/UsageAnalysis/files/TestTrack/Chem/activity-cliffs-spec.ts` (per
  `existing-test-index.yaml` line 32321: `test_name: "Chem: Activity Cliffs"`, `category: Chem`,
  `layer: playwright`, helpers `spec-login`, `features_covered: [chem.activity-cliffs]`). The Automator
  will extend it to cover the full enumerated dataset variant matrix per this migrated scenario.
- **Helpers usage.** Standard `loginToDatagrok` + `softStep` + `closeAllViews` from
  `helpers-registry.yaml` cover the section harness. No registered helper currently abstracts the
  Activity Cliffs viewer add / dialog / scatter interaction pattern; candidate helpers surfaced in
  migration report Decisions.
- **Order in chain.** `order: 7` per source JSON.
- **Source-text fixes silently applied during migration.** The original body had a step-numbering gap
  (Step 7 → Step 9, with no Step 8) and a TODO note for `smiles_2_columns.csv` ("add to linked
  datasets") that has been resolved — the file exists at
  `System:AppData/Chem/tests/smiles_2_columns.csv` (verified per repo search at
  `public/packages/Chem/files/tests/smiles_2_columns.csv`, per the calculate.md migration finding); the
  original JSON footer's
  `System:AppData/UsageAnalysis/test_datasets/smiles_2_columns.csv` path is corrected here to the
  canonical Chem package path. Steps renumbered cleanly into the migrated body's 12-step linear
  sequence. See migration report Decisions § Source-text fixes.
