---
feature: chem
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: []
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

Multi-format Activity Cliffs viewer walk: open a dataset, run **Chem | Analyze | Activity Cliffs** with
default parameters, exercise the Show-only-cliffs toggle, the cliff-count link, the cliffs-grid row click,
and the scatter-plot zoom / line-click interaction, then re-run with custom parameters. Repeated across
SMILES, molV2000, molV3000, and 2-column SMILES dataset variants — smiles.csv is the representative case,
and the remaining variants exercise the dialog's column-detection and notation-handling branches.

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
   is registered (package source `Chem/src/package.ts`). The custom editor (`ActivityCliffsEditor`)
   must surface similarity threshold + method + metric + preprocessing inputs.
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
   opens. For Variant D (two molecule columns) select the first molecule column in the column picker;
   defaults are acceptable for the remaining inputs.
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

- **Deferred.** The bounded-row confirmation dialog (which appears above 5000 / 2000 rows depending on
  method) is not exercised here — the linked datasets are all below that threshold, so this scenario
  stays a happy-path walk. The dialog-level empty-input / invalid-input failure mode is covered
  separately by sibling `chemical-space.md`, not by this scenario.
- **Sibling spec.** A Playwright spec already exists at `activity-cliffs-spec.ts` ("Chem: Activity
  Cliffs"); it is extended to cover the full dataset-variant matrix described here.
