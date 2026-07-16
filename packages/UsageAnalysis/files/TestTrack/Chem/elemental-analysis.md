---
feature: chem
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: []
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/elemental-analysis.md
migration_date: 2026-05-11
source_text_fixes:
  - step-numbering-glitch-1-2-2-3
  - dataset-list-extension
  - explicit-verification-step-added
candidate_helpers:
  - helpers.playwright.chem.openChemTopMenuItem
  - helpers.playwright.chem.toggleAllCheckboxes
  - helpers.playwright.chem.assertColumnsAppended
unresolved_ambiguities:
  - turn-all-checkboxes-on-exact-checkbox-set-on-the-elemental-analysis-dialog
  - expected-column-name-pattern-after-elemental-analysis
  - per-variant-column-count-delta
  - cell-ordering-inside-the-matrix
  - original-make-sure-the-result-is-correct-implicit-verification
scope_reductions:
  - id: SR-01
    check: A-STRUCT-02
    rationale: "A-STRUCT-02 carryforward (chain-level edge/perf coverage)"
related_bugs: []
---

# Chem | Analyze | Elemental Analysis

Multi-format happy-path walk of the **Chem | Analyze | Elemental Analysis** calculator. For each
linked dataset variant (smiles / molV2000 / molV3000), open the dataset, run the Elemental
Analysis top-menu item, toggle every per-element checkbox on in the dialog, click **OK**, and
verify that the per-element atom-count columns are appended to the active table view's grid.

## Setup

1. **Provision linked datasets.** All three files are bundled in the platform's System file shares
   (no external provisioning required):
   - `System:DemoFiles/chem/smiles.csv` — SMILES-notation molecule column (Variant A — canonical
     SMILES).
   - `System:AppData/Chem/mol1K.sdf` — molV2000 SDF (Variant B — V2000 molblock notation).
   - `System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf` — molV3000 SDF (Variant C — V3000 molblock
     notation).
2. **Confirm Chem package is loaded** so that the **Chem | Analyze | Elemental Analysis...**
   top-menu item (`elementalAnalysis`, source `public/packages/Chem/src/package.ts#L961`) is
   registered. The Charts and PowerGrid packages should also be loaded so that the optional Radar
   viewer + per-row Radar grid-cell follow-ons are available; their absence surfaces a clear
   user-visible warning but does not block the per-element atom-count column append covered by
   this scenario.
3. **No fixture consumed.** Per chain YAML `depends_on: []`. Each dataset is opened fresh inside
   the scenario.

## Scenarios

### Elemental Analysis per format

Data-driven walk: for each dataset variant (D1–D3), open the dataset, run Elemental Analysis with
all per-element checkboxes toggled on, click **OK**, and verify per-element atom-count columns
appear in the grid.

Dataset matrix:

| id | dataset path | format |
|----|--------------|--------|
| D1 | `System:DemoFiles/chem/smiles.csv` | smiles |
| D2 | `System:AppData/Chem/mol1K.sdf` | molV2000 |
| D3 | `System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf` | molV3000 |

For each dataset `D`:

1. Open the dataset `D` (close any previously open views first to start from a clean state). Wait
   for the table view to render with the molecule column populated and the RDKit cell renderer
   applied.
2. From the top menu, run **Chem > Analyze > Elemental Analysis...**. The Elemental Analysis
   dialog opens with the molecule column auto-detected and a list of per-element checkboxes.
3. Toggle every per-element checkbox in the dialog **on** (select all elements).
4. Click **OK** to run the calculator with these parameters (default behavior for any remaining
   inputs).
5. Verify: per-element atom-count columns are appended to the active table view's grid (one
   column per element that was checked, with integer atom counts per row). The grid row count is
   unchanged. No console errors fire during dialog open, OK click, or column-append.
6. Close the active view before moving to the next dataset — Elemental Analysis appends columns
   to the active table, so a fresh open keeps each variant isolated.

Implicit cross-variant invariants:

- Elemental Analysis is invoked independently per dataset — no inter-variant state coupling.
- Across the 3 variants, dialog open + OK click never crashes the console; if the calculator
  cannot parse a notation variant the dialog surfaces a user-visible error balloon, not a silent
  failure (assert per-variant that EITHER per-element columns are appended OR a user-visible
  error balloon is shown — never silent).

## Notes

- **Optional Radar viewer / per-row Radar grid cells.** Charts (Radar viewer) and PowerGrid
  (per-row Radar grid cells) are optional add-on surfaces enabled when those packages are
  present. The per-element atom-count column append (step 5) is the canonical verification here —
  the optional Radar surfaces are not exercised by this scenario.
- **Sibling spec.** A Playwright spec already exists at `elemental-analysis-spec.ts` ("Chem:
  Elemental Analysis"); it likely covers a single format today and is extended to the full
  3-format walk described here.
