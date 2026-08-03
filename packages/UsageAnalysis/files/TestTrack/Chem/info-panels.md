---
feature: chem
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: []
realizes: []
realized_as:
  - info-panels-spec.ts
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/info-panels.md
migration_date: 2026-05-11
source_text_fixes: []
candidate_helpers:
  - helpers.playwright.chem.openMoleculeColumnContextPanel
  - helpers.playwright.chem.openMoleculeCellContextPanel
  - helpers.playwright.chem.assertChemPanelTabsPresent
unresolved_ambiguities:
  - exact-tab-inventory-across-phase-a-step-4-and-phase-a-step-12
  - a-structure-for-the-highlight-panel-phase-a-step-8
  - block-1-vs-block-2-dataset-reuse-for-the-smiles-format
scope_reductions:
  - id: SR-01
    check: A-STRUCT-02
    rationale: "A-STRUCT-02 carryforward (chain-level edge/perf coverage)"
    verdict_status: SCOPE_REDUCTION
related_bugs: []
---

# Info Panels

End-to-end exercise of the Chem Context Panel surface. Two phases:

- **Phase A** — open `smiles.csv`, click the molecule column header, walk
  every Chemistry / Biology / Structure tab on the Context Panel. Then
  open `chembl-scaffolds.csv` to drive the **Chemistry > Rendering** panel
  (scaffold-column alignment + highlight) and the **Chemistry > Highlight**
  panel (custom substructure highlight at a distinct color, layered over
  the scaffold highlight). Finally click the first molecule cell in
  `smiles.csv` and expand every tab on the molecule Context Panel.
- **Phase B** — multi-format coverage walk. For each of four molecule
  notation formats (SMILES / molV2000 / molV3000 / SMARTS), open the
  format's reference dataset, confirm RDKit cell-rendering, click a
  structure cell, confirm the Chem Context Panel populates with the
  applicable tabs, expand each tab and verify content renders without
  errors.

## Setup

1. **Provision linked datasets.** The scenario consumes five Datagrok
   files; all are bundled in the platform's System file shares (no
   external provisioning required):
   - `System:DemoFiles/chem/smiles.csv` — primary SMILES dataset
     (canonical_smiles column, used in Phase A column-header walk and
     Phase B SMILES format slice).
   - `System:AppData/Chem/chembl-scaffolds.csv` — scaffolds dataset with
     a `Scaffold` column, drives the Rendering panel scaffold-column
     alignment + the Highlight panel layered highlight in Phase A.
   - `System:AppData/Chem/mol1K.sdf` — molV2000 format reference
     (Phase B molV2000 slice).
   - `System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf` — molV3000 format
     reference (Phase B molV3000 slice).
   - `System:AppData/UsageAnalysis/test_datasets/SMARTS_example_temp.csv`
     — SMARTS format reference (Phase B SMARTS slice).
2. **Confirm Chem package is loaded** so the Context Panel tabs
   (Chemistry / Biology / Structure) are registered. The RDKit cell
   renderer is the default Molecule renderer; package property `Renderer`
   should be RDKit (default) for this scenario.
3. **No fixture consumed.** Per chain YAML `depends_on: []`. Each
   dataset is opened fresh inside the scenario steps.

## Scenarios

### Phase A — Context Panel walk on smiles.csv + Rendering / Highlight on chembl-scaffolds

Walks the column-header Context Panel surface on `smiles.csv`, drives
the Rendering and Highlight panels on `chembl-scaffolds.csv`, then walks
the cell-level Context Panel surface on the first molecule of
`smiles.csv`.

1. Open `System:DemoFiles/chem/smiles.csv`.
2. Verify the `canonical_smiles` column auto-detects as a Molecule
   column (RDKit renderer renders cells; SMILES auto-detection populates
   units / semType / cell renderer).
3. Click the `canonical_smiles` column header to drive Context Panel
   into column-context mode.
4. In the Context Panel, walk every Chemistry / Biology / Structure
   info panel (tabs): **Chemistry > Rendering**, **Chemistry > Highlight**,
   **Chemistry > Descriptors**, **Biology > Drug Likeness**,
   **Chemistry > Properties**, **Biology > Structural Alerts**,
   **Biology > Pharmacophore Features**, **Structure > Identifiers**,
   **Structure > 2D Structure**, **Structure > 3D Structure**,
   **Biology > Toxicity**. Each tab is expandable; tab content renders
   with no console errors.
5. Open `System:AppData/Chem/chembl-scaffolds.csv` to exercise the
   Rendering panel scaffold-alignment surface.
6. With `chembl-scaffolds.csv` active and the molecule column header
   selected, open the **Chemistry > Rendering** panel in the Context
   Panel. Choose `Scaffold` as the Scaffold column. Check the
   **Highlight scaffold** checkbox.
7. Verify: molecules in the Smiles column (the molecule cells) are
   visually aligned by the chosen scaffold AND the scaffold substructure
   is highlighted in the rendered cells.
8. With the **Chemistry > Highlight** panel still in the Context Panel,
   add a structure to the Highlight list and choose a color distinct
   from the Scaffold highlight color.
9. Verify: both highlightings are visible simultaneously on the molecule
   cells in the Smiles column — the scaffold-alignment highlight (from
   step 6-7) AND the custom-structure highlight (from step 8) layered at
   different colors, no clobbering, no rendering errors.
10. Switch back to `smiles.csv`. Click the first molecule cell in the
    `canonical_smiles` column.
11. Verify: Context Panel switches from column-context to molecule-context
    (cell-level data — single-molecule info panels become available).
12. Expand every tab on the molecule Context Panel
    (**Chemistry > Descriptors**, **Biology > Drug Likeness**,
    **Chemistry > Properties**, **Biology > Structural Alerts**,
    **Biology > Pharmacophore Features**, **Structure > Identifiers**,
    **Structure > 2D Structure**, **Structure > 3D Structure**,
    **Biology > Toxicity**). Each tab populates its widget; no errors
    appear in the console or on the panel.

### Phase B — Context Panel walk across four notation formats (data-driven)

Parameterized walk across the four molecule notation formats supported by
the Chem cell renderer and Context Panel. Each row of the table below is
one iteration of the same step sequence (Steps 1-6 below) against the
named dataset.

| format    | dataset                                                                          |
|-----------|------------------------------------------------------------------------------------|
| smiles    | `System:DemoFiles/chem/smiles.csv`                                               |
| molV2000  | `System:AppData/Chem/mol1K.sdf`                                                  |
| molV3000  | `System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf`                                 |
| smarts    | `System:AppData/UsageAnalysis/test_datasets/SMARTS_example_temp.csv`             |

For each row in the table above, execute the following steps in order
on the dataset for that row:

1. Open the dataset listed for the current `format`.
2. Verify the molecule column renders — structures of molecules appear
   in cells via the RDKit cell renderer; no rendering errors. (Implicit
   format edge case: SMARTS-only cells must render via the SMARTS branch
   of the renderer; molV2000 / molV3000 cells must round-trip through
   `MolfileHandler` and render equivalently to SMILES cells.)
3. Click a cell containing a structure to drive the Context Panel into
   molecule-context mode.
4. Verify: all applicable Chem info panels appear on the Context Panel
   (Chemistry / Biology / Structure tabs registered against the molecule
   cell).
5. Expand each tab on the Context Panel.
6. Verify: the content of each info panel is displayed correctly for the
   current format — no errors, widget populated, descriptors / 2D / 3D /
   identifiers / toxicity / etc. tab content renders for the format's
   molecule.

## Notes

- **Renderer assumption.** This scenario assumes the default Chem package property `Renderer` is
  set to RDKit. If the environment has switched to OpenChemLib, the scaffold-highlight and
  substructure-highlight assertions may differ under the OCL renderer branch — that's out of
  scope here and is a separate coverage surface.
- **Tab inventory caveat.** The exact set of tabs visible in the Context Panel is
  environment-dependent — packages like Bio, Peptides, or ChemDraw may inject additional tabs not
  authored by Chem. This scenario verifies only the Chem-authored panel set, not the union of
  every tab any installed package contributes.
- **Bugs not covered here.** Two related bugs are covered by dedicated bug-focused specs, not this
  scenario: the RDKit cell renderer crashing in non-Chem viewer tooltips (GROK-16870) — this
  scenario only walks Chem Context Panel rendering, not other viewers' tooltips — and the Convert
  Notations duplicate-registration leak (GROK-17964) — this scenario never triggers that action.
