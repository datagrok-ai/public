---
feature: chem
sub_features_covered: [chem.panels, chem.panels.rendering, chem.panels.highlights, chem.panels.descriptors, chem.panels.drug-likeness, chem.panels.properties, chem.panels.structural-alerts, chem.panels.pharmacophore, chem.panels.identifiers, chem.panels.structure-2d, chem.panels.structure-3d, chem.panels.toxicity, chem.rendering, chem.rendering.molecule-cell, chem.rendering.rdkit-renderer, chem.notation.detect-smiles]
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/info-panels.md
migration_date: 2026-05-11
migration_report: info-panels-migration-report.md
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

Per chain YAML (`scenario-chains/chem.yaml` rev 1): independent scenario,
`classification: medium`, `pyramid_layer: integration`, target_layer
`playwright`, strategy `simple`. UI coverage owned (`chem-info-panel-*`,
`context-panel-column-header-click`, `context-panel-cell-click`) — not
delegated. No cross-file fixtures.

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
   column (RDKit renderer renders cells; SMILES auto-detection via
   `chem.notation.detect-smiles` populated units / semType / cell
   renderer).
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
|-----------|----------------------------------------------------------------------------------|
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

- **Renderer assumption.** This scenario assumes the default Chem
  package property `Renderer` is set to RDKit. If the environment has
  switched the package property to OpenChemLib, the assertions on
  scaffold-highlight and substructure-highlight may differ in the OCL
  renderer branch — out of scope for this scenario; OCL coverage is
  parallel and lives in the `chem.rendering.ocl-renderer` test surface
  (no dedicated current test in section).
- **Tab inventory caveat.** The exact set of tabs visible in the Context
  Panel is environment-dependent — packages like Bio, Peptides, or
  ChemDraw may inject additional tabs not authored by Chem. Phase A
  step 4 / Phase B step 4 verify "all applicable Chem info panels" —
  the Chem-authored set per atlas `chem.panels.*` (12 panel ids) — not
  the universal-and-superset of every tab any installed package
  contributes.
- **No JS API substitution.** Phase A steps 3-9 and Phase B step 3
  exercise the Context Panel mode-switch (column-header click vs. cell
  click) which is itself the UI behavior under test (`context-panel-
  column-header-click`, `context-panel-cell-click` per chain
  `ui_coverage_responsibility`). UI driving is required for those
  steps; JS API direct invocation of panel computation functions would
  not exercise the mode-switch or the tab-render lifecycle.
- **Cross-cutting bug awareness — GROK-16870.** Per chain
  `bug_focused_candidates[]`, `chem-grok-16870-spec.ts` is the bug-
  focused candidate spec for the RDKit cell renderer crash on null
  DataFrame in tooltip context of non-Chem viewers (Box Plot / Scatter
  Plot / Histogram / Grid). info-panels.md exercises Chem context-panel
  rendering across Chem panels but does NOT exercise non-Chem viewer
  tooltip rendering of molecule cells — that invariant lives in the
  parallel bug-focused spec. `related_bugs: []` in the frontmatter
  because the bug-focused candidate owns the invariant; awareness only
  here.
- **Cross-cutting bug awareness — GROK-17964.** Per chain
  `bug_focused_candidates[]`, `chem-grok-17964-spec.ts` is the bug-
  focused candidate for the Convert Notations action-registration leak
  invariant. info-panels.md walks the Chem Context Panel surface but
  does NOT trigger a Convert Notations error to assert exactly-once
  registration after handler error + project reload. Awareness only;
  not added to `related_bugs` (parallel-coverage candidate owns the
  invariant).
- **Helpers usage.** No registered helper in `helpers-registry.yaml`
  currently abstracts the "open dataset → click cell → walk Chem context
  panel tabs" pattern. Candidate helper surfaced in migration report
  Decisions section; spec implementation can inline the pattern until a
  helper lands.
- **Order in chain.** `order: 1` per source JSON; tie-broken
  lexicographically before `Advanced/scaffold-tree-functions.md`. Runs
  early in section. No `must_run_last` constraint.
