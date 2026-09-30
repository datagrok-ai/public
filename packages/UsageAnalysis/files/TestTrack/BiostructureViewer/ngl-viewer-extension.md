---
feature: biostructureviewer
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [biostructureviewer.int.grid-cell-to-viewer]
realizes: [biostructureviewer.viewer.ngl]
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - ngl-viewer-extension-spec.ts
---

# BiostructureViewer — NGL viewer: opening from a Molecule3D cell, empty state, settings

The NGL viewer does not take its structure from a Molecule3D column. It gets it from the grid
cell menu **Show > NGL** (the cell value is passed to the viewer), from the table's `.pdb` tag,
or from a local file picked with its **Open...** link. Its **Ligand Column Name** accepts
Molecule (small-molecule) columns only.

NGL file previews and double-click opening are covered by
`biostructureviewer-file-open-and-preview.md`; the grid cell menu itself by
`biostructureviewer-bug-grok-14552.md`; the **PDB id viewer** panel by
`biostructureviewer-network-ui.md`.

## Setup

- Logged in; BiostructureViewer installed.
- Open **Browse > Files > App Data > BiostructureViewer** and double-click `pdb_data.csv`
  (6 rows; `pdb` is a Molecule3D column; no Molecule column).

## Scenarios

### Scenario 1 — Open NGL from a Molecule3D cell and change its settings

Steps:

1. Right-click the `pdb` cell of row 1 and choose **Show > NGL**.

   * Expected result: an NGL viewer is docked in the table view; it holds an NGL host
     (`.d4-ngl-viewer` with a `canvas`) and no **Open...** link. No error balloon.

2. Open the NGL viewer settings (gear icon on its title bar), category **Style**.

   * Expected result: **Representation** offers `cartoon`, `backbone`, `ball+stick`,
     `licorice`, `hyperball`, `surface`.

3. Set **Representation** to `ball+stick`.

   * Expected result: the settings show `ball+stick`. No error balloon, no console error.

4. In the category **Data**, open the **Ligand Column Name** choice.

   * Expected result: the `pdb` column is not offered (the property accepts Molecule columns
     only).

### Scenario 2 — NGL added to a table without a structure source shows Open...

Steps:

1. Close the NGL viewer from Scenario 1. Click **Add viewer** in the toolbox, search `NGL` and
   click it.

   * Expected result: an NGL viewer is docked showing an **Open...** link and no structure,
     because it does not read the Molecule3D column. No error balloon.

## Automation notes

- Scenario 1 asserts the viewer state only; what is drawn in the canvas (the structure in the
  chosen style, ligand overlays) is not asserted.
- The **Open...** link is a native file input for a local file; loading a file through it is
  exercised in `biostructureviewer-data-file-persistence.md`.
