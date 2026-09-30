---
feature: biostructureviewer
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: [biostructureviewer.cp.viewer-add-and-render-pdb]
realizes: [biostructureviewer.viewer.biostructure]
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/BiostructureViewer/biostructure-viewer.md
migration_date: '2026-06-04'
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
related_bugs:
  - GROK-11759
realized_as:
  - biostructure-viewer-spec.ts
---

# BiostructureViewer — Smoke: Biostructure viewer on a structure table

The Biostructure viewer (Mol\* engine) added to a table that has a structure column picks the
column by itself, shows the current row's structure, switches rendering styles from its
settings and exports the structure from its context menu. Everything here runs on the
package's own data, without outside network services.

Opening structure files from the Files browser is covered by
`biostructureviewer-file-open-and-preview.md`; the overlay buttons of the Mol\* viewport by
`molstar-overlay-extension.md`; structures loaded by PDB ID from RCSB by
`biostructureviewer-network-ui.md`.

## Setup

- Logged in; the BiostructureViewer package is installed.
- File: `System:AppData/BiostructureViewer/pdb_data.csv` — 6 rows; `pdb_id` holds PDB IDs
  (detected as PDB_ID), `pdb` holds full PDB texts (detected as Molecule3D). Row 1 is 1QBS
  (HIV-1 protease with an inhibitor).
- Close all views before the scenarios.

## Scenarios

### Scenario 1 — Add the viewer to a structure table; it wires the columns by itself

Steps:

1. Open **Browse > Files > App Data > BiostructureViewer** and double-click `pdb_data.csv`.

   * Expected result: a table view with 6 rows opens.

2. Click **Add viewer** in the toolbox, type `Biostructure` in the search box and click
   **Biostructure**.

   * Expected result: a Biostructure viewer is docked. It holds the Mol\* engine (a
     `.msp-plugin` element with a `.msp-viewport` inside the viewer) and shows no
     **Data file** input. No error balloon.

3. Open the viewer settings (gear icon on the viewer title bar).

   * Expected result: in the **Data** category, **Biostructure Id Column Name** is `pdb`;
     in the **Style** category, **Representation** is `cartoon`.

4. Make row 3 current by clicking it in the grid.

   * Expected result: no error balloon; the viewer still holds the Mol\* engine.

### Scenario 2 — Switch the representation (GROK-11759)

Steps:

1. With the viewer from Scenario 1 and its settings open, change **Style > Representation**
   to `ball-and-stick`.

   * Expected result: the settings show `ball-and-stick`. No error balloon, no console error.

2. Change **Representation** to `molecular-surface`.

   * Expected result: the settings show `molecular-surface`. No error balloon, no console
     error.

3. Change **Representation** back to `cartoon`.

   * Expected result: the settings show `cartoon`. No error balloon, no console error.

### Scenario 3 — Reset Camera button

Steps:

1. In the Mol\* viewport of the viewer, click the overlay button with the tooltip
   **Reset Camera** (`title="Reset Camera"`).

   * Expected result: the button exists; after the click no error balloon and no console
     error; the viewer still holds the Mol\* engine.

### Scenario 4 — Download the structure from the viewer menu

Steps:

1. Right-click inside the Biostructure viewer.

   * Expected result: a context menu with a group **Download** holding **As CIF** and
     **As PDB**.

2. Choose **Download > As PDB**.

   * Expected result: the browser downloads `file.pdb`. No error balloon.

3. Right-click inside the viewer again and choose **Download > As CIF**.

   * Expected result: the browser downloads `file.cif`. No error balloon.

## Automation notes

- Property labels are written as the settings panel is expected to show them (camel-case
  names split into words: `biostructureIdColumnName` → **Biostructure Id Column Name**); match
  them once on the stand.
- The representation choices come from the Mol\* built-in representation list
  (`cartoon`, `ball-and-stick`, `molecular-surface`, …); the exact list shown in the
  dropdown is to be matched on the stand.
- Scenario 1 step 4 and Scenario 3 check only that the viewer keeps working; what is drawn in
  the viewport is not asserted.
- The viewer also sets **Ligand Column Name** to `pdb` (the same Molecule3D column it takes the
  structure from); this is not asserted here.

---
{
  "order": 1,
  "datasets": [
    "System:AppData/BiostructureViewer/pdb_data.csv"
  ]
}
