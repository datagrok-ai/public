---
feature: biostructureviewer
target_layer: playwright
coverage_type: regression
priority: p1
realizes: [biostructureviewer.viewer.biostructure, biostructureviewer.viewer.ngl]
produced_from: ticket-review
related_bugs:
  - GROK-17485
  - GROK-17967
  - GROK-16143
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as: []
---

# BiostructureViewer — Structure loaded into an empty viewer: current-row ligand and project persistence

A Biostructure viewer added to a table with no structure column shows a **Data file** input
(GROK-16143). The chosen file is kept in the viewer's `dataJson` property, so it must survive
saving and reopening the project (GROK-17485). With a ligand table (`1bdq-obs-pred.sdf`) the
viewer must show only the current row's ligand over the protein, not all ligands at once
(GROK-17967).

## Setup

- Logged in; BiostructureViewer and Chem installed.
- Files:
  - `System:AppData/BiostructureViewer/samples/1bdq-obs-pred.sdf` — 22 ligand poses with
    Observed / Predicted values.
  - `System:AppData/BiostructureViewer/samples/1bdq.pdb` — the protein.
  - A local copy of the same protein for the NGL **Open...** link: in the repository,
    `public/packages/BiostructureViewer/files/samples/1bdq.pdb`.
- Cleanup: delete the project `BsvDataFilePersistence` if it exists, before and after.

## Scenarios

### Scenario 1 — Empty Biostructure viewer loads a structure from Data file (GROK-16143)

Steps:

1. Open **Browse > Files > App Data > BiostructureViewer > samples** and double-click
   `1bdq-obs-pred.sdf`.

   * Expected result: a table view with 22 rows; the molecule column has the Molecule
     semantic type.

2. Click **Add viewer** in the toolbox, search `Biostructure` and click it.

   * Expected result: a Biostructure viewer is docked. It shows a **Data file** input
     (`.bsv-viewer-splash`) and no Mol\* engine yet. In the viewer settings (gear icon),
     **Data > Ligand Column Name** names the molecule column.

3. In **Data file**, choose `1bdq.pdb` from **App Data > BiostructureViewer > samples**.

   * Expected result: the **Data file** input disappears; the viewer holds the Mol\* engine
     (`.msp-plugin`). No error balloon.

### Scenario 2 — Only the current row's ligand is shown (GROK-17967)

Steps:

1. With the viewer from Scenario 1, click row 1 in the grid, then row 5.

   * Expected result: after each click exactly one ligand structure (the current row's) is
     loaded besides the protein — not 22. No error balloon.

2. In the viewer settings, turn **Behaviour > Show Current Row Ligand** off.

   * Expected result: no ligand structure is loaded besides the protein.

3. Turn **Show Current Row Ligand** back on. Add an **NGL** viewer (**Add viewer**, search
   `NGL`). In the NGL viewer click **Open...** and pick the local `1bdq.pdb`. Click row 3 in
   the grid.

   * Expected result: the NGL viewer shows the protein and loads exactly one ligand (row 3)
     besides it. No error balloon.

### Scenario 3 — Structure survives save and reopen of the project (GROK-17485)

Steps:

1. Close the NGL viewer. Click **SAVE** on the toolbar; in the **Save project** dialog name the
   project `BsvDataFilePersistence` and click **OK**.

   * Expected result: the project is saved; no error balloon.

2. Close all views.
3. Open **Browse > Dashboards**, find `BsvDataFilePersistence` and open it.

   * Expected result: the table view has 22 rows; the Biostructure viewer holds the Mol\*
     engine (`.msp-plugin`) and does **not** show the **Data file** input. No error balloon
     "Parsed object is empty".

## Automation notes

- Scenario 2 needs a reading of how many structures the viewer has loaded; the page shows no
  element for it. A package-level reading is required (Mol\*: the number of structures in the
  plugin's structure hierarchy; NGL: the number of components on the stage).
- The NGL **Open...** link is a native file input that takes a file from the local disk, not
  from Datagrok file shares.
- The cross-user leg of GROK-17485 (share the project, reopen as another user) is not part of
  this file; the structure travels in the layout the same way.
