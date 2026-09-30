---
feature: biostructureviewer
target_layer: playwright
coverage_type: regression
priority: p1
realizes: [biostructureviewer.import.pdb, biostructureviewer.import.pdbqt, biostructureviewer.import.xyz, biostructureviewer.import.with-ngl, biostructureviewer.preview.biostructure, biostructureviewer.preview.ngl-structure, biostructureviewer.preview.ngl-surface, biostructureviewer.preview.ngl-density]
produced_from: ticket-review
related_bugs:
  - GROK-14442
  - GROK-16968
  - GROK-17654
  - GROK-18999
  - GROK-13650
  - CLAUDE-33
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as: []
---

# BiostructureViewer — Opening and previewing structure files from the Files browser

Double-clicking a structure file runs the package file handler chosen by the extension;
single-clicking shows a preview. `.pdb`, `.mmcif` and `.xyz` open in a separate view titled
**Mol\*** that holds the Mol\* engine (no table, no grid). `.pdbqt` with docking poses opens
as a table of poses. NGL-only formats preview with the NGL engine. A `.pdb` file also has the
context menu item **Open table residues**.

Regression guard for GROK-14442 (`.pdb` was routed to the `.pdbqt` importer), GROK-16968
(structure files stopped opening), GROK-17654 (PDB/CIF preview broken), GROK-18999 (empty
preview name for `.pdb` / `.pdbqt`), GROK-13650 (NGL formats failed to preview) and
CLAUDE-33 (closing an unrelated view while a structure preview was open raised
`Cannot read properties of undefined (reading 'children')`).

## Setup

- Logged in; BiostructureViewer installed.
- Files used:
  - `System:AppData/BiostructureViewer/samples/1bdq.pdb`, `1RQ9.mmcif`,
    `1rq9-assembly1.cif`, `caffeine.xyz`, `1bdq.autodock-gpu.pdbqt` (two docking poses, no
    receptor).
  - `System:DemoFiles/bio/ngl-formats/1blu.mmtf`, `1crn.ply`, `1crn.obj`, `1lee.ccp4`,
    `3pqr.cns`.
  - `System:DemoFiles/demog.csv`.
- Close all views before each scenario.

## Scenarios

### Scenario 1 — Preview PDB, mmCIF, CIF and XYZ files (GROK-17654, GROK-18999)

Steps:

1. Open **Browse > Files** and go to **App Data > BiostructureViewer > samples**.
2. Click `1bdq.pdb` once.

   * Expected result: the preview shows the Mol\* engine (`.msp-plugin` inside the preview)
     and is titled `1bdq.pdb` (not empty). No error balloon starting with "Preview file".
     No console error.

3. Click `1RQ9.mmcif` once, then `1rq9-assembly1.cif` once, then `caffeine.xyz` once.

   * Expected result: each preview shows the Mol\* engine and is titled with the file name.
     No error balloon.

4. Click `1bdq.autodock-gpu.pdbqt` once.

   * Expected result: the preview is titled `1bdq.autodock-gpu.pdbqt`. No error balloon.

### Scenario 2 — Preview NGL-only formats (GROK-13650)

Steps:

1. In **Browse > Files** go to **Demo > bio > ngl-formats**.
2. Click each of `1blu.mmtf`, `1crn.ply`, `1crn.obj`, `1lee.ccp4`, `3pqr.cns` once, in turn.

   * Expected result: for every file the preview holds an NGL host (`.d4-ngl-viewer` with a
     `canvas`), not the Mol\* engine. No error balloon starting with "Preview file". No
     console error containing `ext '' unknown`.

### Scenario 3 — Double-click a PDB opens it in a Mol* view (GROK-14442, GROK-16968)

Steps:

1. In **App Data > BiostructureViewer > samples**, double-click `1bdq.pdb`.

   * Expected result: a new view **Mol\*** becomes current and holds the Mol\* engine
     (`.msp-plugin`, `.msp-viewport`). The view has no grid. No **Open file** dialog appears
     (that dialog belongs to the `.pdbqt` importer). No error balloon.

2. Close the view. Double-click `1RQ9.mmcif`; after closing that view, double-click
   `caffeine.xyz`.

   * Expected result: each time a **Mol\*** view with the Mol\* engine opens. No error
     balloon.

### Scenario 4 — Double-click a PDBQT opens the poses as a table (GROK-14442, reverse direction)

Steps:

1. Double-click `1bdq.autodock-gpu.pdbqt`.

   * Expected result: a table view opens with **2** rows and a column named `molecule`. A
     dialog **Open file** appears with the text "Docking target structure required to display
     ligand poses from pdbqt data." No **Mol\*** view is opened.

2. Click **CANCEL** in the dialog.

   * Expected result: the dialog closes; the table view stays with 2 rows. No error balloon.

### Scenario 5 — Closing an unrelated view while a structure preview is open (CLAUDE-33)

Steps:

1. In **App Data > BiostructureViewer > samples**, click `1bdq.pdb` once so its preview is
   shown.
2. Open `System:DemoFiles/demog.csv` (**Demo > demog.csv**, double-click), then close its
   view with the close button on its tab.

   * Expected result: no error balloon; no console error containing `reading 'children'`.

3. Go back to **Browse > Files**, click `1bdq.pdb` once again, open **Help** from the
   sidebar and close it.

   * Expected result: same as step 2; the `1bdq.pdb` preview can still be opened.

### Scenario 6 — Double-click NGL-only formats opens them in an NGL view

Steps:

1. In **Demo > bio > ngl-formats**, double-click `1blu.mmtf`.

   * Expected result: a new view **NGL** becomes current and shows the structure in an NGL
     host (`.d4-ngl-viewer` with a `canvas`). No error balloon; no console error containing
     `ext '' unknown`.

2. Close the view and double-click `1lee.ccp4`.

   * Expected result: same as step 1 (density map in an NGL host).

### Scenario 7 — Open table residues from a PDB file's context menu

Steps:

1. In **App Data > BiostructureViewer > samples**, right-click `1bdq.pdb`.

   * Expected result: the context menu has the item **Open table residues**.

2. Click **Open table residues**.

   * Expected result: a table view opens with the columns `code`, `compId`, `seqId`,
     `label`, `seq`, `frame`, one row per residue; an NGL viewer is docked on the right. The
     `compId` column holds residue names (the first rows are `PRO`, `GLN`, `ILE`, `THR`,
     `LEU` — the start of chain A in `1bdq.pdb`). No error balloon.

## Automation notes

- `.pdbqt` opening is known to call the handler twice (GROK-14438, won't fix): assert "a
  table view with 2 rows exists", not "exactly one view".
- Scenario 1: where the Files browser shows the preview's name is to be matched on the stand;
  the NGL previews of Scenario 2 are created without a name, so their title is not asserted.
- Scenario 6 is a suspected defect, verify on the stand: by code reading, the double-click
  handler for NGL-only formats (`viewNglUI`) creates the **NGL** view but does not put the NGL
  host into it, and loads the file without its extension, so an empty **NGL** view or an
  `ext '' unknown` error is likely instead of the expected result.
- Scenario 7 is a suspected defect, verify on the stand: by code reading, `compId` is created
  as an integer column while three-letter residue names are written into it, so the column
  may come out empty.
