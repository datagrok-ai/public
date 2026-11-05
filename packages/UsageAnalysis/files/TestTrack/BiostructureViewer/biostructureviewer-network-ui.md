---
feature: biostructureviewer
target_layer: manual-only
coverage_type: regression
priority: p2
realizes_atlas: [biostructureviewer.cp.pdb-id-data-provider-roundtrip, biostructureviewer.int.prolif-cross-cell-type]
realizes: [biostructureviewer.viewer.biostructure, bio.menu.transform.fetch-pdb-sequences, biostructureviewer.panel.pdb-information, biostructureviewer.panel.pdb-id-viewer, biostructureviewer.panel.protein-ligand-interactions]
produced_from: ticket-review
manual_only_reason: |
  Every scenario depends on a service outside the platform: RCSB (structure download,
  GraphQL metadata and sequences) or the server-side Python script behind ProLIF. Run
  manually on a stand with outbound access to rcsb.org and a working Python script
  environment.
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as: []
---

# BiostructureViewer — Features that need RCSB or server-side Python (manual)

Structures loaded by PDB ID through the RCSB data providers, **Bio > Transform > Fetch PDB
Sequences...**, the PDB_ID context panels, and the **Protein-Ligand Interactions** (ProLIF)
panel.

## Setup

- Logged in; BiostructureViewer and Bio installed; the stand can reach rcsb.org; server-side
  Python scripts run.
- Files: `System:AppData/BiostructureViewer/pdb_id.csv` (column `pdb_id` with `1QBS`, `1ZP8`,
  `2BDJ`, `1IAN`, `4UJ1`, detected as PDB_ID) and
  `System:AppData/BiostructureViewer/pdb_data.csv` (`pdb_id` PDB_ID and `pdb` Molecule3D).
- Close all views before each scenario.

## Scenarios

### Scenario 1 — Biostructure viewer loads structures by PDB ID from RCSB

Steps:

1. Open **Browse > Files > App Data > BiostructureViewer** and double-click `pdb_id.csv`.
2. Click **Add viewer** in the toolbox, search `Biostructure` and click it.

   * Expected result: the viewer loads the structure of the current row (`1QBS`) and holds the
     Mol\* engine. In the settings, **Data > Biostructure Id** is `pdb_id` and
     **Data > Biostructure Data Provider** names one of the RCSB providers (RCSB PDB, RCSB
     mmCIF, RCSB bCIF). No error balloon.

3. Click row 2 (`1ZP8`) in the grid.

   * Expected result: the viewer loads `1ZP8`. No error balloon.

4. In the settings, pick another RCSB provider in **Biostructure Data Provider**.

   * Expected result: the viewer reloads the current structure. No error balloon.

### Scenario 2 — Fetch PDB Sequences adds chain sequence columns

Steps:

1. Open `pdb_id.csv`.
2. From the top menu choose **Bio > Transform > Fetch PDB Sequences...**.

   * Expected result: a dialog opens with a table input (the current table) and a column input
     limited to PDB_ID columns (`pdb_id`).

3. Keep the defaults and click **OK**.

   * Expected result: a progress indicator "Extracting protein sequences via GraphQL..." shows
     while sequences are fetched. Then columns **Chain 1**, **Chain 2**, … are appended (as
     many as the largest number of protein chains among the structures); each cell holds the
     amino-acid sequence of that chain; rows with fewer chains leave the extra cells empty. The
     new columns are shown as Macromolecule sequences (sequence renderer), not plain text. No
     error balloon.

### Scenario 3 — Running Fetch PDB Sequences again does not overwrite columns

Steps:

1. On the table from Scenario 2, run **Bio > Transform > Fetch PDB Sequences...** again on
   `pdb_id` and click **OK**.

   * Expected result: a second set of chain columns is added with new names
     (**Chain 1 (2)**, **Chain 2 (2)**, …). `pdb_id` and the first set of chain columns are
     unchanged. No error balloon.

### Scenario 4 — PDB_ID cell: PDB Information and PDB id viewer panels

Steps:

1. Open `pdb_data.csv` and click the `pdb_id` cell of row 1 (`1QBS`).
2. In the context panel, expand **PDB Information**.

   * Expected result: after a loader, the section shows metadata for `1QBS` fetched from RCSB.
     No console error.

3. Expand **PDB id viewer**.

   * Expected result: the section holds an NGL host (`.d4-ngl-viewer`) that loads `1QBS` from
     `files.rcsb.org`. No warning balloon "Can't load structure for PDB ID".

4. Click the `pdb_id` cell of row 2 (`1ZP8`).

   * Expected result: both sections are rebuilt for `1ZP8`. No console error.

### Scenario 5 — Protein-Ligand Interactions (ProLIF)

Steps:

1. In `pdb_data.csv`, click the `pdb` cell of row 1 (1QBS, which has a non-water ligand).
2. In the context panel, expand **Protein-Ligand Interactions**.

   * Expected result: after the server script finishes, the section shows an interaction
     diagram between the protein and its ligand. No error balloon.

3. Click the `pdb_id` cell of row 1.
4. Expand **Protein-Ligand Interactions**.

   * Expected result: the structure is fetched from RCSB and the section shows an interaction
     diagram; no text "Could not fetch PDB 1QBS". No error balloon.
