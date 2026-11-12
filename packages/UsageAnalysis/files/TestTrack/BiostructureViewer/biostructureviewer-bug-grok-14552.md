---
feature: biostructureviewer
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [GROK-14552]
realizes: [biostructureviewer.cell.molecule3d, biostructureviewer.viewer.biostructure, biostructureviewer.viewer.ngl]
produced_from: ticket-review
related_bugs:
  - GROK-14552
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - biostructureviewer-bug-grok-14552-spec.ts
---

# BiostructureViewer — Grid context menu on Molecule3D cells and on row whitespace (GROK-14552)

Right-clicking the empty area of a grid row past the last column used to raise
`Cannot read properties of null (reading 'semType')` from the package's context-menu hook.
On a Molecule3D cell the package adds **Copy**, **Download** and a **Show** group with
**Biostructure** and **NGL**.

## Setup

- Logged in; BiostructureViewer installed.
- Open **Browse > Files > App Data > BiostructureViewer** and double-click `pdb_data.csv`
  (6 rows; `pdb` is a Molecule3D column, `pdb_id` is a PDB_ID column).
- Make the grid wider than its columns, so there is empty space to the right of the last
  column.

## Scenarios

### Scenario 1 — Right-click the row whitespace past the last column

Steps:

1. Right-click row 2 in the empty area to the right of the last column.

   * Expected result: no error balloon; no console error containing `semType`. Whatever menu
     opens has no **Show > Biostructure** item.

2. Press Escape.

### Scenario 2 — Right-click a Molecule3D cell

Steps:

1. Right-click the `pdb` cell of row 1.

   * Expected result: the menu has the items **Copy** and **Download** and a group **Show**
     with **Biostructure** and **NGL**. No console error.

2. Click **Copy**.

   * Expected result: info balloon "Value copied to clipboard".

3. Right-click the `pdb` cell of row 2 and choose **Show > Biostructure**.

   * Expected result: a Biostructure viewer holding the Mol\* engine (`.msp-plugin`) is
     docked in the table view. No error balloon.

4. Right-click the `pdb` cell of row 2 and choose **Show > NGL**.

   * Expected result: an NGL viewer (`.d4-ngl-viewer` with a `canvas`) is docked in the table
     view. No error balloon.

5. Right-click the `pdb_id` cell of row 1.

   * Expected result: the menu has no **Show > Biostructure** or **Show > NGL** items (the
     package adds them to Molecule3D cells only).

## Automation notes

- Scenario 2 step 2: the clipboard needs a secure origin (https or localhost); on plain http
  the package shows the warning "The clipboard functionality requires a secure origin, either
  HTTPS or localhost." instead of the info balloon.
