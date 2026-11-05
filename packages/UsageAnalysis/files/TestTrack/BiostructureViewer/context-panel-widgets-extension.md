---
feature: biostructureviewer
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: []
realizes: [biostructureviewer.panel.3d-structure, biostructureviewer.panel.pdb-information]
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - context-panel-widgets-extension-spec.ts
---

# BiostructureViewer — Context panel for a Molecule3D cell (3D Structure / PDB Information)

When the current grid cell is a Molecule3D value, the context panel shows the package panels
**3D Structure** (an embedded Biostructure viewer for that cell) and **PDB Information**
(header data read from the PDB text itself, without any network call).

Panels that need outside services — **PDB Information** and **PDB id viewer** on a PDB_ID
cell, **Protein-Ligand Interactions** (ProLIF) — are in `biostructureviewer-network-ui.md`.

## Setup

- Logged in; BiostructureViewer installed.
- Open **Browse > Files > App Data > BiostructureViewer** and double-click `pdb_data.csv`
  (6 rows; `pdb` is a Molecule3D column; row 1 is 1QBS, row 2 is 1ZP8).
- The context panel is open (toggle it from the sidebar if it is collapsed).

## Scenarios

### Scenario 1 — 3D Structure panel for a Molecule3D cell

Steps:

1. Click the `pdb` cell of row 1.
2. In the context panel, expand **3D Structure**.

   * Expected result: after a short loader, the section holds a viewer root with the class
     `bsv-container-info-panel` and the Mol\* engine (`.msp-plugin`) inside it. No console
     error.

3. Click the `pdb` cell of row 2.

   * Expected result: the **3D Structure** section is rebuilt for the new cell: the loader
     shows, then a `bsv-container-info-panel` root with the Mol\* engine. No console error.

### Scenario 2 — PDB Information panel for a Molecule3D cell

Steps:

1. Click the `pdb` cell of row 1.
2. In the context panel, expand **PDB Information**.

   * Expected result: under **General**, **Classification** is `ASPARTYL PROTEASE`,
     **Description** is `HIV-1 PROTEASE INHIBITORS WIIH LOW NANOMOLAR POTENCY` (the typo is in
     the source file), **PDB URL** is the link `rcsb.org/structure/1QBS`.

3. Click the `pdb` cell of row 2.

   * Expected result: **Classification** is `HYDROLASE`, **Description** is
     `HIV PROTEASE WITH INHIBITOR AB-2`, **PDB URL** is `rcsb.org/structure/1ZP8`. No console
     error.

## Automation notes

- A Molecule3D cell whose structure has a non-water ligand also gets the
  **Protein-Ligand Interactions** section, which runs a server-side Python script; it is not
  expanded or asserted here.
- Scenario 2 reads only the **General** table; the collapsible sub-sections below it are not
  asserted.
