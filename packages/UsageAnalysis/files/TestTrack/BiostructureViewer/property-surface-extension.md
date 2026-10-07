---
feature: biostructureviewer
target_layer: playwright
coverage_type: edge
priority: p2
realizes_atlas: [biostructureviewer.int.binding-site-overlay]
realizes: [biostructureviewer.viewer.biostructure]
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - property-surface-extension-spec.ts
---

# BiostructureViewer — Settings: Binding Site whole residues and the Controls group

Covers two parts of the Biostructure viewer settings not exercised elsewhere: the
**Binding Site Whole Residues** switch and the **Controls** category (**Show Welcome Toast**,
**Show Import Controls**). The **Layout Show Controls** switch is covered in
`molstar-overlay-extension.md`.

## Setup

1. Open **Browse > Files > App Data > BiostructureViewer** and double-click `pdb_data.csv`.
2. Click **Add viewer** in the toolbox, search `Biostructure` and click it. Wait until the
   viewer holds the Mol\* engine (`.msp-plugin`).
3. Open the viewer settings (gear icon on the viewer title bar).

## Scenarios

### Scenario 1 — Binding Site Whole Residues

Steps:

1. In the category **Binding Site**, check **Show Binding Site**.

   * Expected result: **Show Binding Site** is checked; **Binding Site Whole Residues** is
     checked (default). No error balloon.

2. Uncheck **Binding Site Whole Residues**.

   * Expected result: the settings show it unchecked. No error balloon, no console error.

3. Check **Binding Site Whole Residues** again, then uncheck **Show Binding Site**.

   * Expected result: the settings show the new values. No error balloon, no console error.

### Scenario 2 — Controls category

Steps:

1. Open the category **Controls**.

   * Expected result: it holds **Show Welcome Toast** and **Show Import Controls**, both
     unchecked.

2. Check **Show Import Controls**, then uncheck it.

   * Expected result: the switch follows each click. No console error.

## Automation notes

- Scenario 1: the difference between whole residues and single atoms is visible only in the
  drawing; only the property values and the absence of errors are asserted.
- Scenario 2: what Mol\* shows when **Show Import Controls** is on is not checked; it is up to
  Mol\*.
