---
feature: biostructureviewer
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: []
realizes: [biostructureviewer.viewer.biostructure]
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - molstar-overlay-extension-spec.ts
---

# BiostructureViewer — Mol* viewport overlay buttons (Screenshot / Toggle Controls / Selection Mode / Settings / Binding site / Expanded viewport)

The Mol\* viewport of every Biostructure viewer carries a row of overlay buttons
(`.msp-viewport-controls-buttons`). Each button is found by its tooltip (`title`) and shows its
state through the classes `msp-btn-link-toggle-on` / `msp-btn-link-toggle-off`. The package
adds its own **Binding site** button, linked both ways with the **Binding Site** settings, and
lets Escape leave the expanded viewport. **Reset Camera** is covered by the smoke test.

## Setup

1. Open **Browse > Files > App Data > BiostructureViewer** and double-click `pdb_data.csv`
   (6 rows: `pdb_id` detected as PDB_ID, `pdb` detected as Molecule3D).
2. Click **Add viewer** in the toolbox, search `Biostructure` and click it. The viewer takes
   the `pdb` column by itself (settings **Data > Biostructure Id Column Name** = `pdb`) and
   shows row 1 (1QBS, HIV-1 protease with an inhibitor).
3. Wait until the overlay button row (`.msp-viewport-controls-buttons`) is present in the
   viewer.

Each scenario starts from this state; close whatever the previous scenario opened.

## Scenarios

### Scenario 1 — Screenshot / State Snapshot

Steps:

1. Click the overlay button with `title="Screenshot / State Snapshot"`.

   * Expected result: the button class becomes `msp-btn-link-toggle-on`; a Mol\* panel
     (an `.msp-` element absent before the click) opens with the screenshot and state
     snapshot options.

2. Click the same button again.

   * Expected result: the button class returns to `msp-btn-link-toggle-off`; the panel is
     gone. No console error.

### Scenario 2 — Toggle Controls Panel

Steps:

1. Click the overlay button with `title="Toggle Controls Panel"`.

   * Expected result: the button class becomes `msp-btn-link-toggle-on`; the Mol\* control
     panels (structure tree and controls) appear next to the viewport.

2. Click the same button again.

   * Expected result: the button class returns to `msp-btn-link-toggle-off`; the control
     panels are hidden and only the 3D viewport is shown. No console error.

3. Open the viewer settings (gear icon), category **Layout**, and turn
   **Layout Show Controls** on.

   * Expected result: the Mol\* control panels appear; the **Toggle Controls Panel** button
     class is `msp-btn-link-toggle-on`.

4. Turn **Layout Show Controls** off.

   * Expected result: the control panels are hidden. No console error.

### Scenario 3 — Toggle Selection Mode

Steps:

1. Click the overlay button with `title="Toggle Selection Mode"`.

   * Expected result: the button class becomes `msp-btn-link-toggle-on`.

2. Click the same button again.

   * Expected result: the button class returns to `msp-btn-link-toggle-off`. No console
     error.

### Scenario 4 — Settings / Controls Info

Steps:

1. Click the overlay button with `title="Settings / Controls Info"`.

   * Expected result: the button class becomes `msp-btn-link-toggle-on`; a Mol\* settings
     panel (an `.msp-` element absent before the click) opens over the viewport. It is not
     the Datagrok settings panel of the viewer.

2. Click the same button again.

   * Expected result: the button class returns to `msp-btn-link-toggle-off`; the panel is
     gone. No console error.

### Scenario 5 — Binding site button stays in sync with the Binding Site settings

Steps:

1. Find the overlay button with `title="Binding site"`.

   * Expected result: the button exists and is enabled (its title is "Binding site", not
     "No ligand detected"); its class is `msp-btn-link-toggle-off`.

2. Click it.

   * Expected result: a popover opens with the header "Binding Site", a **Show side chains**
     checkbox (unchecked) and a **Radius** slider showing `5.0 Å`.

3. In the popover, check **Show side chains**.

   * Expected result: the button class becomes `msp-btn-link-toggle-on`.

4. Open the viewer settings (gear icon), category **Binding Site**.

   * Expected result: **Show Binding Site** is checked.

5. In the settings, set **Binding Site Radius** to `8`.
6. Click the **Binding site** button to open the popover again if it has closed.

   * Expected result: the popover shows `8.0 Å`.

7. In the settings, uncheck **Show Binding Site**.

   * Expected result: the popover checkbox **Show side chains** is unchecked; the button
     class is `msp-btn-link-toggle-off`.

8. With the popover open, press Escape.

   * Expected result: the popover closes. No console error.

### Scenario 6 — Toggle Expanded Viewport, then Escape

Steps:

1. Click the overlay button with `title="Toggle Expanded Viewport"`.

   * Expected result: the Mol\* content carries the class `msp-layout-expanded`; the button
     class is `msp-btn-link-toggle-on`.

2. Press Escape.

   * Expected result: `msp-layout-expanded` is gone; the button class is
     `msp-btn-link-toggle-off`; the table view is still open (Escape closed nothing else).
     No console error.

## Automation notes

- Scenario 2: by code, the **Layout Show Controls** property drives the Mol\* layout, but no
  code writes the button state back into the property; whether the property follows the
  button is not verified on the stand and is not asserted.
- Scenario 3: what the selection mode highlights in the structure is not asserted, only the
  button state.
- Scenario 5: the popover is attached to `document.body` (class `bsv-bs-popover`), not to the
  viewer; an outside click closes it, which is why step 6 may need to reopen it.
- Scenario 6: the viewer does not pass the Mol\* "show expand" option, so the button relies on
  the Mol\* default; an earlier live check listed five overlay buttons without
  **Toggle Expanded Viewport**. Its presence is to be checked on the stand.
- The button tooltips come from the bundled `@rcsb/rcsb-molstar` version; match them once on
  the stand after a Mol\* upgrade.
