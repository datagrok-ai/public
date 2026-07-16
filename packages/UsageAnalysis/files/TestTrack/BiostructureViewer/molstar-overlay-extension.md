---
feature: biostructureviewer
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: []
realizes: []
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T14:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T14:30:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T17:47:00Z
    spec_runs:
      - spec: molstar-overlay-extension-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 53
        failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T18:30:00Z
    review_round: 1
    failure_keys: []
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
      - check: A-STRUCT-MECH-02
        status: PASS
      - check: A-STRUCT-MECH-03
        status: PASS
      - check: A-STRUCT-MECH-04
        status: PASS
      - check: A-STRUCT-MECH-05
        status: PASS
      - check: A-STRUCT-MECH-06
        status: PASS
      - check: A-STRUCT-03
        status: PASS
      - check: A-STRUCT-04
        status: PASS
      - check: A-LAYER-ALIGN-01
        status: PASS
      - check: A-CONT-01
        status: PASS
      - check: A-BUG-01
        status: PASS
      - check: A-MERIT-01
        status: PASS
      - check: A-MERIT-02
        status: PASS
realized_as:
  - molstar-overlay-extension-spec.ts
---

# BiostructureViewer — Mol* viewport overlay buttons extension (Screenshot / Toggle Controls / Selection Mode / Settings)

Covers the Mol\* viewport overlay button row that floats over the 3D
viewport at the bottom right of every Biostructure viewer instance. The
row carries five buttons:

- `Reset Camera` — restores default camera orientation/zoom (already
  covered by the main smoke test; not part of this scenario).
- `Screenshot / State Snapshot` — opens the screenshot / state-snapshot
  panel for image export and state persistence.
- `Toggle Controls Panel` — shows/hides the Mol\* left controls panel
  (structure tree, etc.); mirrors the `layoutShowControls` Datagrok
  property.
- `Settings / Controls Info` — opens the Mol\* settings / controls info
  panel.
- `Toggle Selection Mode` — enters/exits atom/residue selection mode for
  in-viewport picking.

This scenario covers the four buttons not already exercised by the smoke
test; each opens a distinct overlay panel or mode.

## Setup

1. Sign in as a standard test user with rights to use the Files browser.
2. Open a fresh empty workspace (no stale viewers / tables / panels docked).
3. Load any small DataFrame with a `Molecule3D` column whose values include
   at least one protein-with-ligand PDB string (`1QBS` is the canonical
   reference; any PDB whose `HETATM` lines include a non-water ligand
   works). The BiostructureViewer package ships sample Molecule3D fixtures
   under `tables/` / `files/` — reuse those or a small ad-hoc DataFrame
   built via `DG.DataFrame.fromColumns` with a Molecule3D-tagged string
   column.
4. Add the Biostructure viewer to the current table view — either by
   selecting the `Biostructure` entry from the **Add viewer** dropdown (via
   the `+` Toolbox button), or programmatically via
   `tv.addViewer('Biostructure')`. The viewer docks over the table; its
   `[name="viewer-Biostructure"]` container is the scope root for the
   selectors below.
5. Wait for the structure to render before exercising overlay buttons —
   await `viewer.awaitRendered(timeoutMs)` (or poll for `.msp-viewport
   canvas` plus a settle delay). A blank dark viewport with only the
   axis gizmo means "not rendered yet" or "parse failed"; the overlay
   buttons themselves render before the structure, so they ARE clickable
   when blank — but Screenshot would capture an empty canvas, which is
   not the intended assertion.

The four scenarios below are independent — each starts from the rendered
viewer (Setup 1–5) and the others' DOM side-effects are not assumed.
Between scenarios, close any panel the prior scenario opened and clear
any toggled mode (this keeps assertions independent).

## Scenarios

### Scenario 1: `Screenshot / State Snapshot` overlay button — opens the screenshot panel

Steps:
1. With the Biostructure viewer rendered (Setup 1–5), locate the Mol*
   viewport overlay button row at the bottom-right corner of the
   `.msp-viewport` container (selector
   `.msp-viewport-controls-buttons button`).
2. Find the button with `title="Screenshot / State Snapshot"` (exact
   match against the `title` attribute; the universal-refdoc table is
   authoritative).
3. Click the button.
4. Observe the Mol* UI for a screenshot / state-snapshot panel mount —
   Mol* opens a screenshot configuration panel inline (typically
   below the viewport overlay row) with image-export options and a
   state-snapshot section.
5. Close the panel by clicking the same overlay button a second time
   (toggle off; the button uses `msp-btn-link-toggle-on|off` class
   pairing) or by clicking the panel's close `×` if present.
6. Observe the panel dismisses.

Expected:
- The `[title='Screenshot / State Snapshot']` button is present inside the
  `.msp-viewport-controls-buttons` row.
- Clicking the button mounts a Mol* screenshot / state-snapshot panel
  (look for a panel DOM node carrying `msp-` class prefixing — e.g.
  `.msp-control` containers — that is absent before the click and present
  after).
- The button class flips from `msp-btn-link-toggle-off` to
  `msp-btn-link-toggle-on` (or analogous toggled state) on first click.
- Clicking the button a second time toggles the panel back off (button
  class returns to off-state; panel DOM node detaches).
- No JS console errors during open / close.
- The structure underneath remains rendered (the screenshot panel is an
  overlay, not a replacement; `.msp-viewport canvas` still carries non-zero
  pixel content).

### Scenario 2: `Toggle Controls Panel` overlay button — toggles the Mol* left controls panel

Steps:
1. With the Biostructure viewer rendered (Setup 1–5), and the Mol* left
   controls panel in its default state (collapsed by default for new
   Biostructure viewer instances per the documented
   "side panels collapsed (3D viewport only)" default — verify
   `.msp-layout-show-left` is NOT present on the plugin root, or
   `.msp-layout-hide-left` is present).
2. Locate the overlay button with `title="Toggle Controls Panel"`.
3. Click the button.
4. Observe the Mol* left controls panel mounts — the structure tree and
   Mol*-native control sections appear in a left-side panel of the plugin
   root.
5. Click the button a second time.
6. Observe the left panel collapses back.
7. (Cross-check Datagrok property mirroring) Open the Datagrok settings
   panel (click `[name="icon-font-icon-settings"]` on the viewer title
   bar, or press F4) and inspect the `Layout` category — the
   `layoutShowControls` property reflects the same on/off state that the
   overlay button is in.

Expected:
- The `[title='Toggle Controls Panel']` button is present in the overlay
  row.
- Clicking the button flips a class on the plugin root that controls left
  panel visibility (e.g. `.msp-layout-show-left` toggles on, or the
  layout regions reflow to expose the left controls column).
- The left controls panel content (structure tree, builtin Mol* controls)
  renders inside the newly visible panel.
- Clicking the button again collapses the panel back to the default
  3D-viewport-only layout.
- The Datagrok `Layout` category property `layoutShowControls` mirrors the
  overlay button state.
- No JS console errors during open / close.

### Scenario 3: `Toggle Selection Mode` overlay button — enters/exits atom/residue selection mode

Steps:
1. With the Biostructure viewer rendered (Setup 1–5), and no Mol* mode
   currently active (the viewport's default interaction mode is camera
   navigation — left-drag rotates the camera, no selection occurs on
   click).
2. Locate the overlay button with `title="Toggle Selection Mode"`.
3. Click the button.
4. Observe the button class flips to the on-state
   (`msp-btn-link-toggle-on`).
5. Left-click on a residue or atom in the rendered structure within
   `.msp-viewport canvas`.
6. Observe the clicked residue / atom highlights (Mol* paints a
   selection-mode highlight on the picked element).
7. Click the overlay button a second time to exit selection mode.
8. Observe the button class flips back to off-state; clicks on the
   viewport no longer pick residues, the camera-navigation default
   resumes.

Expected:
- The `[title='Toggle Selection Mode']` button is present in the overlay
  row.
- Clicking the button toggles its class between off and on states
  (`msp-btn-link-toggle-off` ↔ `msp-btn-link-toggle-on`).
- In selection mode (on-state), clicking inside `.msp-viewport canvas`
  picks the underlying residue / atom and highlights it (look for Mol*
  selection paint — typically a green/highlight color overlay on the
  picked element).
- In camera mode (off-state), clicking inside the viewport rotates /
  pans the camera as the default Mol* behaviour without picking anything.
- Toggling back to off-state restores camera navigation; any prior
  selection highlight either clears or persists per Mol*'s default
  selection-retention behaviour (assertion is on the mode flip, not on
  the selection persistence policy).
- No JS console errors during enter / exit.

### Scenario 4: `Settings / Controls Info` overlay button — opens the Mol* settings panel

Steps:
1. With the Biostructure viewer rendered (Setup 1–5), locate the overlay
   button with `title="Settings / Controls Info"`.
2. Click the button.
3. Observe the Mol* settings / controls info panel mounts as an overlay
   on the viewport (a Mol*-native panel listing setting categories and
   controls reference — distinct from the Datagrok settings panel
   reachable via `[name="icon-font-icon-settings"]` on the title bar).
4. Inspect the panel content — it lists Mol*-native settings (e.g.
   render quality, background color, layout flags) and a controls-info
   reference (keyboard/mouse bindings inside the Mol* viewport).
5. Click the same overlay button a second time (or click the panel's
   close `×` if present) to dismiss.

Expected:
- The `[title='Settings / Controls Info']` button is present in the
  overlay row.
- Clicking the button mounts a Mol*-native settings / controls-info
  panel inline (look for a `.msp-`-prefixed panel DOM node that is
  absent before the click and present after).
- The mounted panel is the Mol* engine's own settings UI — distinct
  from the Datagrok property panel (`[name="viewer-Biostructure"]
  [name="icon-font-icon-settings"]` opens the Datagrok panel; this
  overlay button opens the Mol* panel).
- The button class flips to the on-state on click and back to off-state
  on the second click.
- The structure underneath continues to render (the settings panel is an
  overlay, not a replacement; `.msp-viewport canvas` still carries
  non-zero pixel content).
- No JS console errors during open / close.

## Notes

- Note that `biostructure.overlay.reset-camera` (the fifth overlay
  button) is intentionally excluded from this scenario — it is already
  covered by the main smoke test.
- Source citations:
  - See: .claude/skills/grok-browser/references/viewers/biostructureviewer.md#L75
    (Mol* viewport overlay buttons reference, with the `title`-attribute
    table for all five buttons).
  - See: public/packages/BiostructureViewer/src/viewers/molstar-viewer/molstar-viewer.ts#L258 —
    `layoutShowControls` Datagrok property (cross-referenced by
    Scenario 2 to assert overlay-button-to-Datagrok-property mirroring).
