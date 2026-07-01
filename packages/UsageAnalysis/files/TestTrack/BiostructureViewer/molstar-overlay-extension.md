---
feature: biostructureviewer
sub_features_covered:
  - biostructure.overlay.screenshot
  - biostructure.overlay.toggle-controls
  - biostructure.overlay.selection-mode
  - biostructure.overlay.settings-info
target_layer: playwright
coverage_type: regression
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

Coverage-extension scenario authored under cycle
`2026-06-04-biostructureviewer-migrate-01` Gate F SR routing. Extends the
section beyond the Mol*-centric smoke, the five bug-focused regression guards,
the NGL-viewer breadth extension, the Property-surface + edge extension, and
the Context-panel widgets extension to the **Mol* viewport overlay button
family** that floats over the 3D viewport at the bottom right of every
Biostructure viewer instance. The overlay button row carries five buttons
addressed by their `title` attribute (per the universal-refdoc table at
`.claude/skills/grok-browser/references/viewers/biostructureviewer.md#L75-L82`):

- `Reset Camera` — restores default camera orientation/zoom (ALREADY COVERED
  by the smoke scenario via `biostructure.overlay.reset-camera`; not part of
  this extension's scope)
- `Screenshot / State Snapshot` — opens the screenshot / state-snapshot panel
  for image export and state persistence
- `Toggle Controls Panel` — shows/hides the Mol* left controls panel
  (structure tree, etc.); mirrors the `layoutShowControls` Datagrok property
- `Settings / Controls Info` — opens the Mol* settings / controls info panel
- `Toggle Selection Mode` — enters/exits atom/residue selection mode for
  in-viewport picking

The four overlay buttons exercised here all share the Mol* DOM idiom
(`.msp-viewport-controls-buttons button[title='<NAME>']`, classes
`msp-btn msp-btn-icon msp-btn-link-toggle-on|off`) but each opens a
distinct overlay panel / mode that no prior scenario asserts.

All 4 sub_features below are net-new to the live covered-union of 45 (none
of `overlay.screenshot`, `overlay.toggle-controls`, `overlay.selection-mode`,
`overlay.settings-info` appears in any of the section's nine prior
scenarios' `sub_features_covered:` lists; the only `overlay.*` id already
covered is `overlay.reset-camera`, which lives in the smoke and is
intentionally excluded from this scenario's list). `round_net_new = +4`;
the section's projected covered-union climbs from 45 → 49 after this
scenario lands.

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
   axis gizmo means "not rendered yet" or "parse failed" per the
   documented common pitfall (see atlas `edge_cases[6]`); the overlay
   buttons themselves render before the structure, so they ARE clickable
   when blank — but Screenshot would capture an empty canvas, which is
   not the intended assertion.

The four scenarios below are independent — each starts from the rendered
viewer (Setup 1–5) and the others' DOM side-effects are not assumed.
Between scenarios, close any panel the prior scenario opened and clear
any toggled mode (this keeps assertions independent).

## Scenarios

### Scenario 1: `Screenshot / State Snapshot` overlay button — opens the screenshot panel

Exercises `biostructure.overlay.screenshot`.

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

Exercises `biostructure.overlay.toggle-controls`.

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
  overlay button state (atlas declares this mirroring under
  `biostructure.prop.layout`'s interactions: "click Mol* overlay 'Toggle
  Controls Panel' to flip layoutShowControls").
- No JS console errors during open / close.

### Scenario 3: `Toggle Selection Mode` overlay button — enters/exits atom/residue selection mode

Exercises `biostructure.overlay.selection-mode`.

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

Exercises `biostructure.overlay.settings-info`.

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

- **target_layer rationale**: All four scenarios are UI-driven — they
  click overlay DOM buttons inside `.msp-viewport-controls-buttons`,
  assert overlay panel mounts (Mol*-native DOM-attached widgets),
  assert WebGL canvas content (camera vs selection mode), and
  cross-check Datagrok property mirroring (Scenario 2's `layoutShowControls`
  cross-check). None of these surfaces is reachable from `apitest`-only
  JS-API calls: the overlay buttons are not exposed as a Datagrok
  function, and the picking / mode toggle behaviours live entirely in
  the Mol* engine's DOM-and-WebGL surface. Hence `target_layer: playwright`.
- **coverage_type rationale**: `regression`. The four scenarios are
  general coverage of the Mol* viewport overlay button family; they do
  not map onto any of the atlas's 7 `edge_cases[]` entries (which are
  bug-anchored to file-handler / project-persistence / multi-ligand /
  whitespace-right-click / view-removal surfaces, or to the two
  documentation pitfalls). They are not a critical_path / golden path
  (no `priority: p0`/`p1` mapping in atlas `critical_paths[]`; the
  section's existing smoke at `biostructure-viewer.md` already covers
  the camera-navigation surface and the `Reset Camera` overlay button),
  so `smoke` is inappropriate. They are not boundary-value or
  atlas-edge_case, so `edge` is inappropriate. They are not stress /
  latency-sensitive, so `perf` is inappropriate. `regression` is the
  STEP E heuristic default for general coverage of a thematic feature
  surface.
- **Bug coverage**: `related_bugs: []`. Gate F's F-BUG-COVERAGE-01 has
  already closed in this cycle (all 5 known bugs are realized as
  bug-focused scenarios under the section). This is a pure
  breadth-extension scenario; no bug-library entry overlaps the four
  overlay-button sub_features above.
- **net_new**: Against `live_covered_union` (45 ids) and atlas
  `manual_only[]` (empty), all 4 listed sub_features are net-new — none
  appears in any of the section's nine prior scenarios'
  `sub_features_covered:` lists. `round_net_new = +4`. The STEP C
  net-new refusal is satisfied (the breadth-loop progress bound holds:
  covered_union 45 → 49 strictly rises by 4). Note that
  `biostructure.overlay.reset-camera` is INTENTIONALLY excluded from
  this scenario's `sub_features_covered:` — it is already covered by
  the smoke; including it here would not contribute net-new value and
  would also create an unnecessary scope overlap with the smoke
  scenario's overlay coverage.
- **breadth-loop residual**: After this scenario, the section's
  covered-union projects to 49/72 = 68.1%, still just under the 70%
  threshold. Remaining high-yield candidates per the chain
  `gate_f_verdict.gaps[sub-feature-coverage-gap]` (post-this-scenario):
  Biotrack +2 (`biostructure.biotrack-viewer` + `.props`), Grid
  context-menu extension +1 (`biostructure.grid-context-menu.open-pdb-residues`),
  File-preview Mol* +4 (`file-preview` + `file-preview.molstar-structure` +
  `file-preview.molstar-topology` + `file-preview.molstar-density`),
  File-open extension +1 (`file-open.importXYZ`), Data-provider
  extension +2 (`data-provider.rcsb-pdb` + `data-provider.rcsb-bcif`),
  JS API extension +2 (`api.viewPdbById` + `api.viewPdbByData`),
  `cell-renderer.pdb-id` +1. Cumulative residual net-new after this
  scenario lands: +13. Coverage projects to (49+13)/72 = 86.1% — well
  above the 70% threshold. Next breadth round can keep sequencing
  these.
- **Atlas citations**:
  - See: .claude/skills/grok-browser/references/viewers/biostructureviewer.md#L75
    (universal-refdoc Mol* viewport overlay buttons section with the
    five-row `title`-attribute table; the four buttons this scenario
    addresses are rows L78–L81. Resolved via atlas `ui_reference_doc:`
    since the path lives off the default convention under `viewers/`).
  - See: .claude/skills/grok-browser/references/viewers/biostructureviewer.md#L78 —
    `Screenshot / State Snapshot` overlay button row (atlas anchor for
    `biostructure.overlay.screenshot`).
  - See: .claude/skills/grok-browser/references/viewers/biostructureviewer.md#L79 —
    `Toggle Controls Panel` overlay button row (atlas anchor for
    `biostructure.overlay.toggle-controls`).
  - See: .claude/skills/grok-browser/references/viewers/biostructureviewer.md#L80 —
    `Settings / Controls Info` overlay button row (atlas anchor for
    `biostructure.overlay.settings-info`).
  - See: .claude/skills/grok-browser/references/viewers/biostructureviewer.md#L81 —
    `Toggle Selection Mode` overlay button row (atlas anchor for
    `biostructure.overlay.selection-mode`).
  - See: public/packages/BiostructureViewer/src/viewers/molstar-viewer/molstar-viewer.ts#L258 —
    `layoutShowControls` Datagrok property (cross-referenced by
    Scenario 2 to assert overlay-button-to-Datagrok-property
    mirroring; declared on atlas under `biostructure.prop.layout`).
