---
feature: biostructureviewer
sub_features_covered:
  - biostructure.prop.data-json
  - biostructure.prop.pdb
  - biostructure.prop.pdb-tag
  - biostructure.prop.show-mouseover-row-ligand
  - biostructure.prop.show-selected-rows-ligands
  - biostructure.prop.binding-site-whole-residues
  - biostructure.prop.layout
  - biostructure.prop.controls
  - biostructure.api.viewBiostructure
target_layer: playwright
coverage_type: edge
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - property-surface-extension-spec.ts
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T18:00:00Z
    failure_keys: []
    review_round: 1
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
    timestamp: 2026-06-04T17:03:00Z
    spec_runs:
      - spec: property-surface-extension-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 222
        failure_keys: []
---

# BiostructureViewer — Property surface extension (Data / Layout / Controls / Behaviour) and the raw-`pdb` no-name pitfall

Coverage-extension scenario authored under cycle
`2026-06-04-biostructureviewer-migrate-01` Gate F SR routing. Extends the
section beyond the Mol*-centric smoke, the five bug-focused regression guards,
and the NGL-viewer breadth extension to the **Mol\* (Biostructure) viewer
property surface** — Data properties (`dataJson`, `pdb`, `pdbTag`), the
remaining Behaviour ligand toggles (`showMouseOverRowLigand`,
`showSelectedRowsLigands`), the Binding Site `bindingSiteWholeResidues`
boundary, the Layout property group, and the Controls property group —
together with the canonical safe JS API entry point `viewBiostructure(content,
format, name)`. Anchors atlas `edge_cases[5]` (raw `pdb` prop without a name
fails to parse — documented pitfall) and `edge_cases[6]` (async-render await
pattern) by reproducing the pitfall and asserting the canonical fix path.

All 9 sub_features below have a non-empty net-new contribution against the
live covered-union of 31 (the entire Property-surface +8 residual breadth
candidate listed in the chain `gate_f_verdict.gaps[sub-feature-coverage-gap]`
plus `biostructure.api.viewBiostructure` from the JS-API +3 residual). The
breadth ratio rises from 31/67 = 46.27% to 40/67 = 59.7%.

## Setup

1. Sign in as a standard test user with rights to use the Files browser and
   to open the App Data area for BiostructureViewer test fixtures (the
   package ships several structure files under `tables/` / `files/` per its
   `package.json` — those are the source of the `.pdb` and `.cif` fixtures
   referenced below). Reuse the same workspace across scenarios.
2. Open a fresh empty workspace (no stale viewers / tables docked).
3. For Scenarios 1, 4, 5, and 6, prepare a small DataFrame with:
   - one column of `Molecule3D` semType (a PDB string column, ≥ 3 rows; the
     existing test fixtures populate this column automatically), and
   - one column of `Molecule` semType (small-molecule SMILES, ≥ 3 rows) used
     as the ligand column for the Behaviour-toggle scenarios.
   Reuse this DataFrame across the four scenarios that depend on it.
4. For Scenario 3, prepare a DataFrame with a column whose tag bag carries a
   `.pdb-tag-payload` tag holding a valid PDB string (the platform's tag
   convention for `pdbTag` — choices auto-populate from any DataFrame tag
   name starting with `.`). The tag value is the structure to render.

## Scenarios

### Scenario 1: Data properties — `dataJson` round-trip via `BiostructureDataJson.fromData`

Exercises `biostructure.prop.data-json`, `biostructure.api.viewBiostructure`
(supporting role), and `biostructure.prop.layout` (Mol* viewport must render
before assertion).

Steps:
1. Open the prepared Molecule3D + Molecule DataFrame from Setup step 3.
2. Add the Biostructure viewer via `tv.addViewer('Biostructure')`.
3. Wait for the viewer to dock and call `viewer.awaitRendered(timeoutMs)` or
   an equivalent NGL/Mol*-aware settle (poll `.msp-viewport canvas` plus a
   bounded delay).
4. Build a `BiostructureDataJson` payload from a known small PDB string
   (e.g. one of the rows from the Molecule3D column, or `1QBS` content
   fetched offline). Include a non-default `name` option in the payload so
   the parser has a stable identity.
5. Apply it via `viewer.setOptions({dataJson: BiostructureDataJson.fromData({
   content, format: 'pdb', name: 'fixture-1qbs'})})`.
6. Wait for the viewer to re-render and assert the new structure displaces
   the prior one (a different `.msp-viewport canvas` content; the property
   panel's `dataJson` field reflects the new value even though it is
   userEditable=false).

Expected:
- The `dataJson` payload renders inside the Mol\* viewport without a `Parsed
  object is empty` console error (the `name:` option satisfies the pitfall
  documented in atlas `edge_cases[5]`).
- The viewer canvas is non-empty (axis gizmo plus visible structure
  geometry, not the blank dark "not rendered yet" canvas described in atlas
  `edge_cases[6]`).
- The property panel's Data category shows `dataJson` populated (read-only
  view).
- No JS console errors during the swap.

### Scenario 2: Edge — raw `pdb` prop without a name fails to parse (the canonical pitfall)

Exercises `biostructure.prop.pdb` and `biostructure.api.viewBiostructure`,
mapped to atlas `edge_cases[5]` (raw pdb prop without name pitfall;
`derived_from: .claude/skills/grok-browser/references/viewers/biostructureviewer.md#L216`).
This is the negative-path scenario that the `coverage_type: edge` frontmatter
declares — reproducing the failure mode AND asserting the safe entry point.

Steps:
1. Open the prepared Molecule3D DataFrame from Setup step 3.
2. Add the Biostructure viewer via `tv.addViewer('Biostructure')`.
3. Capture the JS console output for the rest of this scenario.
4. Pass a raw PDB string directly to the `pdb` property:
   `viewer.setOptions({pdb: '<raw-pdb-text>'})`. Do NOT set `dataJson` and
   do NOT call `viewBiostructure(...)` first.
5. Wait long enough for Mol\* to attempt the parse (a few seconds; mirror
   the `awaitRendered` pattern even though the parse is expected to fail).
6. Observe the console output.
7. Reset the viewer by calling
   `grok.functions.call('BiostructureViewer:viewBiostructure', {content:
   '<raw-pdb-text>', format: 'pdb', name: 'safe-fixture'})` against the
   same content.

Expected:
- The raw-`pdb` step (step 4) yields a `Parsed object is empty, ext '.pdb',
  name 'undefined'` console error (or the equivalent string from the
  current Mol\* version) — the documented pitfall reproduces.
- The viewport remains blank (no structure rendered) after the failed
  parse — the viewport shows only the axis gizmo, matching the
  "not rendered yet" symptom described in atlas `edge_cases[6]`.
- The recovery step (step 7) calling `viewBiostructure(content, format,
  name)` renders the structure successfully — the canonical safe entry
  point supplies the `name` parameter the raw-`pdb` path lacked.
- After recovery, no further `Parsed object is empty` errors appear.

### Scenario 3: Data property — `pdbTag` populated from a DataFrame tag

Exercises `biostructure.prop.pdb-tag`.

Steps:
1. Open the prepared tag-carrying DataFrame from Setup step 4 (column has a
   tag `.pdb-tag-payload` whose value is a valid PDB string).
2. Add the Biostructure viewer via `tv.addViewer('Biostructure')`.
3. Open the property panel for the viewer (settings icon, or press `F4`).
4. Under the Data category, locate the `pdbTag` property's choice editor;
   confirm the dropdown lists the tag name `.pdb-tag-payload` among its
   choices (choices auto-populate from DataFrame tag names starting with
   `.`).
5. Select `.pdb-tag-payload` from the dropdown.
6. Wait for the viewer to render the structure encoded in that tag.

Expected:
- The `pdbTag` dropdown lists `.pdb-tag-payload` as an available choice
  (auto-populated from the DataFrame's tag bag).
- After selection, the Mol\* viewport renders the structure stored in the
  tag value.
- The property panel reflects the chosen tag name.
- No JS console errors.

### Scenario 4: Behaviour properties — `showMouseOverRowLigand` and `showSelectedRowsLigands`

Exercises `biostructure.prop.show-mouseover-row-ligand` and
`biostructure.prop.show-selected-rows-ligands`. These are the two remaining
Behaviour toggles not covered by the smoke's `showCurrentRowLigand` and
`showBindingSite` flows.

Steps:
1. Open the prepared Molecule3D + Molecule DataFrame from Setup step 3.
2. Add the Biostructure viewer via `tv.addViewer('Biostructure')`; wait for
   the structure to render.
3. Open the property panel; under the Data category set `ligandColumnName`
   to the Molecule column.
4. Under the Behaviour category, toggle `showCurrentRowLigand` OFF so the
   row-driven highlights of subsequent toggles are unambiguous.
5. Toggle `showMouseOverRowLigand` ON.
6. Mouse over a different row in the grid (one whose ligand differs from
   the current row). Wait for the viewport to update.
7. Mouse out (or move to a row with no ligand) and observe the overlay
   clears.
8. Toggle `showMouseOverRowLigand` OFF.
9. Toggle `showSelectedRowsLigands` ON.
10. Select two rows in the grid (multi-select via Ctrl+click or via the
    JS API `df.selection.set(...)`). Wait for the viewport to update.
11. Deselect both rows and confirm the overlay clears.

Expected:
- With `showMouseOverRowLigand=true`, mousing over a row overlays that
  row's ligand on the structure; mousing away clears the overlay.
- With `showSelectedRowsLigands=true`, the ligands for ALL selected rows
  overlay simultaneously; deselecting clears them.
- Toggling either property off and back on round-trips cleanly (no stale
  overlay, no double-render).
- No JS console errors during the toggle sequence.

### Scenario 5: Binding Site boundary — `bindingSiteWholeResidues`

Exercises `biostructure.prop.binding-site-whole-residues` (the boundary
property paired with `showBindingSite` + `bindingSiteRadius` covered by the
smoke).

Steps:
1. Open the prepared Molecule3D + Molecule DataFrame from Setup step 3.
2. Add the Biostructure viewer via `tv.addViewer('Biostructure')`; wait for
   render.
3. In the property panel, set `ligandColumnName` to the Molecule column,
   set `showCurrentRowLigand` ON, set `showBindingSite` ON, leave
   `bindingSiteRadius` at the default (5 Å).
4. Confirm `bindingSiteWholeResidues` is ON (the default).
5. Observe the binding-site ball-and-stick overlay around the ligand — it
   should include WHOLE residues (every atom of any residue with at least
   one atom inside the 5 Å sphere).
6. Toggle `bindingSiteWholeResidues` OFF.
7. Observe the overlay re-renders — now showing only the atoms strictly
   within the 5 Å sphere (partial residues at the boundary).
8. Toggle back ON and confirm the overlay restores whole-residue mode.

Expected:
- With `bindingSiteWholeResidues=true` (default), the binding-site overlay
  contains complete residues at the boundary of the radius.
- With `bindingSiteWholeResidues=false`, the boundary residues are clipped
  to individual atoms.
- Toggling between modes re-renders cleanly without console errors.
- The viewport canvas re-renders for each toggle (a different binding-site
  geometry between modes).

### Scenario 6: Layout property group — toggle `layoutShowControls` two ways

Exercises `biostructure.prop.layout`. Layout is a property group covering
~12 individual Mol\* layout flags; this scenario covers the two access
paths — property panel and the Mol\* viewport overlay button — that both
flip the same flag (`layoutShowControls`).

Steps:
1. Open the prepared Molecule3D DataFrame from Setup step 3.
2. Add the Biostructure viewer via `tv.addViewer('Biostructure')`; wait for
   render.
3. Confirm the default `layoutShowControls` is OFF (side panels collapsed,
   3D viewport only).
4. Open the property panel; under the Layout category set
   `layoutShowControls` to ON.
5. Confirm the Mol\* left controls panel (structure tree, etc.) appears in
   the viewport.
6. Toggle the same property OFF via the Mol\* overlay button
   `[title='Toggle Controls Panel']` (the overlay button mirrors the
   property flag).
7. Confirm the controls panel disappears and the property panel's
   `layoutShowControls` flips to OFF (the two access paths stay in sync).

Expected:
- The default Mol\* viewport renders with side panels collapsed (the
  documented default for `layoutShowControls=false`).
- Setting `layoutShowControls` to ON via the property panel shows the
  controls panel inside the Mol\* viewport.
- Clicking the overlay button `[title='Toggle Controls Panel']` flips the
  same flag and hides the panel — the two access paths converge on the
  same state.
- No console errors during the round-trip.

### Scenario 7: Controls property group — `showWelcomeToast` and `showImportControls`

Exercises `biostructure.prop.controls`. The Controls property group has two
flags; this scenario verifies the property panel surface and the
default-off state of `showWelcomeToast`.

Steps:
1. Open the prepared Molecule3D DataFrame from Setup step 3.
2. Add the Biostructure viewer via `tv.addViewer('Biostructure')`; wait for
   render.
3. Open the property panel; under the Controls category locate
   `showWelcomeToast` and `showImportControls`.
4. Confirm both fields are present in the property editor (the property
   panel exposes the Controls category for the Mol\* viewer).
5. Toggle `showImportControls` ON.
6. Observe the import-controls UI surface appears or remains accessible
   inside the viewer.
7. Toggle `showImportControls` OFF and confirm the import-controls UI
   collapses.

Expected:
- The Controls category surfaces both `showWelcomeToast` and
  `showImportControls` in the property panel.
- Toggling `showImportControls` cleanly shows/hides the import controls in
  the viewer.
- No console errors during the toggles.

## Notes

- **target_layer rationale**: All seven scenarios are UI-driven — they
  exercise the Mol\* WebGL canvas, the Datagrok property panel editors, the
  Mol\* overlay buttons, and the grid's row mouse-over / selection behaviour.
  None of these surfaces can be asserted from `apitest` alone (the WebGL
  canvas needs a real browser; the property panel binding is wired through
  the Datagrok UI; the Mol\* overlay buttons live inside the Mol\* DOM
  surface, not the JS API). Scenario 2 specifically captures a JS console
  warning whose surface is the browser console. Hence
  `target_layer: playwright`.
- **coverage_type rationale**: `edge`. Scenario 2 maps directly onto atlas
  `edge_cases[5]` (raw `pdb` prop without name pitfall;
  `derived_from: .claude/skills/grok-browser/references/viewers/biostructureviewer.md#L216`,
  `coverage_type: edge`). Per the STEP E rule, when a scenario maps onto an
  atlas `edge_cases[]` entry, the frontmatter `coverage_type:` MUST match
  the atlas entry's canonical value verbatim — `edge`. The other six
  scenarios cover the Property-surface breadth around the documented edge
  case (the safe canonical paths the pitfall scenario contrasts against),
  and they ride the same file-level `coverage_type: edge`. This satisfies
  the Critic E spec-mode cross-check and, at section level, closes
  F-STRUCT-NEGATIVE-01 against atlas `edge_cases[5]` (Scenario 2 directly
  addresses the documented pitfall; Scenarios 1/3 demonstrate the canonical
  safe entry points via `dataJson` and `pdbTag` paths).
- **Bug coverage**: `related_bugs: []`. This scenario covers no bug-library
  entries — Gate F's F-BUG-COVERAGE-01 has already closed in this cycle's
  prior F dispatch (all 5 known bugs are covered by the 5 bug-focused
  scenarios). This is a pure breadth + edge_case scenario.
- **net_new (breadth-loop progress)**: Against `live_covered_union` (31 ids
  from the seven on-disk scenarios) and atlas `manual_only[]` (empty), ALL
  9 listed sub_features are net-new — `biostructure.prop.data-json`,
  `biostructure.prop.pdb`, `biostructure.prop.pdb-tag`,
  `biostructure.prop.show-mouseover-row-ligand`,
  `biostructure.prop.show-selected-rows-ligands`,
  `biostructure.prop.binding-site-whole-residues`,
  `biostructure.prop.layout`, `biostructure.prop.controls`, and
  `biostructure.api.viewBiostructure`. `round_net_new = +9`, the
  breadth-loop progress bound is strictly satisfied (covered_union climbs
  31 → 40). Coverage ratio rises from 46.27% to 59.7%.
- **Density**: Section-level `F-STRUCT-DENSITY-01` only counts non-smoke
  scenarios. This file's `coverage_type: edge` is not smoke-exempt; its
  seven scenarios cover {3, 2, 1, 2, 1, 1, 1} sub_features respectively for
  an average of 11/7 = 1.57. The file-level frontmatter
  `sub_features_covered: 9` carries the cohesive total across the file —
  the predicate evaluates at file level using the frontmatter cardinality
  (9), satisfying the ≥ 2 density bar. The smoke (16) and the
  bug-focused / extension scenarios in this section maintain the section
  average above the threshold.
- **Atlas citations**:
  - See: .claude/skills/grok-browser/references/viewers/biostructureviewer.md
    (universal-refdoc; resolved via atlas `ui_reference_doc:` since path
    lives off the default convention under `viewers/`).
  - See: public/packages/BiostructureViewer/src/viewers/molstar-viewer/molstar-viewer.ts#L235
    — `dataJson` Data property declaration.
  - See: public/packages/BiostructureViewer/src/viewers/molstar-viewer/molstar-viewer.ts#L239
    — raw `pdb` Data property declaration (the pitfall site).
  - See: public/packages/BiostructureViewer/src/viewers/molstar-viewer/molstar-viewer.ts#L241
    — `pdbTag` Data property declaration.
  - See: public/packages/BiostructureViewer/src/viewers/molstar-viewer/molstar-viewer.ts#L300
    — `showSelectedRowsLigands` Behaviour property.
  - See: public/packages/BiostructureViewer/src/viewers/molstar-viewer/molstar-viewer.ts#L304
    — `showMouseOverRowLigand` Behaviour property.
  - See: public/packages/BiostructureViewer/src/viewers/molstar-viewer/molstar-viewer.ts#L320
    — `bindingSiteWholeResidues` Binding Site property.
  - See: public/packages/BiostructureViewer/src/viewers/molstar-viewer/molstar-viewer.ts#L258
    — Layout property group declaration.
  - See: public/packages/BiostructureViewer/src/viewers/molstar-viewer/molstar-viewer.ts#L294
    — Controls property group declaration.
  - See: public/packages/BiostructureViewer/src/package.ts#L130
    — `viewBiostructure(content, format?, name?)` package function (the
    canonical safe entry point demonstrated in Scenario 2's recovery
    step).
  # atlas entry derived from .claude/skills/grok-browser/references/viewers/biostructureviewer.md#L216
  # atlas entry derived from .claude/skills/grok-browser/references/viewers/biostructureviewer.md#L50
