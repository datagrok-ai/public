---
feature: biostructureviewer
target_layer: playwright
coverage_type: edge
priority: p2
realizes_atlas: [biostructure-ligand-overlay-row-driven, biostructure-binding-site-overlay]
realizes: [biostructureviewer.biostructure]
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

Covers the Mol\* (Biostructure) viewer's property surface beyond what the
main smoke test exercises: Data properties (`dataJson`, `pdb`, `pdbTag`),
the remaining Behaviour ligand toggles (`showMouseOverRowLigand`,
`showSelectedRowsLigands`), the Binding Site `bindingSiteWholeResidues`
boundary, the Layout property group, and the Controls property group —
together with the canonical safe JS API entry point
`viewBiostructure(content, format, name)`. Scenario 2 reproduces a
documented pitfall: passing a raw PDB string to the `pdb` property without
a name fails to parse; the scenario reproduces the failure and demonstrates
the safe recovery path.

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
  reproduced in Scenario 2).
- The viewer canvas is non-empty (axis gizmo plus visible structure
  geometry, not the blank dark "not rendered yet" canvas).
- The property panel's Data category shows `dataJson` populated (read-only
  view).
- No JS console errors during the swap.

### Scenario 2: Edge — raw `pdb` prop without a name fails to parse (the canonical pitfall)

This is the negative-path scenario: reproducing the documented failure mode
and asserting the safe entry point that avoids it.

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
  "not rendered yet" symptom.
- The recovery step (step 7) calling `viewBiostructure(content, format,
  name)` renders the structure successfully — the canonical safe entry
  point supplies the `name` parameter the raw-`pdb` path lacked.
- After recovery, no further `Parsed object is empty` errors appear.

### Scenario 3: Data property — `pdbTag` populated from a DataFrame tag

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

These are the two remaining Behaviour toggles not covered by the smoke
test's `showCurrentRowLigand` and `showBindingSite` flows.

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

This is the boundary property paired with `showBindingSite` +
`bindingSiteRadius`, which are covered by the smoke test.

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

Layout is a property group covering roughly a dozen individual Mol\* layout
flags; this scenario covers the two access paths — property panel and the
Mol\* viewport overlay button — that both flip the same flag
(`layoutShowControls`).

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

The Controls property group has two flags; this scenario verifies the
property panel surface and the default-off state of `showWelcomeToast`.

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

- Scenario 2 is the documented edge case; Scenarios 1 and 3 demonstrate
  the canonical safe entry points (`dataJson` and `pdbTag` paths) that
  contrast with it. The other scenarios cover the rest of the property
  surface as general breadth around that edge case.
- Source citations:
  - See: .claude/skills/grok-browser/references/viewers/biostructureviewer.md
    (viewers reference doc).
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
