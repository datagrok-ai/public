---
feature: biostructureviewer
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [biostructure-prolif-cross-cell-type]
realizes: [biostructureviewer.biostructure]
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions:
  - id: SR-01
    check: E-SCENARIO-RUNTIME-ALIGNMENT
    rationale: |
      Scenario 1's `non-blank canvas` pixel-content Expected bullet is NOT
      pixel-asserted because the recon environment fails to create a WebGL
      rendering context (Mol* logs "Could not create a WebGL rendering
      context" — same pattern as sibling biostructure-viewer /
      property-surface / ngl-viewer specs). Structural mount IS asserted:
      `.bsv-container-info-panel` present + data-source `Biostructure
      Viewer:3D Structure` exact + accordion aria-expanded flips on header
      click.
    verdict_status: SCOPE_REDUCTION
  - id: SR-02
    check: E-SCENARIO-RUNTIME-ALIGNMENT
    rationale: |
      Scenario 4 Path B (Molecule3D + isAutoDockPose -> dockingInteractions
      Widget) is deferred per the scenario Deferral note: the fixture corpus
      has no AutoDock pose Molecule3D row, and synthesizing one requires a
      co-installed Docking package + receptor AppData under
      System:AppData/Docking/targets/. The spec asserts the registration's
      existence at the DG.Func.find surface (the isAutoDockPose condition
      predicate), capturing the three-condition-gated registration invariant.
    verdict_status: SCOPE_REDUCTION
  - id: SR-03
    check: E-SCENARIO-RUNTIME-ALIGNMENT
    rationale: |
      Scenario 1's re-mount sub-step is asserted at the JS-API panel-function
      re-invocation level (a fresh structure3D widget for a different cell
      value yields a fresh widget root) rather than at the live-DOM accordion
      re-render level — re-mount is WebGL-uncertain in the recon environment,
      same root cause as SR-01.
    verdict_status: SCOPE_REDUCTION
  - id: SR-04
    check: E-SCENARIO-RUNTIME-ALIGNMENT
    rationale: |
      Scenario 5's "selecting a SMILES column in the picker updates the
      widget's internal state" Expected bullet is asserted at the
      picker-structure level (input-host-SMILES-column present + select
      options filtered to the Molecule semType column `ligand`) rather than
      via a UI selection event, per the scenario's Scenario 5 note scoping the
      test to the direct JS-API call surface (the contract).
    verdict_status: SCOPE_REDUCTION
realized_as:
  - context-panel-widgets-extension-spec.ts
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
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T16:30:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T17:36:00Z
    spec_runs:
      - spec: context-panel-widgets-extension-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 198
        failure_keys: []
---

# BiostructureViewer — Context-panel widgets extension (3D Structure / PDB Information / PDB id viewer / ProLIF / Link With Molecule Column)

Covers the right-hand context panel widgets that appear when the current
grid cell has semantic type `Molecule3D` or `PDB_ID`, extending beyond the
main smoke test and the bug-focused regression guards. Five panel widgets
are exercised:

- `3D Structure` (Molecule3D) — inline Mol\* render of the structure.
- `PDB Information` (Molecule3D) — header / file info derived from the PDB
  string.
- `PDB Information` (PDB_ID) — async metadata fetched for the PDB ID.
- `Protein-Ligand Interactions` (ProLIF) — renders differently depending on
  the cell: a ligand-bound Molecule3D structure, an AutoDock pose, or a
  PDB_ID (fetched from RCSB).
- `Link With Molecule Column` — a cross-package widget that lets the
  Docking package's AutoDock info panel inject a "Link SMILES column"
  checkbox.

## Setup

1. Sign in as a standard test user with rights to use the Files browser and
   open the App Data area. The BiostructureViewer package ships sample PDB /
   Molecule3D fixtures under `tables/` / `files/`; reuse those (or any small
   PDB string) for the Molecule3D rows.
2. Open a fresh empty workspace (no stale viewers / tables / panels docked).
3. Prepare a small DataFrame with three columns that exercise the panel
   semTypes:
   - one column with `semType=Molecule3D` (PDB strings; ensure at least one
     row's PDB has a non-water HETATM ligand — most ligand-bound structures
     such as `1QBS` qualify);
   - one column with `semType=PDB_ID` (e.g. `1QBS`, `1BNA`, `2J1X`);
   - one column with `semType=Molecule` (SMILES strings) used by Scenario 5
     for `Link With Molecule Column`.
   Datagrok auto-detects Molecule3D / PDB_ID via the BiostructureViewer
   detectors. The same DataFrame is reused across all five scenarios.
4. Ensure the right-hand context panel is expanded (toggle from the gutter
   if collapsed). The context panel is where every widget in this scenario
   mounts.

## Scenarios

### Scenario 1: `3D Structure` panel — inline Mol* render for a Molecule3D cell

Steps:
1. Open the prepared DataFrame from Setup step 3.
2. Click any cell in the Molecule3D column to make it the current cell.
3. Observe the right-hand context panel.
4. Locate the `3D Structure` panel section (collapsible header).
5. Wait for the inline Mol* render to complete (the panel widget mounts a
   `PdbGridCellRendererBack`-backed viewer; allow the awaitRendered settle
   pattern documented for the main Biostructure viewer to apply here as
   well — poll for the `.msp-viewport canvas` node inside the panel root,
   or use a short settle delay).
6. Click a different Molecule3D cell.
7. Observe the panel.

Expected:
- The `3D Structure` panel section appears in the right context panel for
  the current Molecule3D cell.
- Inside the panel section, a Mol* viewport mounts (a `<canvas>` element
  inside an `.msp-viewport` container; the panel root carries the
  `bsv-container-info-panel` class added by `structure3D` after render).
- The render is non-blank — the canvas has non-zero pixel content, NOT the
  blank dark viewport with only the axis gizmo (that means "not rendered
  yet" / "parse failed" per the documented common pitfall).
- Switching the current cell to a different Molecule3D value re-mounts the
  panel against the new cell value (the loader appears briefly, then the
  new structure renders).
- No JS console errors during mount / re-mount.

### Scenario 2: `PDB Information` panel on a Molecule3D cell — header from PDB string

Steps:
1. From the prepared DataFrame, click any Molecule3D cell.
2. In the right context panel, locate the `PDB Information` section.
3. Expand the section if collapsed.
4. Observe the rendered widget content (header / file-info derived from the
   PDB string by `pdbFileInfoWidget(molecule.value)`).
5. Click a different Molecule3D cell.

Expected:
- A `PDB Information` panel section is present in the context panel for the
  current Molecule3D cell.
- The widget renders some non-empty content derived from the PDB string
  (header records, file-info such as HEADER / TITLE / COMPND lines, or
  whatever `pdbFileInfoWidget` extracts).
- Switching to a different Molecule3D cell refreshes the widget against the
  new PDB string.
- No JS console errors.
- IMPORTANT: this panel registration uses the same display name
  (`PDB Information`) as the PDB_ID-targeted panel in Scenario 3. The two
  registrations are differentiated by the param semType — assert by panel
  position relative to the current cell semType, not by name match.

### Scenario 3: `PDB Information` panel on a PDB_ID cell — async `pdbInfoWidget`

Steps:
1. From the prepared DataFrame, click a cell in the PDB_ID column (e.g. the
   row whose PDB_ID value is `1QBS`).
2. In the right context panel, locate the `PDB Information` section (now
   the PDB_ID-flavoured registration — `pdbInfoPanel` ->
   `pdbInfoWidget(pdbId)`).
3. Wait for the async widget body to assemble (the widget call is `async`;
   the panel briefly shows a loader before metadata appears).
4. Observe the rendered metadata.
5. Click a different PDB_ID cell.

Expected:
- A `PDB Information` panel section is present in the context panel for the
  current PDB_ID cell.
- The async widget resolves and renders metadata for the cell's PDB ID
  (whatever the upstream `pdbInfoWidget(pdbId)` produces — typically header
  / classification / chain count fields).
- Switching to a different PDB_ID cell triggers the widget's fetch /
  re-render against the new ID.
- No `Could not fetch PDB <id>` error string in the widget body (that would
  indicate the upstream fetch failed; treat as a network-blip flake and
  retry once if observed transiently).
- No JS console errors.

### Scenario 4: `Protein-Ligand Interactions` (ProLIF) — three condition-gated registrations resolve per cell type

This widget has three registrations under the same panel name, each gated
by a condition on the current cell. The scenario exercises all three
resolution paths from the single panel name.

Steps:
1. From the prepared DataFrame, click a Molecule3D cell whose PDB string
   contains a non-water HETATM ligand (the gating predicate
   `BiostructureViewer:hasNonWaterHetatm(molecule)` returns true).
2. In the right context panel, locate the `Protein-Ligand Interactions`
   section.
3. Wait for the ProLIF widget body to render (the widget builds a LigNetwork
   via `makeProlifWidget`; allow async render to settle).
4. Observe the rendered interactions diagram.
5. Switch the current cell to a Molecule3D cell whose value is an AutoDock
   pose (PDB containing a `REMARK ... binding energy` line — the gating
   predicate `BiostructureViewer:isAutoDockPose(molecule)` returns true).
   If the prepared DataFrame does not contain an AutoDock pose, defer this
   sub-step to the cited deferral below.
6. Observe the `Protein-Ligand Interactions` panel — confirm it is the
   Docking-flavour registration (receptor pre-fetched from
   `System:AppData/Docking/targets/<receptor>/`).
7. Switch the current cell to a PDB_ID cell (e.g. `1QBS`).
8. Observe the `Protein-Ligand Interactions` panel — confirm it is the
   PDB_ID-flavour registration (`pdbIdInteractionsWidget` fetches the PDB
   from RCSB via `grok.dapi.fetchProxy`).

Expected:
- For the Molecule3D + hasNonWaterHetatm cell: the ProLIF panel renders an
  interactions LigNetwork diagram via `pdbInteractionsWidget` →
  `makeProlifWidget({protein})`.
- For the AutoDock-pose cell (if exercised): the ProLIF panel renders the
  Docking-flavour widget via `dockingInteractionsWidget` →
  `makeDockingProlifWidget({protein: receptor, ligand: pose})`. The
  receptor is pre-fetched from Docking AppData on mount; the inline diagram
  shows the pose-against-receptor interactions.
- For the PDB_ID cell: the ProLIF panel renders via
  `pdbIdInteractionsWidget` after a `dapi.fetchProxy` to
  `https://files.rcsb.org/download/<id>.pdb`. If the fetch is non-ok, the
  widget body reads `Could not fetch PDB <id>`; treat that as flake-retry
  once.
- Each registration mounts under the same display name
  (`Protein-Ligand Interactions`) but the source registration differs per
  cell semType + condition. The condition gating is the load-bearing check.
- No JS console errors during any of the three registrations' mounts.

Deferral (cited): if the test fixture corpus does not include an AutoDock
pose Molecule3D row (a pose PDB with a `REMARK ... binding energy` line and
a corresponding receptor under `System:AppData/Docking/targets/`), the
Docking-flavour sub-step is deferred — the Docking package owns the
test-fixture setup for AutoDock poses (receptor deposit happens on Docking
package install). Real technical dependency: the fixture requires a
co-installed Docking package + receptor AppData files; tests cannot
synthesize an AutoDock pose ad-hoc inside the BiostructureViewer corpus.

### Scenario 5: `Link With Molecule Column` — cross-package widget injection

This is a function-widget exposed at the BiostructureViewer package level
so that the Docking package's AutoDock info panel can inject a "Link SMILES
column" checkbox at the bottom of its own widget (instead of creating a
second top-level panel section with the same name).

Steps:
1. Open the prepared DataFrame (Setup step 3) which carries both a
   Molecule3D column and a Molecule (SMILES) column.
2. Call the function directly via the JS API to exercise its widget shape:
   `await grok.functions.call('BiostructureViewer:mol3dAtomPickerLinkWidget',
   {mol3DCol: <the Molecule3D column ref>})`.
3. Observe the returned `DG.Widget` — embed it in a transient host (e.g.
   dock as a sidebar via `ui.dialog().add(widget.root).show()` or attach
   `widget.root` to an in-test DOM container) for visual inspection.
4. Interact with the "Link SMILES column" UI inside the widget (the
   `getMol3DAtomPickerLinkWidget` helper produces a SMILES-column picker
   tied to the Molecule3D column reference).

Expected:
- `grok.functions.call('BiostructureViewer:mol3dAtomPickerLinkWidget', {...})`
  resolves to a `DG.Widget` (the call returns successfully — the function
  registers at package load with `outputs: [{name: 'result', type:
  'widget'}]`).
- The widget root contains a Molecule-column picker scoped to the SMILES
  semType columns in the host DataFrame.
- Selecting a SMILES column in the picker updates the widget's internal
  state (the helper's contract — the linkage is in-memory, the widget does
  not need to be re-rendered).
- No JS console errors during the call or during picker interaction.
- IMPORTANT: this function is a **cross-package helper** — its primary
  consumer is the Docking package's AutoDock info panel, which calls it
  by package-qualified name. The scenario exercises only the direct
  JS-API call surface (the contract); the Docking-side injection flow is
  out of scope for the BiostructureViewer atlas (it belongs to the
  Docking section's tests when that atlas exists).

## Notes

- Scenario 4's three ProLIF registrations resolve to different widgets
  depending on the current cell (a ligand-bound structure, an AutoDock
  pose, or a PDB_ID) — the condition gating on the cell is the load-bearing
  behavior under test.
- Source citations:
  - See: .claude/skills/grok-browser/references/viewers/biostructureviewer.md#L160
    (Context-Panel Info Panels reference).
  - See: public/packages/BiostructureViewer/src/package.ts#L854 —
    `structure3D` (`3D Structure` panel registration for Molecule3D).
  - See: public/packages/BiostructureViewer/src/package.ts#L878 —
    `pdbFileInfoPanel` (`PDB Information` panel for Molecule3D).
  - See: public/packages/BiostructureViewer/src/package.ts#L247 —
    `pdbInfoPanel` (`PDB Information` async panel for PDB_ID).
  - See: public/packages/BiostructureViewer/src/package.ts#L262 —
    `pdbInteractionsWidget` (ProLIF panel for Molecule3D +
    hasNonWaterHetatm).
  - See: public/packages/BiostructureViewer/src/package.ts#L295 —
    `dockingInteractionsWidget` (ProLIF Docking-flavour for AutoDock
    pose).
  - See: public/packages/BiostructureViewer/src/package.ts#L317 —
    `pdbIdInteractionsWidget` (ProLIF panel for PDB_ID; fetches PDB
    via `dapi.fetchProxy` from RCSB).
  - See: public/packages/BiostructureViewer/src/package.ts#L893 —
    `mol3dAtomPickerLinkWidget` (`Link With Molecule Column`
    cross-package helper).
