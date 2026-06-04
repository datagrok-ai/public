---
feature: biostructureviewer
sub_features_covered:
  - biostructure.panel.structure-3d
  - biostructure.panel.pdb-file-info
  - biostructure.panel.pdb-info
  - biostructure.panel.prolif
  - biostructure.panel.link-molecule-column
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
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
    scope_reduction_proposal: |-
      Accept four documented scope reductions, each grounded in a real technical limitation cited in the scenario .md or by established sibling-spec precedent. (SR-01) Scenario 1's `non-blank canvas` pixel-content Expected bullet is NOT pixel-asserted because the dev recon environment fails to create a WebGL rendering context (Mol* engine logs `Could not create a WebGL rendering context` — same pattern as sibling biostructure-viewer-spec.ts / property-surface-extension-spec.ts SR-01 / ngl-viewer-extension-spec.ts SR-04). Structural mount IS deterministically asserted: `.bsv-container-info-panel` container class present + data-source=`Biostructure Viewer:3D Structure` exact + accordion-pane aria-expanded flips false to true on header click. (SR-02) Scenario 4 Path B (Molecule3D + isAutoDockPose -> dockingInteractionsWidget) is deferred per the scenario .md explicit Deferral note: the fixture corpus has no AutoDock pose Molecule3D row, and synthesizing one requires a co-installed Docking package + receptor AppData files under System:AppData/Docking/targets/ — a real cross-package technical dependency. The spec asserts the registration's existence at the DG.Func.find function-registry surface (panel decorator processed at package load), capturing the bug-invariant of `three condition-gated registrations under the same panel name` via the isAutoDockPose condition predicate string match. (SR-03) Scenario 1's re-mount sub-step is asserted at the JS-API panel-function re-invocation level (a fresh structure3D widget produced for a different cell value yields a fresh widget root) rather than at the live-DOM accordion re-render level — re-mount is WebGL-uncertain in the recon environment, same root cause as SR-01. (SR-04) Scenario 5's `Selecting a SMILES column in the picker updates the widget's internal state` Expected bullet is asserted at the picker-structure level (input-host-SMILES-column present + select options filtered to Molecule semType column `ligand`) rather than via a UI selection event, per the scenario .md Scenario 5 IMPORTANT note that explicitly scopes the test to `the direct JS-API call surface (the contract)`. All four SRs preserve the structural bug-invariant of the respective sub_features. Structural mechanics: E-STRUCT-MECH-01..06 all PASS (file exists; valid TS; one test() block; test name + softStep names traceable to scenario headings; imports from @playwright/test and the sibling `../spec-login` helper module; leading /* --- sub_features_covered: [...] --- */ block lists 5 ids all resolving to atlas sub_features[].id). Content checks: traceability is clean (each of 5 scenarios maps to >=1 softStep with deterministic expect() verifications); selector discipline (every locator [name="div-section--3D-Structure"], [name="div-section--PDB-Information"], [name="div-section--Protein-Ligand-Interactions"], [name="pane-3D-Structure"], [name="pane-PDB-Information"], [name="pane-PDB-id-viewer"], [name="pane-Protein-Ligand-Interactions"], [name="input-host-SMILES-column"], [name="Browse"]) is documented in the leading class-2 recon-notes block with live-MCP-observed-date 2026-06-04 on dev.datagrok.ai including the surface they belong to and how they were reached; helpers loginToDatagrok and softStep are registered, specTestOptions and stepErrors are sibling-module exports from `../spec-login` (co-located with the registered helpers, not reinventions); E-LAYER-COMPLIANCE-01 is satisfied (target_layer: playwright with >=1 DOM-driving call: [name="Browse"] waitFor + three accordion-pane header clicks for Scenarios 1/2/3/4); E-BOUND-01/02 clean (only writes under public/packages/UsageAnalysis/files/TestTrack/BiostructureViewer/). Gate B awareness: scenario carries gate_verdicts.b.verdict: FAIL same-cycle with failure_keys [B-RUN-PASS, B-STAB-01], placing this in retry context per the spec-mode §`Gate B awareness on retry context` rule. The block carries no `flake_evidence[]` payload, so Critic E cannot identify a specific code path Validator FAILed on that the new spec retains; the new spec body shows substantive work addressing assertable structural invariants over WebGL-uncertain pixel assertions (the four documented SRs), defensible JS-API fallbacks for re-mount uncertainty (SR-03), a function-registry probe with the precise condition predicate match for Path B (SR-02), and a bounded fetchProxy wait window for Path C — these are good-faith engineering changes to the test surface, not cosmetic edits. E-RETRY-IGNORES-GATE-B does NOT fire (per mode: `Critic E should NOT enforce convergence — that is Validator's job at Gate B on the next run`).
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

Coverage-extension scenario authored under cycle
`2026-06-04-biostructureviewer-migrate-01` Gate F SR routing. Extends the
section beyond the Mol*-centric smoke, the five bug-focused regression guards,
the NGL-viewer breadth extension, and the Property-surface + edge extension to
the **context-panel widget family** that surfaces under the right-hand context
panel when the current cell has `semType=Molecule3D` or `semType=PDB_ID`. Five
panel widgets register against these semTypes:

- `3D Structure` (Molecule3D) — inline Mol* render via `PdbGridCellRendererBack`
- `PDB Information` (Molecule3D) — header / file info derived from the PDB string
- `PDB Information` (PDB_ID) — async metadata assembled via `pdbInfoWidget(pdbId)`
- `Protein-Ligand Interactions` (ProLIF) — three registrations gated by condition:
  Molecule3D + hasNonWaterHetatm, Molecule3D + isAutoDockPose (Docking-flavour
  with receptor pre-fetch), PDB_ID (fetches PDB via `dapi.fetchProxy` from RCSB)
- `Link With Molecule Column` — cross-package function-widget exposed for
  Docking's AutoDock info panel to inject a "Link SMILES column" checkbox

All 5 sub_features below have a non-empty net-new contribution against the
live covered-union of 40 (none of `panel.structure-3d`, `panel.pdb-file-info`,
`panel.pdb-info`, `panel.prolif`, `panel.link-molecule-column` appears in any
of the section's eight prior scenarios' `sub_features_covered:` lists —
verified against the chain `gate_f_verdict.gaps[sub-feature-coverage-gap]`
listing). `round_net_new = +5`; the section's projected covered-union climbs
from 40 → 45 after this scenario lands.

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

Exercises `biostructure.panel.structure-3d`.

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

Exercises `biostructure.panel.pdb-file-info`.

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

Exercises `biostructure.panel.pdb-info`.

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

Exercises `biostructure.panel.prolif`. This widget has three registrations
under the same panel name, each gated by a `condition:` expression on the
input param. The scenario exercises all three resolution paths from the
single panel name.

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

Exercises `biostructure.panel.link-molecule-column`. This is a function-widget
exposed via `@grok.decorators.func({name: 'Link With Molecule Column'})` at
the BiostructureViewer package level so that the Docking package's AutoDock
info panel can inject a "Link SMILES column" checkbox at the bottom of its
own widget (instead of creating a second top-level panel section with the
same name).

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

- **target_layer rationale**: All five scenarios are UI-driven — they
  exercise the right-hand context panel (Scenarios 1–4) and DOM-attached
  widget shapes (Scenario 5). Scenario 1 needs the browser's WebGL surface
  for the inline Mol* canvas. Scenario 4's three registration paths are
  condition-gated by `BiostructureViewer:hasNonWaterHetatm` /
  `BiostructureViewer:isAutoDockPose` at panel-resolution time — that
  resolution lives on the panel host, not in JS-API space. Scenario 5 is
  the only one with a pure-JS-API callable shape, but its assertion is
  "the returned widget renders a SMILES picker" — a DOM assertion. Hence
  `target_layer: playwright`.
- **coverage_type rationale**: `regression`. The scenarios are general
  coverage of the context-panel widget surface; they do not map onto any
  of the atlas's 7 `edge_cases[]` entries (which are bug-anchored to the
  Mol* / file-handler / project-persistence surfaces or to the two
  documentation pitfalls — raw `pdb` prop without name, async-render
  await). They are not a critical_path / golden path (no `priority: p0`
  mapping; the section's existing smoke at `biostructure-viewer.md`
  already covers the canonical Mol* + Bio top-menu surfaces), so `smoke`
  is inappropriate. They are not boundary-value or atlas-edge_case, so
  `edge` is inappropriate. They are not stress / latency-sensitive, so
  `perf` is inappropriate. `regression` is the STEP E heuristic default
  for general coverage of a thematic feature surface.
- **Bug coverage**: `related_bugs: []`. Gate F's F-BUG-COVERAGE-01 has
  already closed in this cycle (all 5 known bugs are realized as
  bug-focused scenarios under the section). This is a pure
  breadth-extension scenario; no bug-library entry overlaps the five
  panel sub_features above.
- **net_new**: Against `live_covered_union` (40 ids) and atlas
  `manual_only[]` (empty), all 5 listed sub_features are net-new — none
  appears in any of the section's eight prior scenarios'
  `sub_features_covered:` lists. `round_net_new = +5`. The STEP C
  net-new refusal is satisfied (the breadth-loop progress bound holds:
  covered_union 40 → 45 strictly rises by 5).
- **breadth-loop residual**: After this scenario, the section's
  covered-union projects to 45/67 = 67.2%, still under the 70% threshold.
  Remaining high-yield candidates per the chain
  `gate_f_verdict.gaps[sub-feature-coverage-gap]`: Biotrack +2
  (`biostructure.biotrack-viewer` + `.props`), Grid context-menu
  extension +1 (`biostructure.grid-context-menu.open-pdb-residues`),
  File-preview Mol* +4 (`file-preview` +
  `file-preview.molstar-structure` + `file-preview.molstar-topology` +
  `file-preview.molstar-density`), File-open extension +1
  (`file-open.importXYZ`), Mol* overlay extension +4
  (`overlay.screenshot` + `overlay.toggle-controls` +
  `overlay.selection-mode` + `overlay.settings-info`), Data-provider
  extension +2 (`data-provider.rcsb-pdb` + `data-provider.rcsb-bcif`),
  JS API extension +2 (`api.viewPdbById` + `api.viewPdbByData`),
  `cell-renderer.pdb-id` +1. Cumulative residual net-new after this
  scenario lands: +17. Coverage projects to (45+17)/67 = 92.5% — well
  above the 70% threshold. Next breadth round can keep sequencing
  these.
- **Atlas citations**:
  - See: .claude/skills/grok-browser/references/viewers/biostructureviewer.md#L160
    (universal-refdoc Context-Panel Info Panels section; resolved via
    atlas `ui_reference_doc:` since the path lives off the default
    convention under `viewers/`).
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
