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
realized_as:
  - ngl-viewer-extension-spec.ts
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
    timestamp: 2026-06-04T16:44:32Z
    spec_runs:
      - spec: ngl-viewer-extension-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 142
        failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T17:30:00Z
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
---

# BiostructureViewer — NGL viewer extension

Covers the sibling NGL viewer surface, beyond the main smoke test and the
bug-focused guards: file-handler routing for structure formats Mol* cannot
open, NGL preview panes in the Files browser, the NGL-specific grid
context-menu entry, the NGL property surface, and the PDB_ID context-panel
widget that wraps an NGL canvas.

## Setup

1. Sign in as a standard test user with rights to use the Files browser and
   open the App Data area for BiostructureViewer test fixtures (the package
   ships several structure files under `tables/` / `files/` per its
   `package.json` — those are the source of the `.mmtf`, `.ply`, `.ccp4`, `.cns`,
   and `.pdb` fixtures referenced below).
2. Open a fresh empty workspace (no stale viewers / tables docked).
3. For Scenarios 2 and 5, prepare a small DataFrame with one column of
   `Molecule3D` semType (a PDB string column) and one column of `PDB_ID`
   semType (e.g. `1QBS`, `1BNA`) — Datagrok auto-detects both via the
   BiostructureViewer detectors. Reuse this DataFrame across both scenarios.

## Scenarios

### Scenario 1: NGL viewer add via JS API + property surface (Style + Data + Behaviour)

Steps:
1. Open a small table containing a Molecule3D column (PDB strings).
2. Call `tv.addViewer('NGL')` to dock the NGL viewer.
3. Wait for the NGL canvas to render (poll for the NGL viewer's WebGL canvas
   element; a settle delay matches the awaitRendered pattern documented for
   Mol*).
4. Open the property panel for the NGL viewer (click the settings icon on the
   title bar or press `F4`).
5. Under the Style category, change `representation` to a non-default value
   (e.g. `ball+stick`).
6. Under the Data category, set `ligandColumnName` to the Molecule3D column.
7. Under the Behaviour category, toggle `showCurrentRowLigand` off then on.
8. Toggle `showMouseOverRowLigand` off then on.
9. Toggle `showSelectedRowsLigands` on, then select two rows in the grid.

Expected:
- NGL canvas appears with a rendered structure (non-empty WebGL surface).
- Representation change re-renders the structure in the new style without
  console errors.
- Ligand-column wiring resolves: with `showCurrentRowLigand=true`, the current
  row's ligand value is overlaid; toggling redraws the overlay.
- `showSelectedRowsLigands=true` overlays ligands for the two selected rows.
- No JS console errors during any property panel interaction.

### Scenario 2: NGL file-handler routing — extensions Mol* cannot open

Steps:
1. In the Files browser, navigate to a folder containing structure files of
   formats Mol* cannot parse — at least one each of `.mmtf`, `.cns`,
   `.prmtop`, and `.ccp4` (any subset that the test fixture provides; one of
   each is preferred).
2. Double-click the `.mmtf` file.
3. Confirm an NGL viewer opens (not a Biostructure / Molstar viewer).
4. Close the view.
5. Repeat the double-click step for `.cns`, `.prmtop`, and `.ccp4` in turn,
   closing the view between each.

Expected:
- For every Mol*-incapable extension (`.mmtf`, `.cns`, `.prmtop`, `.ccp4`),
  the file-handler search resolves to `importWithNgl` and the NGL viewer is
  the viewer that mounts (NOT the Mol*/Biostructure viewer).
- No `Parsed object is empty` console error fires (that error would indicate
  the file was incorrectly routed to the Mol* handler).
- No file-handler-search collision behaviour resembling GROK-14442 (which
  guards `.pdb` vs `.pdbqt` disambiguation in the Mol* path).

### Scenario 3: NGL file pre-viewers — structure / surface / density

Steps:
1. In the Files browser, single-click a `.mmtf` file (do NOT double-click —
   single-click renders the preview pane, double-click opens via the file
   handler).
2. Observe the preview pane.
3. Single-click a `.ply` file (mesh / surface format).
4. Observe the preview pane.
5. Single-click a `.ccp4` file (density format).
6. Observe the preview pane.

Expected:
- `.mmtf` preview pane shows an NGL structure render via
  `previewNglStructure`.
- `.ply` preview pane shows an NGL surface render via `previewNglSurface`.
- `.ccp4` preview pane shows an NGL density render via `previewNglDensity`.
- Each preview engages NGL (not Mol*) because the file extension falls in the
  NGL-preview routing table.
- No console error during preview engagement.

### Scenario 4: Grid context menu — Show NGL Viewer on a Molecule3D cell

This is the NGL twin of the "Show Biostructure Viewer" flow covered by
`biostructureviewer-bug-grok-14552.md`.

Steps:
1. Open the prepared Molecule3D table from Setup step 3.
2. Right-click a populated Molecule3D cell to open the grid-cell context menu.
3. Click "Show NGL Viewer".
4. Wait for the NGL canvas to render with the cell's structure.
5. Close the docked NGL viewer.

Expected:
- The NGL viewer is docked over the cell's structure.
- The structure that renders matches the right-clicked cell's value (not the
  current row's structure if the two differ).
- No console error.
- The `Show Biostructure Viewer` menu item (covered separately by the
  GROK-14552 regression guard) and `Show NGL Viewer` are both present in the
  cell-context-menu for Molecule3D cells, side-by-side (the two registrations
  coexist).

### Scenario 5: PDB id context-panel widget — PDB id viewer (NGL wrapped)

Steps:
1. Use the prepared DataFrame from Setup step 3 (PDB_ID column with at least
   one valid PDB ID, e.g. `1QBS`).
2. Click a cell in the PDB_ID column to make it the current cell.
3. Open the right-hand context panel (if collapsed).
4. Wait for the `PDB id viewer` panel widget to render.
5. Click a different PDB_ID cell.

Expected:
- The `PDB id viewer` panel widget renders an NGL widget UI scoped to the
  current PDB_ID cell value.
- Switching the current cell to a different PDB_ID re-renders the widget for
  the new ID (or triggers the widget's fetch / refresh sequence).
- No console error during the widget mount or re-mount.

## Notes

- Source citations:
  - See: .claude/skills/grok-browser/references/viewers/biostructureviewer.md
    (viewers reference doc).
  - See: public/packages/BiostructureViewer/src/package.ts#L142 — `importPdb`
    family registration; sibling `importWithNgl` at #L168.
  - See: public/packages/BiostructureViewer/src/package.ts#L190 — NGL preview
    file-viewer family registration.
  - See: public/packages/BiostructureViewer/src/utils/context-menu.ts#L62 —
    `Show NGL Viewer` cell context-menu registration.
  - See: public/packages/BiostructureViewer/src/viewers/ngl-viewer.ts#L63 —
    NGL viewer property surface.
  - See: public/packages/BiostructureViewer/src/package.ts#L238 — PDB id
    viewer panel widget.
