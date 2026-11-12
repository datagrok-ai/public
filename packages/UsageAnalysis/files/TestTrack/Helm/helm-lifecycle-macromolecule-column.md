---
feature: helm
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [macromolecule_helm_column]
realizes: [helm.cell.helm, helm.cell-editor.molecule, helm.panel.properties]
realized_as:
  - helm-lifecycle-macromolecule-column-spec.ts
related_bugs: []
---

# Helm — Macromolecule HELM column: render, edit, compute properties & convert

Runs the full lifecycle of a grid column detected as Macromolecule/HELM: the
grid renders each HELM sequence as connected monomer structures (with correct
re-rendering after scrolling and cache invalidation), the cell editor lets a
user edit and save a sequence back to the cell, the Properties context panel
computes formula / molecular weight / extinction coefficient (including a
safeguard for very long sequences), and the column converts to molfiles.
This targets the column lifecycle after the 2026 HELM editor rewrite
(JSDraw2/Dojo replaced by a native SVG editor); `getMolfiles` now returns an
HWE pseudo-molfile instead of a V2000/V3000 molfile, and the helper's `parse` /
`removeGaps` (the former Scenarios 5 and 6) are gone from `IHelmHelper`. See
grok-browser `references/helm.md`.

## Setup

1. Authenticate to Datagrok as the test user.
2. Dataset: `System.AppData/Helm/samples/HELM.csv` (platform-shipped
   540-row HELM peptide sample, with the numeric `Activity` column).
   **Open it directly via its instance-derived file URL** —
   `<instance>/file/System.AppData/Helm/samples/HELM.csv?browse=files` (e.g.
   `https://dev.datagrok.ai/file/...`) — so the platform and dataset open in a
   single navigation; do NOT open the platform first and then read the CSV.
   Opening the file triggers
   `grok.data.detectSemanticTypes`; the Bio Macromolecule detector
   classifies the `HELM` column with `semType=Macromolecule`,
   `units=helm`, `cell.renderer=helm`, `quality=Macromolecule` so the
   grid auto-applies `HelmGridCellRenderer`. The Helm package `@init`
   (`initHelm`) runs at platform startup; the new SVG editor has no
   Dojo/JSDraw2 cold-load, so opens are fast.
3. Datagrok monomer library: default platform monomer library loaded by
   Bio is sufficient. Do NOT swap the monomer lib during the run.

## Scenarios

### Scenario 1: render_helm_cell — renderer attaches on the Macromolecule HELM column and invalidates on tag change

Steps:
1. Open `HELM.csv` directly via its file URL (see Setup).
2. Wait for the table view to render. Confirm the `HELM` column shows
   monomer structures via `HelmGridCellRenderer` and the `Activity` column
   shows numeric values.
3. Scroll the grid by ~50 rows and back to row 0 to force a re-render of
   the same cells.
4. Re-apply the column tag via the API path
   (`column.setTag(DG.TAGS.CELL_RENDERER, 'helm')`) to force renderer
   cache invalidation, then scroll again.

Expected:
- The `HELM` column initially renders monomer structures in every visible
  row; no raw `PEPTIDE1{...}$$$$` text leaks through.
- After scrolling, re-rendered cells still show monomer structures — the
  `HelmGridCellRendererBack` LRU cache rehydrates without console errors.
- After the tag re-apply forces cache invalidation, cells re-render cleanly
  (the renderer reattaches; no orphaned bitmap state).

### Scenario 2: edit_helm_cell — cellEditor opens, edits round-trip through HELM string and commit to the cell

Steps:
1. With `HELM.csv` open from Scenario 1, double-click a non-empty HELM cell.
2. Confirm the full-screen HELM Web Editor dialog opens
   (`.d4-dialog.d4-dialog-full-screen`, `[data-testid="app-root"]`) with
   the sequence drawn as connected monomer structures
   (`[data-testid="editor-svg"]`). The bottom tabs read `Sequence`, `HELM`,
   `Properties` (no `Structure View`).
3. Click the bottom `HELM` tab (`[data-testid="tab-helm"]`) and note the
   raw HELM string in `[data-testid="notation-pane-content"]`.
4. Remove the trailing monomer from the raw HELM text and commit
   (input/Enter/blur) — this re-draws the structure from raw text inline
   without committing to the grid (there is no `Apply` button).
5. Click footer `OK` (`button[name="button-OK"]`).
6. Verify the grid cell now renders the edited HELM string (one fewer
   monomer than before); the committed value carries the `V2.0` suffix.
7. Cancel out of a re-opened dialog to confirm CANCEL discards pending edits.

Expected:
- Double-click activates `editMoleculeCell` (`role: cellEditor`, tag-gated
  on `quality=Macromolecule, units=helm`).
- Editing the notation pane redraws the dialog but does NOT mutate the grid
  (Pitfall #3).
- Footer `OK` calls `cell.setValue(helm)`; the grid cell renders the edited
  sequence on dialog close.
- `CANCEL` on a re-opened dialog leaves the cell unchanged.

### Scenario 3: compute_properties — Properties panel formula / MW / extinction coefficient with the >1000-char guard

Steps:
1. Single-click a non-empty HELM cell so it becomes the current cell.
2. Open the Datagrok Context Pane (right side); locate the `Properties`
   panel and expand it.
3. Read the formula, molecular weight, and extinction coefficient from the
   panel's key/value table
   (`[name="pane-Properties"] table[data-source="Helm:Properties"]`).
4. Click a different HELM cell; confirm the Properties panel updates to the
   new sequence's values.
5. Synthesize a >1000-char HELM string (e.g. by chaining
   `r(A)p.r(C)p.r(G)p.r(U)p` ~60 times into an `RNA1{…}$$$$V2.0` wrapper)
   and assign it to a single cell of a scratch Macromolecule HELM column
   via the JS API; reopen the Properties panel against that cell.

Expected:
- Properties panel shows three key/value rows: `formula` (e.g.
  `C101H140N23O31P`), `molecular weight`, and `extinction coefficient`.
- The panel updates whenever the current cell changes.
- For a >1000-char sequence, the panel renders a `Too long sequence`
  warning (`[data-source="Helm:Properties"]:not(table)`) instead of the
  key/value table — no UI freeze.

### Scenario 4 (bundled — convert_helm_to_molfile): batch HELM → molfile conversion across the column

Steps:
1. With `HELM.csv` open, call `Helm:getMolfiles` against the `HELM` column
   via the JS API:
   `await grok.functions.call('Helm:getMolfiles', { col: helmColumn })`.
2. Wait for the call to resolve and inspect the returned column.
3. Wait > 1 second after the call returns, then re-issue the same call to
   confirm the cached off-screen editor evicts and re-instantiates without
   error.

Expected:
- The returned object is a string column with one molfile entry per row;
  row count equals the source column's row count.
- Each non-null entry is an `HWE pseudo-molfile` (header line
  `HWE pseudo-molfile`). **Do NOT assert a V2000/V3000 header** — the new
  build returns a pseudo-molfile. Empty source rows yield empty / null
  entries.
- Re-issuing the call after the 1-second eviction window succeeds; no
  console error.
