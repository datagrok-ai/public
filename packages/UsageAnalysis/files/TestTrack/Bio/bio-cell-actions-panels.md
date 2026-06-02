---
feature: bio
sub_features_covered:
  - bio.actions.copy-as
  - bio.editors.get-region
  - bio.editors.split-to-monomers
  - bio.panels.composition-analysis
  - bio.panels.monomer-info
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - bio-cell-actions-panels-spec.ts
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T10:15:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T09:51:00Z
    spec_runs:
      - spec: bio-cell-actions-panels-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 123
        failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T11:00:00Z
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

# Bio | Macromolecule cell context actions and per-cell info panels

Regression scenario for the per-cell interaction surface of a
Macromolecule column — the right-click "Copy as ..." action set
(`addCopyMenu`, atlas `bio.actions.copy-as`), the two custom
editors that drive the Bio | Calculate | Get Region... and Bio |
Transform | Split to Monomers... top-menus (atlas
`bio.editors.get-region`, `bio.editors.split-to-monomers`), and
the two info-panel widgets that appear on the Datagrok Context
Pane when a Macromolecule cell is the current selection
(`compositionAnalysisWidget`, atlas
`bio.panels.composition-analysis`; `monomerInfoPanel`, atlas
`bio.panels.monomer-info`).

This is a breadth-coverage scenario authored by Gate F
coverage-extension. It complements the existing top-menu /
viewer-docking scenarios in the section (`analyze.md`,
`search.md`, `convert.md`, `bio-similarity-search.md`,
`bio-diversity-search.md`) by covering the cell-level and
context-panel surfaces that the umbrella scenarios do not
exercise.

## Setup

- Open `System.AppData/Bio/tests/filter_FASTA.csv` from the
  Files browser. The Macromolecule detector classifies the
  sequence column synchronously on open; the table view opens
  with the Macromolecule cell renderer painting cells.
- Verify the table has ≥ 5 rows so the cell-level operations
  below have non-degenerate input.

## Scenarios

### Scenario 1: Right-click Copy as ... entries on a Macromolecule cell

Steps:
1. Right-click any non-empty cell in the Macromolecule column.
   The cell-level context menu opens (atlas
   `bio.actions.copy-as`, registered in
   `public/packages/Bio/src/package.ts#L1527` via
   `addCopyMenu(cell, menu)`).
2. Verify the menu contains four "Copy as ..." entries:
   `Copy as FASTA`, `Copy as SEPARATOR`, `Copy as HELM`,
   `Copy as BILN`.
3. Click `Copy as FASTA`.

Expected:
- The menu shows all four notation conversion targets, even
  though the source cell is in FASTA notation — the entries
  must be unconditionally present for round-trip and
  cross-notation copy.
- After clicking `Copy as FASTA`, the system clipboard contains
  the FASTA-form representation of the clicked cell (assertion
  is on the action being dispatched without error; clipboard
  content read-back is best-effort given Playwright clipboard
  permissions).
- No error balloon appears.

### Scenario 2: Composition analysis and Monomer info panels on a Macromolecule cell

Steps:
1. Single-click any non-empty cell in the Macromolecule column
   so it becomes the current cell.
2. Open the Datagrok Context Pane (the right-side panel; helper
   `waitForElement` against the context-pane container).
3. Locate the `Composition analysis` panel section (atlas
   `bio.panels.composition-analysis`, registered in
   `public/packages/Bio/src/package.ts#L403` via
   `compositionAnalysisWidget`). Expand if collapsed.
4. Locate the `Monomer` panel section (atlas
   `bio.panels.monomer-info`, registered in
   `public/packages/Bio/src/package.ts#L412` via
   `monomerInfoPanel`). Expand if collapsed.

Expected:
- The Composition analysis panel renders a monomer-composition
  table with color-coded counts for the selected sequence cell
  (cell-scope is the single current cell, not the whole column).
- The Monomer panel renders monomer info for the current cell
  context (semType `Monomer` for individual monomer cells, or
  the cell-level summary for the parent Macromolecule cell —
  the panel must surface non-empty content).
- No error balloon appears in either panel.

### Scenario 3: Get Region editor dialog opens for the active sequence column

Steps:
1. On the menu ribbon, open **Bio** > **Calculate** > **Get
   Region** (atlas `bio.calculate.get-region.top-menu` —
   covered by `bio-calculate-scoring.md` for the dispatch path;
   this step opens the editor surface).
2. The `GetRegionEditor` custom editor dialog opens (atlas
   `bio.editors.get-region`, registered in
   `public/packages/Bio/src/package.ts#L213`).
3. Inspect the editor for the sequence-column selector input
   bound to the active Macromolecule column.

Expected:
- The editor dialog is visible and addressable
  (`waitForElement` against the dialog container).
- The sequence-column selector input is populated with the
  current view's Macromolecule columns — at least the FASTA
  column from setup is present and selectable.
- Cancelling the dialog (Escape or `Cancel` button) closes it
  without firing a balloon error.

### Scenario 4: Split to Monomers editor dialog opens for the active sequence column

Steps:
1. On the menu ribbon, open **Bio** > **Transform** > **Split
   to Monomers...**. The `SplitToMonomersEditor` custom editor
   dialog opens (atlas `bio.editors.split-to-monomers`,
   registered in
   `public/packages/Bio/src/package.ts#L226`).
2. Inspect the editor for the sequence-column selector input.

Expected:
- The editor dialog is visible and addressable.
- The sequence-column selector input is populated with the
  current view's Macromolecule columns.
- Cancelling the dialog closes it without firing a balloon
  error.

## Notes

- Scenario produced by Test Designer Gate F coverage-extension
  (cycle 2026-06-01-bio-migrate-02). Net-new sub_features against
  the executor-derived live covered union (68 ids): 5
  (`bio.actions.copy-as`, `bio.editors.get-region`,
  `bio.editors.split-to-monomers`,
  `bio.panels.composition-analysis`,
  `bio.panels.monomer-info`). Projected effect on
  F-STRUCT-COVERAGE-01: union grows 68 → 73; non-manual_only
  covered 66 → 71; ratio 66.7% → ~71.7% — crosses the 70%
  threshold per the gap_description projection.
- target_layer rationale: `playwright` — the verification surface
  is cell-context-menu DOM dispatch (Pattern 5 context-menu
  shape), Context Pane info-panel widget rendering, and custom
  editor dialog DOM. apitest cannot exercise the right-click cell
  menu nor the Context Pane widget paint; the bound editor
  surfaces are dialog-only. Consistent with sibling
  `bio-renderer-dispatch.md`, `bio-manage-libraries-crud.md`,
  `analyze.md`, all `target_layer: playwright`.
- coverage_type rationale: `regression` — none of the five
  covered ids appear in any atlas `critical_paths[]` entry as a
  smoke-eligible golden path; none map to an atlas `edge_cases[]`
  boundary entry; none are stress / large-data shaped. General
  feature-shape coverage of the cell-level + info-panel surfaces
  → `regression` per STEP E heuristic.
- Manual-only subset: none of the five covered ids appear in
  atlas `manual_only[]` (verified against atlas rev 3 list —
  `bio.viewers.web-logo`, `bio.viewers.vd-regions`,
  `bio.rendering.column-header`,
  `bio.rendering.macromolecule-difference`, the demo entries,
  `bio.panels.{structure-3d, atomic-level, tooltip}`;
  `bio.panels.composition-analysis` and `bio.panels.monomer-info`
  are NOT on that list).
- This scenario carries 5 sub_features and contains 4 scenarios
  with `sub_features_used` spread across all four; cardinality
  satisfies `F-STRUCT-DENSITY-01` floor: 2 and
  `F-STRUCT-INTERACTION-01` floor: 3.
- Sibling sub_feature pairs: the related top-menu paths
  (`bio.calculate.get-region.top-menu`,
  `bio.transform.split-to-monomers`) are already covered by
  `bio-calculate-scoring.md` and the Transform top-menu surface
  in `convert.md` / `bio-transform-atomic-level.md`. This
  scenario adds the editor-dialog and per-cell-action surfaces
  that those top-menu scenarios do not exercise.
- The section has no Bio/grok-browser ref-doc verb-form H2
  matching the F-UI-COVERAGE-01 citation regex, so `## Notes`
  citations reference atlas source-line anchors only; no
  `See: <path>#<heading>` citation form applies.
- Related-bug context: no curated `bug-library/bio.yaml` entry
  has any of the five covered ids in its `affects[]` set
  (GROK-12164 covers `bio.rendering.separator` etc.; GROK-16596
  covers `bio.rendering.macromolecule-difference`; the cell-level
  Copy actions and the per-cell info panels are not currently
  bug-tracked); `related_bugs: []`.
- Deferrals: none. Cell context-menu DOM, info-panel widget DOM
  on the Context Pane, and custom editor dialog DOM are all
  observable via Playwright; no pixel-precision or
  non-deterministic surface is involved. The Composition
  analysis panel uses color-coded counts but the assertion is on
  non-empty table content, not pixel color.

---
{
  "order": 25
}
