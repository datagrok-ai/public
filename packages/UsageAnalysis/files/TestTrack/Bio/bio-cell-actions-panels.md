---
feature: bio
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: []
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
---

# Bio | Macromolecule cell context actions and per-cell info panels

Covers the per-cell interaction surface of a Macromolecule column:
the right-click "Copy as ..." menu (FASTA / SEPARATOR / HELM / BILN),
the custom editor dialogs behind **Bio | Calculate | Get Region...**
and **Bio | Transform | Split to Monomers...**, and the two
info-panel widgets — Composition analysis and Monomer info — that
appear on the Context Panel when a Macromolecule cell is selected.

This complements the top-menu and viewer-docking scenarios elsewhere
in the section (`analyze.md`, `search.md`, `convert.md`,
`bio-similarity-search.md`, `bio-diversity-search.md`) by covering
the cell-level and context-panel surfaces those umbrella scenarios
don't exercise.

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
   The cell-level context menu opens.
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
3. Locate the `Composition analysis` panel section. Expand if
   collapsed.
4. Locate the `Monomer` panel section. Expand if collapsed.

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
   Region** (covered by `bio-calculate-scoring.md` for the
   dispatch path; this step opens the editor surface).
2. The `GetRegionEditor` custom editor dialog opens.
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
   dialog opens.
2. Inspect the editor for the sequence-column selector input.

Expected:
- The editor dialog is visible and addressable.
- The sequence-column selector input is populated with the
  current view's Macromolecule columns.
- Cancelling the dialog closes it without firing a balloon
  error.

## Notes

- Sibling coverage: the underlying top-menu paths for Get Region and
  Split to Monomers are already covered by `bio-calculate-scoring.md`
  and by the Transform top-menu surface in `convert.md` /
  `bio-transform-atomic-level.md`. This scenario adds the
  editor-dialog and per-cell-action surfaces (right-click Copy as,
  per-cell info panels) that those top-menu scenarios don't exercise.
- The Composition analysis panel renders color-coded monomer counts,
  but this scenario only asserts the table has non-empty content —
  it does not check pixel-level color correctness.

---
{
  "order": 25
}
