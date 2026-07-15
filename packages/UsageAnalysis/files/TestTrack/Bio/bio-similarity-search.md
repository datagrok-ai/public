---
feature: bio
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: [bio.cp.similarity-search]
realizes: [bio.search.similarity-search, bio.sequence-similarity-search]
produced_from: atlas-driven
related_bugs: []
realized_as:
  - bio-similarity-search-spec.ts
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T06:00:00Z
    failure_keys: []
    review_round: 1
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T05:30:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T23:59:34Z
    spec_runs:
      - spec: bio-similarity-search-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 63
        failure_keys: []
---

# Bio | Search | Similarity Search — KNN viewer docking and row-click refresh

Checks that the **Bio | Search | Similarity Search** top-menu docks
the Sequence Similarity Search viewer, and that clicking a different
row in the table re-queries the K-nearest neighbours shown in the
docked viewer.

This is the healthy-input complement to `empty-input-row-viewers.md`,
which covers the empty-input rejection case for the same viewer —
together they exercise Similarity Search across both healthy and
degenerate inputs.

## Setup

- Open `System.AppData/Bio/tests/filter_FASTA.csv` from the
  Files browser. The Macromolecule detector classifies the
  sequence column synchronously on open; the table view opens
  with the Macromolecule cell renderer painting sequence cells.
- Verify the table has ≥ 5 rows so KNN K=3..5 is meaningful and
  the click-a-row reaction can be observed without exhausting
  the dataset.

## Scenarios

### Scenario 1: Similarity Search top-menu docks the KNN viewer

Steps:
1. On the menu ribbon, open **Bio** > **Search** > **Similarity
   Search**.
2. The Sequence Similarity Search viewer
   (`SequenceSimilarityViewer`) docks in the active view.
3. Wait for the viewer to compute the initial sparse KNN matrix
   for the current row.

Expected:
- The Sequence Similarity Search viewer is present in the
  active view's viewer list (`isViewerPresent` /
  `findViewer` style assertion against the viewer-name
  `Sequence Similarity Search`).
- The viewer displays a non-empty K-nearest-neighbours panel —
  i.e. at least one row card / list entry beyond the current
  row reference is shown.
- No error balloon appears (positive-path complement to the
  empty-input contract in `empty-input-row-viewers.md`).

### Scenario 2: Clicking a row re-queries the KNN viewer

Steps:
1. With the Similarity Search viewer docked from Scenario 1,
   click a different (non-current) row in the underlying table
   grid.
2. The viewer should re-query KNN against the newly-current
   sequence.

Expected:
- The viewer's row-card / list panel updates to show neighbours
  of the newly-clicked row — visibly different cards or a
  changed top-result identifier compared to the Scenario 1
  state.
- The currentRow indicator in the table tracks the click (the
  assertion here is on the visible viewer reaction, not on
  `grid.dataFrame.currentRowIdx` itself — that's an
  apitest-layer check).
- No error balloon appears.

## Notes

- Sibling coverage: `similarity-diversity-tests.ts` covers the
  underlying KNN-compute / Similarity API path, but not the
  top-menu dispatch or the docked-viewer UI. This scenario adds
  the missing UI-layer coverage.
- The empty-input edge case for this viewer (GROK-16111) is
  covered separately in `empty-input-row-viewers.md`, not here.

---
{
  "order": 12
}
