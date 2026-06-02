---
feature: bio
sub_features_covered:
  - bio.viewers.similarity-search
  - bio.search.similarity
  - bio.search.similarity.top-menu
  - bio.detector
  - bio.rendering
target_layer: playwright
coverage_type: smoke
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

# Bio | Search | Similarity Search — positive-path KNN docking and row-click contract

Smoke scenario for atlas critical_path
`bio.cp.similarity-search` (`priority: p0`): the Sequence
Similarity Search top-menu docks the
`SequenceSimilarityViewer`, the macromolecule detector
classifies the input column synchronously, and clicking a row
in the table re-queries the K-nearest neighbours displayed in
the docked viewer.

This is the positive-path complement to
`empty-input-row-viewers.md`, which covers the `coverage_type:
edge` empty-input rejection contract for the same viewer. The
two scenarios together exercise the
`bio.viewers.similarity-search` surface across both healthy and
degenerate inputs.

## Setup

- Open `System.AppData/Bio/tests/filter_FASTA.csv` from the
  Files browser. The Macromolecule detector (atlas
  `bio.detector`) classifies the sequence column synchronously
  on open; the table view opens with the Macromolecule cell
  renderer (atlas `bio.rendering`) painting sequence cells.
- Verify the table has ≥ 5 rows so KNN K=3..5 is meaningful and
  the click-a-row reaction can be observed without exhausting
  the dataset.

## Scenarios

### Scenario 1: Similarity Search top-menu docks the KNN viewer

Steps:
1. On the menu ribbon, open **Bio** > **Search** > **Similarity
   Search** (atlas `bio.search.similarity.top-menu`, registered
   at `public/packages/Bio/src/package.ts#L1261`).
2. The Sequence Similarity Search viewer
   (`SequenceSimilarityViewer`, atlas
   `bio.viewers.similarity-search`) docks in the active view.
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
- The currentRow indicator in the table tracks the click (atlas
  `bio.search.similarity` viewer-current-row bridge — the
  scenario level assertion is on the visible viewer reaction,
  not on `grid.dataFrame.currentRowIdx` itself; that is
  apitest-layer).
- No error balloon appears.

## Notes

- atlas entry derived from critical_paths[bio.cp.similarity-search]
  (atlas `priority: p0`, propagated to scenario frontmatter
  `coverage_type: smoke` per STEP E p0 → smoke mapping).
- target_layer rationale: `playwright` — top-menu dispatch
  (**Bio | Search | Similarity Search**) plus viewer-docking
  DOM plus row-click re-query bridge is the assertion surface;
  apitest cannot exercise the ribbon path nor inspect the
  docked viewer card list. Consistent with sibling
  `empty-input-row-viewers.md`, `sequence-activity-cliffs.md`,
  `sequence-space.md`, all `target_layer: playwright`.
- Sibling apitest coverage:
  `public/packages/Bio/src/tests/similarity-diversity-tests.ts`
  exists for the underlying KNN compute / Similarity API path,
  but does not cover the top-menu dispatch or the docked-viewer
  DOM. This scenario adds the missing UI-layer assertion
  surface.
- This scenario covers 5 sub_features
  (`F-STRUCT-DENSITY-01` floor: 2;
  `F-STRUCT-INTERACTION-01` floor: 3 in a multi-sub_feature
  scenario — satisfied).
- Manual-only subset: none of the five covered sub_features
  appear in atlas `manual_only[]` (verified against atlas rev 3
  `manual_only[]` list — `bio.viewers.web-logo`,
  `bio.viewers.vd-regions`, `bio.rendering.column-header`,
  `bio.rendering.macromolecule-difference`, the demo entries,
  `bio.panels.{structure-3d, atomic-level, tooltip}`).
- Deferrals: none. The viewer-docking DOM and the row-click
  re-query reaction are both observable via Playwright; no
  pixel-precision or non-deterministic surface is involved.
- The section has no Bio/grok-browser ref-doc verb-form H2
  matching the citation regex, so `## Notes` citation pointers
  reference atlas entries only — no `See: <path>#<heading>`
  citation form applies here.
- Related-bug context: no curated `bug-library/bio.yaml` entry
  has Similarity Search in its `affects[]` set
  (GROK-16111 covers the empty-input edge contract, addressed
  in `empty-input-row-viewers.md`); `related_bugs: []`.

---
{
  "order": 12
}
