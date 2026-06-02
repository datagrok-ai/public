---
feature: bio
sub_features_covered:
  - bio.viewers.diversity-search
  - bio.search.diversity
  - bio.search.diversity.top-menu
  - bio.detector
  - bio.rendering
target_layer: playwright
coverage_type: smoke
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - bio-diversity-search-spec.ts
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:30:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T06:15:00Z
    spec_runs:
      - spec: bio-diversity-search-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 70
        failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T07:00:00Z
    review_round: 1
    failure_keys: []
---

# Bio | Search | Diversity Search — positive-path diversity docking and viewer-reaction contract

Smoke scenario for atlas critical_path
`bio.cp.diversity-search` (`priority: p0`): the Sequence
Diversity Search top-menu docks the
`SequenceDiversityViewer`, the macromolecule detector
classifies the input column synchronously, and the docked
viewer surfaces a maximally diverse subset of sequences
drawn from the full distance matrix (or the random-sampled
subset above the 10k-row threshold).

This is the positive-path complement to
`empty-input-row-viewers.md`, which covers the `coverage_type:
edge` empty-input rejection contract for the same viewer on
its current-row branch. The two scenarios together exercise
the `bio.viewers.diversity-search` surface across both
healthy and degenerate inputs, mirroring the
`bio-similarity-search.md` / `empty-input-row-viewers.md`
pair on the Similarity branch.

## Setup

- Open `System.AppData/Bio/tests/filter_FASTA.csv` from the
  Files browser. The Macromolecule detector (atlas
  `bio.detector`) classifies the sequence column synchronously
  on open; the table view opens with the Macromolecule cell
  renderer (atlas `bio.rendering`) painting sequence cells.
- Verify the table has >= 5 rows so the diversity-subset
  selection is meaningful and visibly different from a
  trivial first-N selection. The dataset's row count is well
  below the 10k random-sampling threshold, so the viewer
  uses the full distance-matrix path.

## Scenarios

### Scenario 1: Diversity Search top-menu docks the diversity viewer

Steps:
1. On the menu ribbon, open **Bio** > **Search** > **Diversity
   Search** (atlas `bio.search.diversity.top-menu`,
   registered at `public/packages/Bio/src/package.ts#L1283`).
2. The Sequence Diversity Search viewer
   (`SequenceDiversityViewer`, atlas
   `bio.viewers.diversity-search`) docks in the active view.
3. Wait for the viewer to compute the full distance matrix
   over the input column and surface the diverse subset.

Expected:
- The Sequence Diversity Search viewer is present in the
  active view's viewer list (`isViewerPresent` /
  `findViewer` style assertion against the viewer-name
  `Sequence Diversity Search`).
- The viewer displays a non-empty diversity-subset panel —
  i.e. at least two row cards / list entries are surfaced.
- No error balloon appears (positive-path complement to the
  empty-input contract in `empty-input-row-viewers.md`).

### Scenario 2: Reopening the viewer with a fresh dataset re-runs diversity selection

Steps:
1. With the Diversity Search viewer docked from Scenario 1,
   close the viewer (right-click viewer header → Close, or
   equivalent dock close).
2. Open a second dataset (`System.AppData/Bio/tests/filter_HELM.csv`)
   from the Files browser. The Macromolecule detector
   classifies the HELM column synchronously on open.
3. On the menu ribbon, open **Bio** > **Search** > **Diversity
   Search** again against the fresh table view.

Expected:
- The viewer re-docks against the new table view and
  surfaces a diversity subset whose row identifiers are
  drawn from the new HELM dataset rather than the prior
  FASTA dataset (visibly distinct sequence cards / list
  entries, confirming the diversity computation re-ran on
  the new column rather than re-displaying the cached
  Scenario-1 subset).
- The Macromolecule detector classified the HELM column
  synchronously on open (atlas `bio.detector`) — the
  Diversity Search top-menu becomes invokable without
  requiring a manual semType assignment.
- No error balloon appears.

## Notes

- atlas entry derived from critical_paths[bio.cp.diversity-search]
  (atlas `priority: p0`, propagated to scenario frontmatter
  `coverage_type: smoke` per STEP E p0 -> smoke mapping).
- target_layer rationale: `playwright` — top-menu dispatch
  (**Bio | Search | Diversity Search**) plus viewer-docking
  DOM plus cross-dataset re-run reaction is the assertion
  surface; apitest cannot exercise the ribbon path nor
  inspect the docked viewer card list. Consistent with
  sibling `bio-similarity-search.md`,
  `empty-input-row-viewers.md`, `sequence-activity-cliffs.md`,
  `sequence-space.md`, all `target_layer: playwright`.
- Sibling apitest coverage:
  `public/packages/Bio/src/tests/similarity-diversity-tests.ts`
  exists for the underlying diversity-compute / distance-
  matrix API path, but does not cover the top-menu dispatch
  or the docked-viewer DOM. This scenario adds the missing
  UI-layer assertion surface for the Diversity branch
  (sibling `bio-similarity-search.md` adds it for the
  Similarity branch).
- This scenario covers 5 sub_features
  (`F-STRUCT-DENSITY-01` floor: 2;
  `F-STRUCT-INTERACTION-01` floor: 3 in a multi-sub_feature
  scenario — satisfied).
- Manual-only subset: none of the five covered sub_features
  appear in atlas `manual_only[]` (verified against atlas rev
  3 `manual_only[]` list — `bio.viewers.web-logo`,
  `bio.viewers.vd-regions`, `bio.rendering.column-header`,
  `bio.rendering.macromolecule-difference`, the demo entries,
  `bio.panels.{structure-3d, atomic-level, tooltip}`).
- Deferrals: none. The viewer-docking DOM and the
  cross-dataset re-run reaction are both observable via
  Playwright; no pixel-precision or non-deterministic
  surface is involved.
- The section has no Bio/grok-browser ref-doc verb-form H2
  matching the citation regex, so `## Notes` citation
  pointers reference atlas entries only — no
  `See: <path>#<heading>` citation form applies here.
- Related-bug context: no curated `bug-library/bio.yaml`
  entry has Diversity Search in its `affects[]` set as a
  primary regression risk (GROK-16111 covers the empty-input
  edge contract across Similarity / Diversity / Activity-
  Cliffs row-viewers, addressed in
  `empty-input-row-viewers.md`); `related_bugs: []`.

---
{
  "order": 19
}
