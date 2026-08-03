---
feature: bio
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: [bio.cp.diversity-search]
realizes: [bio.search.diversity-search, bio.sequence-diversity-search]
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

# Bio | Search | Diversity Search — diverse-subset viewer docking and refresh

Checks that the **Bio | Search | Diversity Search** top-menu docks the
Sequence Diversity Search viewer, and that it surfaces a maximally
diverse subset of sequences (drawn from the full distance matrix, or
a random sample once the dataset is above 10k rows).

This is the healthy-input complement to `empty-input-row-viewers.md`,
which covers the empty-input rejection case for the same viewer.
Together the two scenarios exercise Diversity Search across both
healthy and degenerate inputs, mirroring the `bio-similarity-search.md`
/ `empty-input-row-viewers.md` pair for Similarity Search.

## Setup

- Open `System.AppData/Bio/tests/filter_FASTA.csv` from the
  Files browser. The Macromolecule detector classifies the
  sequence column synchronously on open; the table view opens
  with the Macromolecule cell renderer painting sequence cells.
- Verify the table has >= 5 rows so the diversity-subset
  selection is meaningful and visibly different from a
  trivial first-N selection. The dataset's row count is well
  below the 10k random-sampling threshold, so the viewer
  uses the full distance-matrix path.

## Scenarios

### Scenario 1: Diversity Search top-menu docks the diversity viewer

Steps:
1. On the menu ribbon, open **Bio** > **Search** > **Diversity
   Search**.
2. The Sequence Diversity Search viewer
   (`SequenceDiversityViewer`) docks in the active view.
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
  synchronously on open — the Diversity Search top-menu
  becomes invokable without requiring a manual semType
  assignment.
- No error balloon appears.

## Notes

- Sibling coverage: `public/packages/Bio/src/tests/similarity-diversity-tests.ts`
  covers the underlying diversity-compute / distance-matrix logic at
  the API level, but not the top-menu dispatch or the docked-viewer
  UI. This scenario adds the missing UI-layer coverage for the
  Diversity branch (the sibling `bio-similarity-search.md` does the
  same for Similarity Search).
- The empty-input edge case for this viewer (GROK-16111) is covered
  separately in `empty-input-row-viewers.md`, not here.

---
{
  "order": 19
}
