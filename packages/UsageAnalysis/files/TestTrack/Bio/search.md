---
feature: bio
sub_features_covered:
  - bio.search.subsequence
  - bio.search.subsequence.top-menu
  - bio.search.subsequence.editor
  - bio.search.subsequence.filter
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/bio/search.md
migration_date: 2026-06-01
source_text_fixes:
  - promote-step-1-dataset-open-to-setup-section
  - split-step-3-action-and-verification
  - split-step-4-action-and-verification
  - join-step-3-wrapped-sequence-literal
  - convert-trailing-order-json-to-frontmatter
  - bold-input-field-name-sequence
  - bold-action-reset-filter
candidate_helpers: []
unresolved_ambiguities:
  - editor-dialog-skip-applies-only-to-single-macromolecule-column-datasets
scope_reductions: []
related_bugs: []
realized_as:
  - search-spec.ts
order: 3
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T13:30:00Z
    review_round: 1
    failure_keys: []
  d:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T13:00:00Z
    failure_keys: []
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-02-bio-automate-01
    timestamp: 2026-06-02T12:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-02-bio-automate-01
    timestamp: 2026-06-02T12:22:30Z
    spec_runs:
      - spec: search-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 43
        failure_keys: []
---

## Setup

1. Open the dataset `System.AppData/Bio/tests/filter_FASTA.csv`.

## Scenarios

### Scenario 1 — Subsequence search on the filter panel + Reset Filter restores all rows

1. On the menu ribbon, open **Bio** > **Search** > **Subsequence Search**.
2. On the filter panel, set **Sequence** to
   `RTDEVSNHTHDKPTLTWFEEIFEEYHSP`.
3. Verify there is **1 row**.
4. Click **Reset Filter**.
5. Verify all rows should be present.

## Notes

- Atlas critical path `bio.cp.subsequence-search` (priority p0)
  covers this end-to-end flow: open sequences dataset → Subsequence
  Search top-menu → filter widget docks into Filters panel → filter
  rows by subsequence query → Reset Filter restores all rows.
- Atlas sub-features covered:
  `bio.search.subsequence` (umbrella),
  `bio.search.subsequence.top-menu`,
  `bio.search.subsequence.editor` (Macromolecule-column picker —
  runs non-interactively on single-Macromolecule-column datasets like
  `filter_FASTA.csv`; this is the `unresolved_ambiguities` entry),
  `bio.search.subsequence.filter` (`bioSubstructureFilter` widget).
- Helper coverage: `bio.flow.subsequence-search` in helpers-registry
  (see `.claude/skills/grok-browser/references/bio.md:202`)
  documents the selector surface and the `14 → 1 → 14 trueCount`
  invariant.
- Chain context (`scenario-chains/bio.yaml`): `pyramid_layer:
  integration` (multi-subsystem path top-menu → filter widget
  docking → filter-panel input → grid-filter callback). The
  ui-smoke slot is held by `manage.md` (Rule 1 cardinality), so
  this scenario carries `coverage_type: regression`.
- `related_bugs: []` — no bug in `bug-library/bio.yaml` intersects
  `sub_features_covered`. GROK-16111 (`bio.search.similarity.*`)
  is a sibling-search-mode bug that does not apply to Subsequence
  Search.
