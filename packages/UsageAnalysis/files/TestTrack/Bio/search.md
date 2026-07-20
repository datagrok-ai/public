---
feature: bio
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [bio.cp.subsequence-search]
realizes: [bio.search.subsequence-search, bio.bio-substructure-filter]
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

# Bio | Search | Subsequence Search — filter panel search & reset

Checks the **Bio | Search | Subsequence Search** filter-panel widget:
searching for a subsequence narrows the table to matching rows, and
clicking **Reset Filter** restores all rows.

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

- The Macromolecule-column picker in the search dialog runs
  non-interactively here, since `filter_FASTA.csv` only has one
  Macromolecule column to pick from.
- GROK-16111 (a Similarity Search bug) doesn't apply to Subsequence
  Search.
