---
feature: bio
target_layer: playwright
coverage_type: edge
priority: p1
realizes_atlas: [bio.x.empty-input-on-row-viewers]
realizes: [bio.sequence-similarity-search, bio.sequence-diversity-search, bio.analyze.activity-cliffs]
produced_from: atlas-driven
related_bugs:
  - GROK-16111
realized_as:
  - empty-input-row-viewers-spec.ts
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions:
  - id: SR-01
    check: E-SCENARIO-RUNTIME-ALIGNMENT
    rationale: |
      empty-input rejection balloon deferred: GROK-16111 — the three Bio current-row viewers (Sequence Similarity Search, Sequence Diversity Search, Activity Cliffs) do not reject empty/null current-row input; they silently KNN on the empty cell and surface NO rejection balloon (GROK-16111 status: regression-risk, fixed_in: ''). The scenario assertion (Scenario "Expected" Invariant 1: balloon invocation count > 0 on empty input) is correct; the product is broken. Gate B re-fired deterministically [B-RUN-PASS, B-STAB-01] across cycle 2026-06-01-bio-migrate-02 because the hard expect(probe.balloonCount).toBeGreaterThan(0) caught the real, unfixed bug. Per operator Option A (2026-06-02) the proper routing for an assertion on a documented, unfixed product bug is scope reduction (mirrors PowerPack data-enrichment SR-05..08, GROK-20175): replaced the hard expect() with a guarded console.warn (no assertion) so the run is green while GROK-16111 is open. The live no-crash assertion (Invariant 2: the active dataframe row count is unchanged — no silent zero-row result / silent table rewrite) is RETAINED, not softened. Revert this SR + restore the hard balloon assertion when GROK-16111 fixed_in is set.
    verdict_status: SCOPE_REDUCTION
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T05:00:00Z
    failure_keys: []
    review_round: 1
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-02-bio-automate-01
    timestamp: 2026-06-02T07:30:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-02-bio-automate-01
    timestamp: 2026-06-02T16:16:00Z
    spec_runs:
      - spec: empty-input-row-viewers-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 404
        failure_keys: []
---

# Bio — empty current-row input across Similarity Search, Diversity Search & Activity Cliffs

Checks that the three Bio viewers that operate on the current table
row — Sequence Similarity Search, Sequence Diversity Search, and
Activity Cliffs — reject an empty/null input sequence with a
user-visible error balloon, and never silently produce an empty
result.

This is intentionally narrow: it only asserts the error-balloon
behavior and the absence of a silent zero-row result, not the
embedding/KNN math itself — that's covered by the sibling
positive-path scenarios (`sequence-activity-cliffs.md`,
`bio-similarity-search.md`, `bio-diversity-search.md`).

## Setup

- Open `System.AppData/Bio/tests/filter_FASTA.csv` from the
  Files browser. The Macromolecule detector classifies the
  sequence column synchronously; the table view opens.
- Verify the table has ≥ 2 rows so the "empty current cell"
  state can be constructed without an empty-table degenerate.
- Position the **current row** on a row whose sequence cell
  has been cleared to empty/null (set the cell to empty via
  grid cell-edit → blank → commit, or pick a row that already
  has an empty Macromolecule cell after `Edit | Remove rows`
  on a placeholder-row dataset). The current-row sequence
  cell MUST be empty/null when each of the three viewers
  below is invoked.

## Scenarios

### Scenario 1: Sequence Similarity Search — empty current-row balloon

Steps:
1. With the current row pointing at the empty-sequence cell,
   on the menu ribbon open **Bio** > **Search** > **Similarity
   Search**.
2. The Sequence Similarity Search viewer
   (`SequenceSimilarityViewer`) docks in the active view.
3. Observe the viewer's reaction to the empty current-row
   input.

Expected:
- A user-visible error balloon surfaces, naming the
  empty/null input as the reason the viewer cannot compute
  K-nearest-neighbours (this is a regression-risk area per
  GROK-16111).
- The viewer does NOT silently produce a zero-row result
  table — either the viewer dock content remains empty with
  a placeholder ("no input" / equivalent) OR the viewer
  refuses to dock until a non-empty current row is present.
  Either acceptance shape is contract-compliant; a silent
  empty-rows result is the failure mode.

### Scenario 2: Sequence Diversity Search — empty current-row balloon

Steps:
1. Restore the current row to the empty-sequence cell (if
   Scenario 1 moved it).
2. On the menu ribbon open **Bio** > **Search** > **Diversity
   Search**.
3. The Sequence Diversity Search viewer
   (`SequenceDiversityViewer`) docks.
4. Observe the viewer's reaction.

Expected:
- Error balloon surfaces (same contract as Scenario 1) —
  empty/null current-row input is rejected, not silently
  consumed.
- No zero-row diversity-result table is produced.

### Scenario 3: Activity Cliffs — empty current-row balloon

Steps:
1. Restore the current row to the empty-sequence cell.
2. On the menu ribbon open **Bio** > **Analyze** > **Activity
   Cliffs...**
3. In the `SeqActivityCliffsEditor` dialog, leave defaults and
   click **OK**.
4. Observe the engine's reaction.

Expected:
- Error balloon surfaces (same contract): the empty-sequence
  current row makes embedding / cliff detection ill-defined;
  the engine must surface the rejection via balloon rather
  than complete with an empty cliff overlay.
- No silent empty ScatterPlot-with-cliffs viewer dock is
  produced.

---
{
  "order": 11
}
