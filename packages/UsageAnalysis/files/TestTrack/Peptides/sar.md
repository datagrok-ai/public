---
feature: peptides
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [sar-from-context-panel-launch-button]
realizes: [bio.analyze.sar]
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/peptides/sar.md
migration_date: 2026-05-28
source_text_fixes:
  - open-linked-dataset-singular
  - mcl-viewer-missing-space
  - logo-summary-table-titlecase
  - drop-redundant-peptides-analysis-view-parens
  - clarify-settings-button-on-sar-toolbar
  - reword-step9-table-viewer-content-reloads
  - clarify-mutation-cliffs-pairs-plural
  - tighten-non-empty-cell-deterministic
candidate_helpers:
  - peptides.launchSarFromContextPanel
  - peptides.openSarSettingsDialog
  - peptides.toggleMutationCliffsInvariantMapMode
  - peptides.clickFirstNonEmptyCell
  - peptides.verifyContextPanelHasMutationCliffPairsAndDistribution
unresolved_ambiguities:
  - non-empty-cell-selection-strategy
  - calculation-completion-marker-needed
  - distribution-panel-parameter-arbitrariness
scope_reductions:
  - id: SR-01
    check: A-CONT-01
    rationale: |
      Step 7 "Change all parameters arbitrarily" cannot translate to a
      deterministic Playwright assertion without per-parameter
      expected-state declaration. Defer exhaustive per-parameter
      verification to a future parameterized peptides-settings-dialog
      spec; in this scenario, exercise the OK roundtrip with a single
      representative parameter change (Similarity threshold) to
      establish the settings → viewer-reload contract.
    verdict_status: SCOPE_REDUCTION
related_bugs:
  - GROK-19145
  - GROK-14357
  - GROK-18058
  - github-1549
  - GROK-14461
  - GROK-15934
  - GROK-14298
realized_as:
  - sar-spec.ts
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-05-28-peptides-migrate-01
    timestamp: 2026-05-28T00:00:00Z
    review_round: 1
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-30-peptides-automate-02
    timestamp: 2026-05-30T23:58:00Z
    spec_runs:
      - spec: sar-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 142
        failure_keys: []
  d:
    verdict: EVIDENCE_GAP
    cycle_id: 2026-05-29-peptides-migrate-02
    timestamp: 2026-05-29T18:05:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-05-30-peptides-automate-02
    timestamp: 2026-05-30T00:00:00Z
    failure_keys: []
  f:
    verdict: PASS
    cycle_id: 2026-05-29-peptides-migrate-02
    timestamp: 2026-05-29T00:00:00Z
    failure_keys: []
---

# SAR — Launch and verify viewers (context-panel entry path)

Launching SAR from the Peptides context panel's **Launch SAR** button should attach four viewers (Sequence Variability Map, Most Potent Residues, MCL Viewer, Logo Summary Table). This scenario also checks that a Settings dialog change reloads the viewers, that toggling between Mutation Cliffs and Invariant Map mode and clicking a cell populates the context panel, and that changing Distribution panel parameters re-renders it.

## Setup

- **Dataset:** `System:DemoFiles/bio/peptides.csv` (declared in the
  scenario's trailing JSON metadata; loaded via `Open the linked
  dataset`).
- **Initial table state:** an open Datagrok view holding the peptides
  table with a Macromolecule sequence column and a numerical activity
  column (this is what the atlas validator
  `peptides.workflow.sar-dialog` / `peptidesDialog` gates on).

## Scenarios

### Scenario 1 — Launch SAR from the Peptides context panel, verify four viewers render

This exercises the alternative SAR entry path (context-panel button)
— sister of the top-menu `Bio | Analyze | SAR...` path covered by
`peptide-space.md`. It also doubles as a regression target for
GROK-17557, an init-prerequisite race on this entry path.

1. Open the linked dataset (`System:DemoFiles/bio/peptides.csv`).
2. Click the peptides column title (the Macromolecule sequence
   column header).
3. On the **Context Panel**, expand the **Peptides** panel.
4. Click the **Launch SAR** button. Accept default config in the
   `Analyze Peptides` dialog and click **OK**.
5. Verify that four viewers appear and attach to the active
   TableView:
   - **Mutation Cliffs / Invariant Map** (the
     `Sequence Variability Map` / `MonomerPosition` viewer rendered
     in its default mode).
   - **Most Potent Residues** (`MostPotentResidues` viewer).
   - **MCL Viewer** (may take some time — MCL clustering runs in a
     Web Worker; await completion).
   - **Logo Summary Table** viewer.

### Scenario 2 — Apply a settings change via the SAR Settings dialog and verify viewers reload

This exercises the SAR Settings dialog flow and its propagation to
the attached viewers and the underlying model cache. Exhaustive
per-parameter verification is deferred to a future parameterized
settings-dialog spec; this scenario asserts only the settings →
viewer-reload contract via one representative parameter change.

6. Click the **Settings** button on the SAR analysis toolbar (opens
   the `getSettingsDialog(model)` accordion).
7. Change one representative parameter (e.g. Similarity threshold)
   to a non-default value. _(Exhaustive per-parameter verification is
   deferred — see the spec's per-parameter cases for the full matrix.)_
8. Click **OK**.
9. Verify that the table and viewer content reloads according to the
   applied parameters (the WebLogo column-header rendering, the
   MonomerPosition cells, and the Most Potent Residues grid all
   re-render against the updated `MonomerPositionStats` cache).

### Scenario 3 — Toggle Mutation Cliffs / Invariant Map mode and click a cell to populate Context Panel

This exercises the SAR viewer's Mutation Cliffs / Invariant Map
mode-switching and the cell-click → context-panel-population path.

10. Switch between **Mutation Cliffs** and **Invariant Map** modes
    on the `Sequence Variability Map` viewer (the
    `MonomerPosition` viewer's mode toggle).
11. Click the first non-empty cell at a deterministic position
    (e.g. row index 0 in the visible grid) in **Mutation Cliffs** or
    **Invariant Map**. _(Original wording was "random non-empty
    cell"; a deterministic position is used instead for repeatable
    automation.)_
12. Verify that the **Context Panel** displays the necessary panels:
    **Mutation Cliffs pairs** and **Distribution**.

### Scenario 4 — Adjust Distribution panel parameters and verify rendering

13. Change parameters on the **Distribution** panel to a single
    representative non-default value (a non-default breakdown
    grouping or aggregation). _(Exhaustive parameter-matrix
    verification is deferred; the original "arbitrarily" wording is
    resolved to one representative change until a parameterized
    variant lands.)_

## Notes

- **Entry path.** This scenario uses the Peptides context-panel **Launch SAR** button entry path. The sister entry path (top-menu `Bio | Analyze | SAR...`) is covered by `peptide-space.md`.
- **Related bugs.** This scenario's surface touches 7 curated Peptides bugs. Five (GROK-19145, GROK-14357, GROK-18058, github-1549, GROK-14461) each have a dedicated repro test elsewhere that this scenario complements. Two (GROK-15934, GROK-14298) are listed for awareness but aren't directly reproduced here, since this scenario doesn't include a hover step (GROK-15934) or a filter step (GROK-14298).
- **Related bug — GROK-17557.** The context-panel Launch SAR entry path is the regression target for this init-prerequisite race (viewers must not be requested before component loading completes). Step 4 implicitly exercises this; a dedicated spec covers the init-prerequisite assertion directly.
- **Deferral — calculation completion timing.** Step 5 ("four viewers should appear") and Step 9 (settings-driven reload) both depend on an asynchronous Web Worker compute (MCL/DBSCAN); the exact wait condition is left to the automation rather than pinned down here.

## Original trailing metadata (preserved verbatim from source)

```json
{
  "order": 3,
  "datasets": [
    "System:DemoFiles/bio/peptides.csv"
  ]
}
```
