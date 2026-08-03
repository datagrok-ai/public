---
feature: peptides
sub_features_covered:
  - peptides.workflow.sar-dialog
  - peptides.workflow.analyze-ui
  - peptides.workflow.start-analysis
  - peptides.viewers.monomer-position
  - peptides.viewers.most-potent-residues
  - peptides.viewers.mutation-cliffs
  - peptides.viewers.logo-summary-table
  - peptides.widgets.settings-dialog
  - peptides.widgets.distribution
  - peptides.widgets.mutation-cliffs
  - peptides.panels.peptides
target_layer: playwright
coverage_type: regression
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
    evidence_needed: |
      Independent Gate D re-review (fresh invocation; no prior verdict
      inherited -- the EVIDENCE_GAP already on disk and the matching
      decision-log mirrors were NOT consulted to form this verdict; I
      reached the same conclusion by my own diff attempt and re-verified
      every filesystem fact below this invocation).
      Mechanical / fail-fast checks re-derived against the on-disk
      MIGRATED artifact -- all PASS, so I did NOT fail-fast and proceeded
      to the content checks:
      - D-STRUCT-MECH-03: all 8 required migration fields present
        (feature=peptides; sub_features_covered, 11 entries; target_layer=
        playwright; coverage_type=regression; produced_from=migrated, in
        enum; original_path; migration_date=2026-05-28; related_bugs, 7
        entries). Deprecated migrated_from absent.
      - D-STRUCT-MECH-05: original_path resolves to a real on-disk file,
        so the mechanical "exists" predicate holds.
      - D-FRONTMATTER-PHASE1-01: all four Phase 1 fields present as
        parseable YAML lists.
      - D-FRONTMATTER-PHASE1-02: source_text_fixes (8) and
        unresolved_ambiguities (3) are kebab-case slugs <= 80 chars, no
        duplicates; candidate_helpers (5) are dotted-path names with no
        (args) suffix and no duplicates; scope_reductions holds one
        well-formed entry SR-01 {id, check=A-CONT-01, rationale,
        verdict_status=SCOPE_REDUCTION}.
      The content checks (D-STEP-01, D-STEP-02, D-EDGE-01, D-STRUCT-01,
      D-STRUCT-02, D-SAN-02, D-MERIT-01, D-MERIT-02, D-UI-DELEGATION-01)
      each require a diff against a DISTINCT raw TestTrack original, which
      is not reachable this invocation:
      - The dispatched inputs block carries only scenario_path and
        mode_file; no distinct original_path / original input was supplied.
      - A Glob over TestTrack/[Pp]eptides surfaced only this one sar.md
        (plus sar-run.md, sar-spec.ts, sar.md.migrator.dispatch.yaml); no
        sar-original.md, no *-migration-report.md, no .orig sibling, no
        pre-migration variant is reachable. On the case-insensitive
        Windows filesystem the frontmatter original_path
        (.../peptides/sar.md) and the on-disk capitalized
        .../Peptides/sar.md are the SAME file -- the carrier of THIS
        migrated frontmatter (it holds gate_verdicts.a/d/f). So
        original_path self-references the MIGRATED artifact, not a distinct
        pre-migration source.
      - The decision-log slice for feature=peptides carries multiple
        gate: D / verdict: EVIDENCE_GAP / failure_keys: [] recon mirrors
        for sar.md (and sibling peptides.md), all citing the same
        "no distinct raw original to diff against" tool-reachability
        condition; the failed_attempts[] block holds only legend-feature
        entries, none for peptides, so nothing rules out any verdict. The
        slice corroborates a recognized section-wide condition rather than
        a scenario-specific defect, and supplies no distinct original.
      - The migrator dispatch sidecar (sar.md.migrator.dispatch.yaml)
        records body_write: false / re_migration: true and states the only
        original text preserved on disk is the body's "## Original trailing
        metadata" JSON block (order + datasets). That JSON is NOT the
        original's numbered steps / Expected-result assertions, so it
        cannot serve as the content-diff baseline either.
      Under Gate D independence the migrated body's self-described
      provenance (the source_text_fixes slugs, the preserved trailing JSON
      datasets block, the determinism gaps in unresolved_ambiguities,
      SR-01's settings -> viewer-reload reduction) is the producer's own
      account and, by the gate's anti-self-justification design, cannot
      substitute for an independent diff. This is a tool-reachability gap
      on the ORIGINAL, NOT a confirmed content drop, so the verdict is
      EVIDENCE_GAP rather than FAIL.
      The decision-log shows this scenario has already exhausted its review
      rounds at Gate D within cycle 2026-05-29-peptides-migrate-02 with no
      PASS / SCOPE_REDUCTION. Per migration-mode.md, that routes to the
      autopilot-boundary escalation (trigger #1); re-invoking Gate D
      against the same self-referential original will only loop. The
      orchestrator must take one of:
      1. Retrieve the git-tracked pre-migration baseline blob of
         TestTrack/Peptides/sar.md (the revision prior to the migration
         commit) or an operator-supplied raw original, repoint
         original_path at that raw original (handoff invariant:
         original_path references the raw original, never the migrated
         artifact), and re-dispatch with a distinct original in hand so
         the content checklist can run; OR
      2. Route this scenario through the Batch 6.54 per_scenario_migrate
         Critic D bypass (atlas-driven / already-migrated scenarios skip
         Gate D) rather than re-invoking Gate D against a self-referential
         original.
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
`peptide-space.md`. Matches atlas critical_path
`sar-from-context-panel-launch-button` (p1, GROK-17557 init-
prerequisite regression target).

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

This exercises the SAR Settings dialog flow (atlas
`peptides.widgets.settings-dialog`) and its propagation to the
attached viewers and the underlying PeptidesModel cache. Note: per
`scope_reductions[SR-01]`, exhaustive per-parameter verification is
deferred to a future parameterized peptides-settings-dialog spec;
this scenario asserts only the settings → viewer-reload contract via
one representative parameter change.

6. Click the **Settings** button on the SAR analysis toolbar (opens
   the `getSettingsDialog(model)` accordion).
7. Change one representative parameter (e.g. Similarity threshold)
   to a non-default value. _(Per SR-01, exhaustive per-parameter
   verification is deferred — see the spec's per-parameter cases for
   the full matrix.)_
8. Click **OK**.
9. Verify that the table and viewer content reloads according to the
   applied parameters (the WebLogo column-header rendering, the
   MonomerPosition cells, and the Most Potent Residues grid all
   re-render against the updated `MonomerPositionStats` cache).

### Scenario 3 — Toggle Mutation Cliffs / Invariant Map mode and click a cell to populate Context Panel

This exercises the SARViewer mode-switching and the cell-click →
context-panel-population path (atlas
`peptides.widgets.mutation-cliffs` and the SARViewer mode-switch
surface).

10. Switch between **Mutation Cliffs** and **Invariant Map** modes
    on the `Sequence Variability Map` viewer (the
    `MonomerPosition` viewer's mode toggle).
11. Click the first non-empty cell at a deterministic position
    (e.g. row index 0 in the visible grid) in **Mutation Cliffs** or
    **Invariant Map**. _(Original wording was "random non-empty
    cell"; the deterministic-position substitution is the
    Playwright-translatable equivalent — see
    `unresolved_ambiguities[non-empty-cell-selection-strategy]`.)_
12. Verify that the **Context Panel** displays the necessary panels:
    **Mutation Cliffs pairs** and **Distribution**.

### Scenario 4 — Adjust Distribution panel parameters and verify rendering

13. Change parameters on the **Distribution** panel to a single
    representative non-default value (a non-default breakdown
    grouping or aggregation). _(Per SR-01 spirit, exhaustive
    parameter-matrix verification is deferred. The "arbitrarily"
    wording in the original is surfaced as
    `unresolved_ambiguities[distribution-panel-parameter-arbitrariness]`
    until a parameterized variant lands.)_

## Notes

- **Entry path.** This scenario uses the **Peptides context-panel
  Launch SAR button** entry path. The sister entry path (top-menu
  `Bio | Analyze | SAR...`) is covered by `peptide-space.md` plus
  the atlas critical_path `launch-sar-top-menu-happy-path` (p0).
- **Related bugs.** This scenario's surface intersects 7 curated
  Peptides bugs per `bug-library/peptides.yaml`. Five of them
  (GROK-19145 / GROK-14357 / GROK-18058 / github-1549 / GROK-14461)
  emit as Trigger 1 cross-scenario candidates in the chain analyzer's
  `bug_focused_candidates[]` with `sar.md` in their spans. Two
  (GROK-15934 / GROK-14298) skip per
  `bug_match_attempts_skipped[reproduction_unparseable]` because the
  scenario body does not declare a hover step (GROK-15934) or a
  filter step (GROK-14298). These remain listed in `related_bugs` as
  semantic-affecting bugs but their dedicated repros live in
  bug-focused specs outside this scenario.
- **Init-prerequisite race (GROK-17557).** The context-panel Launch
  SAR button entry path is the regression target for GROK-17557 —
  `PeptideUtils.loadComponents()` must complete before viewers are
  requested. Step 4 implicitly exercises this; for the dedicated
  init-prerequisite assertion, see the bug-focused spec
  `peptides-grok-17557-spec.ts` proposed by the chain analyzer.
- **Calculation completion timing.** Step 5's "four viewers should
  appear" plus Step 9's settings-driven reload both require explicit
  wait conditions in the Playwright translation (Web Worker computes
  for MCL/DBSCAN run async). Surfaced as
  `unresolved_ambiguities[calculation-completion-marker-needed]`.

## Original trailing metadata (preserved verbatim from source)

```json
{
  "order": 3,
  "datasets": [
    "System:DemoFiles/bio/peptides.csv"
  ]
}
```
