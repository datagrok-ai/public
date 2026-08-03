---
feature: peptides
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [peptide-space-sar-with-mcl]
realizes: [bio.analyze.sar, bio.analyze.sequence-space, peptides.active-peptide-selection, peptides.logo-summary-table]
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/peptides/peptide-space.md
migration_date: 2026-05-27
source_text_fixes:
  - step-1-bio-analyze-sar-arrow-notation-normalized-to-pipe-form
  - step-1-implicit-dataset-prerequisite-promoted-to-explicit-setup
  - step-3-wrench-button-clarified-as-sar-toolbar-settings-icon
  - step-4-arbitrary-parameters-narrowed-to-mcl-clustering-representative-change
  - step-5-mcl-different-results-sharpened-to-rerender-and-cluster-count-delta
candidate_helpers:
  - peptides.launchSarFromTopMenu
  - peptides.openSarSettingsViaWrenchButton
  - peptides.adjustMclClusteringParameter
  - peptides.verifyMclViewerRerenderedAfterSettingsChange
  - peptides.captureMclClusterCount
unresolved_ambiguities:
  - step-1-implicit-dataset-prerequisite-not-declared-in-source
  - step-2-calculation-completion-marker-needed
  - step-5-mcl-viewer-different-results-no-objective-criterion
scope_reductions:
  - id: SR-01
    check: A-CONT-01
    rationale: |
      Step 4 "Adjust arbitrary parameters and click OK" cannot translate
      to a deterministic Playwright assertion without per-parameter
      expected-state declaration across the full SAR Settings surface
      (General, Viewers, Columns, Sequence Space, MCL). Defer exhaustive
      per-parameter verification to a future parameterized
      peptides-settings-dialog spec (sister of sar.md SR-01); in this
      scenario, exercise the OK roundtrip with a single representative
      MCL-parameter change to establish the settings -> MCL-Viewer
      re-render contract. The arbitrariness wording is captured in
      unresolved_ambiguities[step-5-mcl-viewer-different-results-no-objective-criterion]
      for the future variant.
    verdict_status: SCOPE_REDUCTION
  - id: SR-02
    check: A-CONT-01
    rationale: |
      Step 5 "MCL Viewer should output different results" lacks an
      objective verification criterion (different from what reference
      state? — initial MCL output before settings adjust is not
      captured in the source). Defer cluster-by-cluster output diff
      to a future spec that captures pre/post artifacts; in this
      scenario, sharpen to a deterministic property — viewer
      re-render observed AND cluster count delta non-zero — that is
      Playwright-translatable without exhaustive output comparison.
      The looseness in the source wording is captured in
      unresolved_ambiguities[step-5-mcl-viewer-different-results-no-objective-criterion].
    verdict_status: SCOPE_REDUCTION
pyramid_layer: integration
ui_coverage_responsibility:
  - top-menu-bio-analyze-sar
  - sar-settings-wrench-button
  - mcl-viewer-parameter-driven-rerender
ui_coverage_delegated_to: info-panels.md
realized_as:
  - peptide-space-spec.ts
related_bugs:
  - GROK-19145
  - GROK-14357
  - GROK-18058
  - github-1549
  - GROK-14461
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: 2026-05-29-peptides-migrate-02
    timestamp: 2026-05-29T00:00:00Z
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-05-29-peptides-migrate-02
    timestamp: 2026-05-29T00:00:00Z
    review_round: 1
    failure_keys: []
  f:
    verdict: PASS
    cycle_id: 2026-05-29-peptides-migrate-02
    timestamp: 2026-05-29T00:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-05-30-peptides-automate-02
    timestamp: 2026-05-31T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-30-peptides-automate-02
    timestamp: 2026-05-31T00:09:30Z
    spec_runs:
      - spec: peptide-space-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 47
        failure_keys: []
---

# Peptide Space — SAR launch with sequence-space dim-reduction and MCL clustering (top-menu entry path)

Launching SAR from the top menu with Sequence Space and MCL clustering enabled should produce a sequence-space scatter plot, an MCL cluster viewer, and a cluster-aware Logo Summary Table. This scenario also checks that adjusting an MCL parameter through the Settings dialog re-renders the MCL viewer with a different cluster outcome.

## Setup

- **Dataset:** `System:DemoFiles/bio/peptides.csv` — the same demo
  table used by `info-panels.md`, `peptides.md`, and `sar.md`. The
  original trailing JSON metadata omits a `datasets:` field; the
  scenario implicitly depends on a pre-loaded Macromolecule table
  because Step 1 invokes the top-menu SAR launcher and atlas
  `peptides.workflow.sar-dialog` (`peptidesDialog`) gates on an
  active Macromolecule column and at least one numerical activity
  column. The implicit prerequisite is promoted to this explicit
  Setup step per the chain analyzer's low-priority
  `unresolved_ambiguity` finding. See
  `unresolved_ambiguities[step-1-implicit-dataset-prerequisite-not-declared-in-source]`.
- **Initial table state:** an open Datagrok view holding the
  peptides table with a Macromolecule sequence column
  (`AlignedSequence`) and a numerical activity column.

## Scenarios

### Scenario 1 — Launch SAR from the top menu and verify the SAR analysis layout with sequence-space + MCL

This exercises the **top-menu** SAR entry path
(`Bio | Analyze | SAR...`) — the sister of `sar.md`'s
context-panel-button entry path. Maps to atlas critical_path
`peptide-space-sar-with-mcl` (p1) and atlas interaction
`sar-with-sequence-space-and-mcl`. With the SAR Settings dialog's
**Sequence Space** and **MCL** params enabled (defaults of the
`peptidesDialog` config), the resulting analysis dock includes the
sequence-space scatter plot, the **MCL Viewer** (cluster max
activity), and the cluster-aware **Logo Summary Table**.

1. Invoke `Bio | Analyze | SAR...` from the top menu. The
   `Analyze Peptides` config dialog opens (atlas
   `peptides.workflow.analyze-ui`).
2. Accept the default config and click **OK**. Wait for calculation
   results to complete — MCL/DBSCAN clustering runs in a Web Worker
   and the `Most Potent Residues` / MCL Viewer panels signal
   completion when the cluster columns populate. _(Original wording
   "Wait for calculation results" lacks an explicit completion
   marker; see `unresolved_ambiguities[step-2-calculation-completion-marker-needed]`
   for the Playwright translation requirement.)_

### Scenario 2 — Adjust an MCL setting via the wrench button and verify the MCL Viewer re-renders

This exercises the **SAR Settings dialog** flow (atlas
`peptides.widgets.settings-dialog`) via the **wrench-icon button**
on the SAR analysis toolbar, and its propagation to the MCL Viewer
(atlas `peptides.viewers.cluster-max-activity` backed by
`peptides.compute.calculate-cluster-statistics`).

Per `scope_reductions[SR-01]`, exhaustive per-parameter
verification across the full SAR Settings surface (General,
Viewers, Columns, Sequence Space, MCL) is deferred to a future
parameterized `peptides-settings-dialog` spec; this scenario
asserts only the **MCL settings -> MCL Viewer re-render** contract
via one representative MCL-parameter change.

3. Click the **wrench-icon Settings** button on the SAR analysis
   toolbar to open the settings dialog (atlas
   `getSettingsDialog(model)` accordion — General / Viewers /
   Columns / Sequence Space / MCL).
4. In the **MCL** section, change one representative parameter to a
   non-default value (e.g. the MCL `inflation` parameter, or
   `maxIterations`). Click **OK**. _(Per SR-01, exhaustive
   per-parameter verification is deferred; the original
   "arbitrary parameters" wording is preserved as the unresolved
   ambiguity that motivates the future parameterized variant —
   see `unresolved_ambiguities[step-5-mcl-viewer-different-results-no-objective-criterion]`.)_
5. Verify that the **MCL Viewer** re-renders after the OK round
   trip — concretely:
   - The MCL Viewer's cluster count differs from its pre-change
     value (cluster count delta non-zero), OR
   - The MCL Viewer's cluster-membership composition changes
     (at least one cluster's member set differs).

   _(Per SR-02, this sharpens the source's "MCL Viewer should
   output different results" to a deterministic two-arm property
   that is Playwright-translatable without capturing pre/post
   cluster-by-cluster output diffs. The full output-diff
   verification remains an unresolved ambiguity for a future
   spec.)_

## Notes

- **Entry path.** This scenario uses the top-menu (`Bio | Analyze | SAR...`) entry path. The sister entry path — the Peptides context-panel **Launch SAR** button — is covered by `sar.md` (also the regression target for GROK-17557, an init-prerequisite race on that path).
- **MCL Viewer naming.** The cluster scatter viewer is registered internally as "Active peptide selection" but is referred to informally as the **MCL Viewer** throughout this scenario and in the Datagrok UI — both names refer to the same viewer.
- **Related bugs.** This scenario's surface touches five curated Peptides bugs, each with its own dedicated repro test elsewhere:
  - GROK-19145 — a post-OK compute crash on an edge-case Similarity threshold value. This scenario exercises the settings-dialog surface generically (one representative MCL parameter change), not the full threshold matrix.
  - GROK-14357 — a settings change silently failing to propagate to the Logo Summary Table. This scenario asserts the MCL-viewer re-renders; the LST-specific repro lives elsewhere.
  - GROK-18058 — a validation gap for an empty Activity input in the SAR dialog. This scenario always uses a valid default config.
  - github-1549 — the Monomer-Position viewer renders empty on different-length (ragged) sequences. This scenario uses the demo dataset, which has uniform-length sequences.
  - GROK-14461 — a saved Peptides project's layout isn't restored on reopen. This scenario doesn't exercise save/reopen (see `sar-save-reopen.md`).
- **Deferral — calculation completion timing.** "Wait for calculation results" (Step 2) and the settings-driven MCL re-render (Step 5) both depend on a Web Worker completing (MCL/DBSCAN can take some time); the exact wait condition is left to the automation, not pinned down as a fixed duration here.
- **Deferral — representative parameter only.** The Settings dialog exposes a full accordion (General, Viewers, Columns, Sequence Space, MCL). Step 4 narrows the original "change arbitrary parameters" wording to a single representative MCL parameter change, so Step 5's re-render assertion has a deterministic cause. Exhaustive per-parameter verification across the whole dialog is not covered here.

## Original trailing metadata (preserved verbatim from source)

```json
{
  "order": 4
}
```
