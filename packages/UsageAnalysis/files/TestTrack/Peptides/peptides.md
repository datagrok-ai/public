---
feature: peptides
sub_features_covered:
  - peptides.panels.peptides
  - peptides.rendering.weblogo-header
  - peptides.widgets.distribution
  - peptides.util.modify-selection
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/peptides/peptides.md
migration_date: 2026-05-28
source_text_fixes:
  - open-linked-datasets-resolved-to-system-demofiles-bio-peptides-csv
  - peptides-column-title-clarified-as-macromolecule-column-header
  - peptides-panel-context-panel-affordance-clarified-as-context-panel-accordion
candidate_helpers:
  - openLinkedDataset
  - expandContextPanelAccordion
  - readPeptidesPanelParameter
  - clickWebLogoMonomerAtPosition
  - assertSelectedRowCount
unresolved_ambiguities:
  - step-4-arbitrary-parameter-change-needs-deterministic-values-for-activity-scaling-clusters
  - step-5-arbitrary-amino-acid-selection-needs-deterministic-target-monomer-position
  - step-6-some-rows-selected-needs-quantitative-assertion-threshold
  - step-5-weblogo-chart-locator-under-peptides-pane-needs-canonical-selector-path
scope_reductions: []
pyramid_layer: integration
ui_coverage_responsibility:
  - peptides-context-panel-peptides-tab
  - peptides-panel-parameters
  - weblogo-amino-acid-click-selection
ui_coverage_delegated_to: info-panels.md
realized_as:
  - peptides-spec.ts
related_bugs:
  - GROK-19145
  - GROK-17557
  - GROK-14298
gate_verdicts:
  b:
    verdict: PASS
    cycle_id: 2026-05-30-peptides-automate-02
    timestamp: 2026-05-30T23:50:00Z
    spec_runs:
      - spec: peptides-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 68
        failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-05-30-peptides-automate-02
    timestamp: 2026-05-30T00:00:00Z
    failure_keys: []
  d:
    verdict: EVIDENCE_GAP
    cycle_id: 2026-05-29-peptides-migrate-02
    timestamp: 2026-05-29T14:10:00Z
    failure_keys: []
    evidence_needed: >
      Independent Gate D re-derivation, node critic_d_dispatch, cycle

      2026-05-29-peptides-migrate-02. Fresh isolated invocation: the
      EVIDENCE_GAP

      block already present in this slot is the artifact I overwrite, not an

      input. Every filesystem, dispatch, and decision-log fact below was

      re-verified directly by this invocation.

      All mechanical / fail-fast checks PASS on the on-disk artifact, so I did

      NOT fail-fast and proceeded to the content checks:

      - D-STRUCT-MECH-03: all 8 required migration fields present (feature,
        sub_features_covered=4, target_layer=playwright,
        coverage_type=regression, produced_from=migrated, original_path,
        migration_date=2026-05-28, related_bugs=3); deprecated migrated_from
        absent.
      - D-STRUCT-MECH-05: original_path target resolves to a real on-disk file,
        so the mechanical exists predicate holds (the self-reference problem
        below is a content-check concern, not a mechanical one).
      - D-FRONTMATTER-PHASE1-01: all four Phase 1 fields present as parseable
        YAML lists (scope_reductions=[] acceptable).
      - D-FRONTMATTER-PHASE1-02: source_text_fixes (3) and
        unresolved_ambiguities (4) are kebab-case slugs <= 80 chars, no
        duplicates; candidate_helpers (5) are bare dotted function names with no
        (args) suffix and no duplicates; scope_reductions=[] vacuously valid.
      The content checks (D-STEP-01, D-STEP-02, D-EDGE-01, D-STRUCT-01,

      D-STRUCT-02, D-SAN-02, D-MERIT-01, D-MERIT-02) each require an independent

      diff against a DISTINCT raw TestTrack original, which is not reachable

      this invocation:

      - A glob over TestTrack/[Pp]eptides/ returns exactly ONE peptides.md on
        disk (the capitalized .../TestTrack/Peptides/peptides.md). Windows is
        case-insensitive, so the lowercase original_path
        (.../TestTrack/peptides/peptides.md) and the on-disk capitalized path
        collapse to a single inode. I read the file at the lowercase original_path
        and it is byte-for-byte the migrated artifact (same migrated frontmatter
        including this gate_verdicts.d slot, same body). original_path therefore
        self-references the migrated file rather than naming a distinct
        pre-migration original.
      - The migrator dispatch (peptides.md.migrator.dispatch.yaml) records
        verdict PASS with no_op=true and "no destructive body Write performed --
        already fully migrated," so no separate pre-migration source text exists
        as a readable artifact to diff against.
      - The feature=peptides decision-log slice independently corroborates the
        unreachability: peptides.md has now logged repeated recon-d EVIDENCE_GAP
        iterations across cycles 2026-05-28-peptides-migrate-02 and
        2026-05-29-peptides-migrate-02, each terminated by recon-chain-analyzer
        LOOP_CAP_EXCEEDED (cap 3, last_verdict EVIDENCE_GAP). No failed_attempts
        entry preserves a peptides original excerpt that could substitute for the
        diff.
      Under Gate D independence, the migrated body's self-described provenance

      (the Source-provenance Note, the preserved trailing JSON datasets block,

      the determinism gaps surfaced as unresolved_ambiguities) is the producer's

      own account and cannot stand in for an independent diff. No silent content

      drop is confirmed, so this is EVIDENCE_GAP, not FAIL.

      Recon to resolve (orchestrator-side, per EVIDENCE_GAP semantics):

      1. Retrieve the pre-migration baseline blob of
         TestTrack/peptides/peptides.md, or an operator-supplied raw original,
         and repoint original_path at that distinct raw original (handoff
         invariant: original_path must reference the raw original, never the
         migrated artifact).
      2. With a distinct original in hand, re-run the content checklist: confirm
         every original numbered step and Expected-result assertion is mapped in
         the migrated body or acknowledged in the four Phase 1 fields, and that
         no source text was silently dropped.
      Because the prior cycles already hit the Gate D recon loop cap on this

      identical single-inode root cause, re-invoking Gate D against a

      self-referential original will not converge. The orchestrator should
      either

      supply a distinct raw original or route this scenario through the

      already-migrated / atlas-driven Critic D bypass (per Batch 6.54

      per_scenario_migrate). This is a tool-reachability gap on the original,
      NOT

      a confirmed content drop.
  a:
    verdict: PASS
    cycle_id: 2026-05-28-peptides-migrate-01
    timestamp: 2026-05-28T09:30:00Z
    review_round: 1
    failure_keys: []
  f:
    verdict: PASS
    cycle_id: 2026-05-29-peptides-migrate-02
    timestamp: 2026-05-29T00:00:00Z
    failure_keys: []
---

## Setup

- **Dataset:** `System:DemoFiles/bio/peptides.csv` (linked dataset).
- **Feature surface exercised:** Peptides context panel (peptides tab), parameter widgets (Activity / Scaling / Clusters), WebLogo column-header renderer on the Macromolecule column, collaborative-selection backbone propagating WebLogo clicks into DataFrame row selection.
- **Atlas anchors:**
  - `peptides.panels.peptides` (Context Panel Peptides tab — atlas L140 `peptidesPanel`).
  - `peptides.rendering.weblogo-header` (atlas L280 `setWebLogoRenderer`; click-to-select wiring).
  - `peptides.widgets.distribution` (atlas L304 `getDistributionWidget`; parameter-driven histogram).
  - `peptides.util.modify-selection` (atlas L519 `modifySelection`; shift/ctrl-aware selection update).
- **Chain context (`scenario-chains/peptides.yaml`):** scenario classifies as `pyramid_layer: integration`; delegates the shared `peptides-context-panel-peptides-tab` flow to `info-panels.md` smoke (per `ui_coverage_delegated_to: info-panels.md`); owns specialty parameter-adjustment + WebLogo selection coverage directly.
- **Existing realisation:** `public/packages/UsageAnalysis/files/TestTrack/Peptides/peptides-spec.ts` (`Peptides — SAR parameters and WebLogo`).

## Scenarios

### Scenario: Peptides panel parameter widgets are wired to the dataframe

**Given** the linked dataset `System:DemoFiles/bio/peptides.csv` is open as the active TableView,
**When** the user clicks the title of the Macromolecule peptides column,
**Then** the Context Panel surfaces the `Peptides` accordion entry (atlas-mapped to `peptides.panels.peptides`).

**Steps**

1. Open the linked dataset `System:DemoFiles/bio/peptides.csv` (per the trailing JSON `datasets` block of the original scenario). Wait for the TableView to render.
2. Click the title cell of the Macromolecule peptides column to focus it; verify the Context Panel updates to reflect a Macromolecule column selection.
3. On the Context Panel, locate the `Peptides` accordion entry and expand it. Verify it renders without error (mapped to `peptides.panels.peptides` — the panel refuses if the dataframe has no numerical column; pre-condition: the linked CSV has a numerical activity column, so the panel must render).

**Expected result**

- The Context Panel exposes a `Peptides` accordion entry when a Macromolecule column is focused.
- Expanding it renders the panel's parameter widgets (Activity / Scaling / Clusters) and a WebLogo-driven preview.

### Scenario: Activity / Scaling / Clusters parameter changes propagate to the panel

**Given** the `Peptides` accordion is expanded on the Context Panel,
**When** the user changes the `Activity`, `Scaling`, and `Clusters` parameters,
**Then** the panel re-renders the preview (distribution / WebLogo) using the updated parameter values.

**Steps**

4. With the `Peptides` panel expanded, change the values of the `Activity`, `Scaling`, and `Clusters` parameter widgets. (Specific values needed for deterministic playwright translation — see `unresolved_ambiguities[step-4-arbitrary-parameter-change-needs-deterministic-values-for-activity-scaling-clusters]`. The original scenario specifies "arbitrarily" — for the existing `peptides-spec.ts` realisation, the deterministic value set is the one already encoded in that spec.)
5. Observe that the panel re-renders after each parameter change (the `getDistributionWidget` histogram and the WebLogo preview update). No console errors fire (atlas-mapped to `peptides.widgets.distribution` re-render path).

**Expected result**

- Parameter widget changes trigger panel re-render without console errors.
- The preview reflects the updated parameter values (distribution histogram redraws; WebLogo preview re-rendered).

### Scenario: Clicking a WebLogo monomer selects matching rows in the DataFrame

**Given** the `Peptides` panel is expanded with a WebLogo preview rendered,
**When** the user clicks an amino acid (monomer) on the WebLogo chart under the peptides pane,
**Then** rows of the underlying DataFrame whose sequence carries that monomer at the clicked position are added to the selection (atlas-mapped to `peptides.rendering.weblogo-header` click handler routing into `peptides.util.modify-selection` → DataFrame `BitSet`).

**Steps**

6. Locate the WebLogo chart inside the expanded `Peptides` panel (locator path needs codifying for deterministic playwright translation — see `unresolved_ambiguities[step-5-weblogo-chart-locator-under-peptides-pane-needs-canonical-selector-path]`).
7. Click an amino acid (monomer) on the WebLogo chart. (Original specifies "arbitrary" monomer — for deterministic playwright translation a specific (monomer, position) target is required; see `unresolved_ambiguities[step-5-arbitrary-amino-acid-selection-needs-deterministic-target-monomer-position]`.)
8. Verify that the underlying DataFrame selection updates — at least one row is now selected (atlas-mapped to `peptides.util.modify-selection` populating the unified Selection map and `peptides.util.get-selection-bitset` projecting it to `DG.BitSet`).
9. Assert the selected-row count is greater than zero. (Original wording is "some rows were selected" — exact count not bounded by the source; see `unresolved_ambiguities[step-6-some-rows-selected-needs-quantitative-assertion-threshold]`. The existing realisation `peptides-spec.ts` should encode the deterministic expected count for the chosen (monomer, position) target.)

**Expected result**

- Clicking a monomer on the column-header WebLogo populates the DataFrame selection (`df.selection` cardinality > 0).
- Selection propagation is observable to any wired sibling viewer / property-panel widget (collaborative-selection backbone — atlas L370 `peptides.model.fire-bitset-changed` is the broadcast point, but for this scenario only the DataFrame `BitSet` mutation is asserted; full multi-viewer propagation is exercised by the higher-level SAR integration scenario `sar.md`).

## Notes

- **Layer choice (playwright).** Step 7 requires DOM interaction with the custom WebLogo column-header renderer (`setWebLogoRenderer` installs a per-column custom renderer with mouse handlers — atlas L280). The collaborative-selection assertion (step 9) is UI-driven; playwright is the right layer. Realisation: existing `peptides-spec.ts` sibling.
- **Coverage-type choice (regression).** Per chain YAML `pyramid_layer: integration` + A-LAYER-ALIGN-01 advisory mapping (`integration -> usually regression may be smoke for API-smoke`). This scenario is UI-driven multi-subsystem coverage (Context Panel + parameter widgets + WebLogo renderer + DataFrame `BitSet` sync) — `regression` is the appropriate test-kind.
- **Bug context.**
  - `GROK-19145` — high similarity threshold leads to a `setTrue` on null BitSet, surfaced via empty WebLogo. Step 4 (parameter changes) is the trigger surface and step 9 (selection-count assertion) is the downstream observation. If the parameter values chosen for step 4 land in the bug's regression-risk envelope, this scenario co-exercises the regression target.
  - `GROK-17557` — Peptides context-panel entry path init-prerequisite race (`SeqHelper is not initialized`). Step 3 (expand Peptides panel) is the trigger; bug status is `fixed in 1.24.0`, so this scenario co-exercises a regression target for the fix.
  - `GROK-14298` — filter broadcast crash + perf combo via `fireBitsetChanged` collaborative-selection backbone (`peptides.model.fire-bitset-changed`). Step 9 verification edges into this surface; full coverage is the dedicated cross-cutting spec (`bug_focused_candidates[GROK-14298]` is skipped per chain `bug_match_attempts_skipped` — no filter step in this scenario body).
- **Delegation.** Per chain `ui_coverage_plan.delegated_scenarios`, `peptides.md` delegates the shared `peptides-context-panel-peptides-tab` flow (entry-point context-panel-rendering smoke) to `info-panels.md`. The Macromolecule-column-title click (step 2) and Peptides-panel-expand (step 3) re-traverse the delegated flow only to reach the parameter-widget and WebLogo surfaces this scenario owns directly.
- **Helpers — registry status.** No existing helper in `helpers-registry.yaml` covers WebLogo monomer click, peptides-panel parameter read/write, or the linked-dataset open shortcut used here. Candidates surfaced in frontmatter `candidate_helpers[]` for Automator-side helper-promotion review; no helper is invented in this scenario body — concrete locators land in `peptides-spec.ts`.
- **Source provenance for "arbitrarily" / "arbitrary".** The original scenario uses non-deterministic phrasing in steps 4, 5, and 6. Migrated form preserves the original intent (any valid parameter combination; any valid monomer; any non-empty selection) and explicitly surfaces the determinism gap as unresolved ambiguities — the existing `peptides-spec.ts` realisation encodes the deterministic choices; reverse-validation against that spec is the right place to close the gap.
- **Chain dependency.** The chain YAML records `depends_on: []` for this scenario — it opens the dataset directly via step 1 (no cross-scenario fixture handoff). No `Setup` precondition imported from another migrated scenario.

{
  "order": 2,
  "datasets": [
    "System:DemoFiles/bio/peptides.csv"
  ]
}
