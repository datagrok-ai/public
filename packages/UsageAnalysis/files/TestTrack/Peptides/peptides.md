---
feature: peptides
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [sar-parameters-and-weblogo-rendering]
realizes: [bio.analyze.sar]
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

# Peptides — Context panel parameter widgets and WebLogo click-to-select

Focusing the Macromolecule peptides column reveals a **Peptides** accordion in the Context Panel with Activity/Scaling/Clusters parameter widgets and a WebLogo preview. This scenario checks that changing those parameters re-renders the preview, and that clicking an amino acid on the WebLogo chart selects the matching rows in the underlying table.

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

- **Related bugs.**
  - GROK-19145 — a high Similarity threshold can lead to a null-BitSet crash, surfacing as an empty WebLogo. Step 4's parameter change is the trigger surface and step 9's selection-count assertion is the downstream observation.
  - GROK-17557 — an init-prerequisite race on the context-panel entry path (`SeqHelper is not initialized`), triggered by expanding the Peptides panel in step 3; fixed in 1.24.0, so this scenario also serves as a regression check for that fix.
  - GROK-14298 — filter-broadcast crash class on the selection backbone; step 9 touches this surface at the edge, but full coverage (including a filter step) lives in a dedicated spec, since this scenario doesn't apply a filter.
- **Delegation.** This scenario relies on the shared "click column title → Peptides panel expands" flow (steps 2-3) but delegates verifying that flow itself to `info-panels.md`; here it's only a stepping stone to reach the parameter-widget and WebLogo surfaces this scenario actually owns.
- **Deferral — non-deterministic source wording.** The original steps 4-6 said "change parameters arbitrarily" / "click an arbitrary monomer" without specifying values. This scenario preserves that intent (any valid combination / any valid monomer / any non-empty selection); the paired spec (`peptides-spec.ts`) encodes the concrete deterministic choices used in automation.

{
  "order": 2,
  "datasets": [
    "System:DemoFiles/bio/peptides.csv"
  ]
}
