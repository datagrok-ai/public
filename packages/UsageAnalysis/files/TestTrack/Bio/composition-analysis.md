---
feature: bio
target_layer: playwright
coverage_type: regression
priority: p0
realizes: [bio.cp.composition-analysis]
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/bio/composition-analysis.md
migration_date: 2026-05-31
realized_as:
  - composition-analysis-spec.ts
source_text_fixes:
  - heading-normalized-to-bio-analyze-composition-integration
  - dataset-list-formatting-trailing-commas-cleaned
  - menu-path-suffix-analysis-dropped-per-atlas-canonical-form
  - implicit-three-dataset-matrix-made-explicit
  - step-2-viewer-presence-promoted-to-explicit-assertion
  - step-3-click-letter-rephrased-as-viewer-to-grid-bridge-assertion
  - dataset-names-substituted-sample-to-filter-per-chain-canonical-bio-fixtures
  - dataset-folder-substituted-samples-to-tests-per-chain-canonical-bio-fixtures
candidate_helpers:
  - bio.flow.composition
unresolved_ambiguities:
  - context-pane-property-checklist-not-defined-in-atlas
  - weblogo-canvas-interactivity-settle-window-not-atlas-declared
  - click-letter-canvas-hit-test-coordinates-not-atlas-declared
scope_reductions:
  - id: SR-01
    check: A-CONT-01
    rationale: |
      Source Step 5 ("change arbitrary properties") is non-deterministic —
      atlas does not enumerate a canonical Context-Pane property checklist
      for the Composition WebLogo viewer (parallels `analyze.md`
      SR-02). The Gear → Context-Pane wiring is preserved as a step
      (the Gear icon opens the Composition viewer's property surface and
      at least one property accepts a user-driven edit) but the
      correctness assertion on the edited-property state is deferred
      until atlas or operator supplies a concrete checklist of
      Context-Pane properties to verify; otherwise the verification
      would silently accept no-op edits.
    verdict_status: SCOPE_REDUCTION
related_bugs: []
gate_verdicts:
  a:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T16:30:00Z
    review_round: 1
    failure_keys: []
  d:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T16:10:00Z
    failure_keys: []
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-02-bio-automate-01
    timestamp: 2026-06-02T13:25:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-02-bio-automate-01
    timestamp: 2026-06-02T16:02:00Z
    spec_runs:
      - spec: composition-analysis-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 140
        failure_keys: []
---

# Bio | Analyze | Composition — composition analysis integration

Covers the **Bio | Analyze | Composition** top menu across three
canonical Macromolecule notations (FASTA, HELM, MSA): top-menu
dispatch → Composition WebLogo viewer docking → viewer-to-grid
selection bridge (clicking a letter selects matching rows) → Gear
icon → Context Panel property editor.

The umbrella runner that exercises Composition alongside Sequence
Space and Activity Cliffs lives in `analyze.md`; this is the deeper,
Composition-specific scenario. WebLogo paint correctness (letter-height
ordering, color, conservation track) is a manual/visual check — this
scenario asserts the user-visible behaviors that drive the viewer
(docking, click-to-select, Gear-to-Context-Panel wiring) without
claiming pixel-level correctness.

## Setup

Composition is exercised on each of three test datasets to verify
notation-agnostic behaviour across the canonical Macromolecule
notations:

- `System.AppData/Bio/tests/filter_FASTA.csv`
- `System.AppData/Bio/tests/filter_HELM.csv`
- `System.AppData/Bio/tests/filter_MSA.csv`

The implicit matrix is **3 datasets × 1 Composition flow = 3 cells**.
Each dataset has a single Macromolecule column so the multi-column
choice dialog ("Multi-column choice dialog when more than one
Macromolecule column is present") does NOT appear in this scenario;
the Composition viewer docks directly on dispatch.

## Scenarios

### Scenario 1 — Open dataset and dispatch Composition

For each dataset in `{filter_FASTA.csv, filter_HELM.csv, filter_MSA.csv}`:

1. Open the dataset from `System.AppData/Bio/tests/`. The Macromolecule
   detector (atlas `bio.detector`) classifies the sequence column
   synchronously; the table view opens.

2. On the menu ribbon, open **Bio** > **Analyze** > **Composition**.
   Verify:
   - A WebLogo viewer docks (`WebLogoViewer`); viewer presence
     verified against `grok.shell.tv.viewers` rather than a
     pixel-level paint assertion (paint correctness is a
     manual/visual check).
   - No multi-column choice dialog appears (the dataset has a single
     Macromolecule column).

### Scenario 2 — Click letter in WebLogo selects rows in the grid

For each dataset (same matrix as Scenario 1):

3. Click a letter cell in the docked Composition WebLogo viewer. Verify:
   - At least one row becomes selected in the source grid (viewer →
     grid selection bridge — the click-handler routes the canvas
     hit to a grid-row selection set on the source dataframe).
   - The selection-row count is positive and consistent with the
     letter's frequency at that position (assertion checks "≥1 row
     selected"; the exact count is dataset-position-dependent and
     not atlas-canonical, so the row-count exact match is not
     asserted).

### Scenario 3 — Gear icon opens Context Pane; properties are editable

After Composition is docked in Scenario 1 on `filter_FASTA.csv`:

4. On the docked Composition WebLogo viewer, click the **Gear** icon.
   Verify:
   - The Context Pane opens with the Composition WebLogo viewer's
     property surface bound (the WebLogo `WebLogoViewer` properties:
     `sequenceColumnName`, `positionWidth`, `startPosition`,
     `endPosition`, `mode`, etc.).
   - The property-editor surface is present (presence of the
     property-grid is asserted; the concrete property checklist is
     deferred).

5. In the Context Pane, edit at least one editable property of the
   Composition WebLogo viewer. Verify:
   - At least one property accepts a user-driven edit (the
     property-grid input commits the edited value back to the
     viewer; the specific property and its expected post-edit value
     are deferred).
   - No error balloon is surfaced from the property edit (the
     baseline assertion: the edit path does not crash; the property
     correctness assertion is deferred).

## Notes

- GROK-18474 (MSA column-header click crash) doesn't apply here —
  this scenario doesn't exercise the MSA column-header click path.
- Cross-scenario context: `analyze.md` also verifies Composition
  viewer docking and the Gear → Context Panel wiring as part of its
  umbrella run across Sequence Space, Activity Cliffs, and
  Composition. This scenario additionally asserts the
  click-letter-selects-rows bridge, which `analyze.md` does not.
- Deferrals: the exact list of properties that must appear on the
  Composition WebLogo viewer's Context Panel, and the correctness of
  an edited property's value, aren't pinned down to a concrete
  checklist yet — the Gear → Context Panel flow and the "at least
  one property is editable" check are still asserted. Likewise, the
  settle time needed after docking before click-to-select reliably
  works, and the exact pixel offsets for hitting a letter in the
  canvas, are empirical values rather than a documented contract.

---
{
  "order": 10
}
