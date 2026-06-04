---
feature: bio
sub_features_covered:
  - bio.analyze.composition
  - bio.viewers.web-logo
target_layer: playwright
coverage_type: regression
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

Integration scenario for the `Bio | Analyze | Composition` top menu
(atlas `bio.analyze.composition` — `package.ts#L1041`) across three
canonical Macromolecule notations (FASTA, HELM, MSA). Verifies the
multi-subsystem path: top-menu dispatch → Composition WebLogo viewer
docking (`bio.viewers.web-logo`, atlas `package.ts#L449`) →
viewer-to-grid selection bridge (click-letter selects rows) → Gear
icon → Context Pane property editor surface.

The umbrella runner that exercises Composition alongside Sequence
Space and Activity Cliffs lives in `analyze.md`; this scenario is the
deeper Composition-specific scenario. WebLogo paint correctness
(letter-height ordering, color, conservation track) is atlas
`manual_only` (`bio.viewers.web-logo`) — this scenario asserts the
user-visible behaviours that drive that viewer (docking, click→select,
Gear→Context-Pane wiring) without claiming pixel-level WebLogo
correctness.

## Setup

Composition is exercised on each of three test datasets to verify
notation-agnostic behaviour across the canonical Macromolecule
notations:

- `System.AppData/Bio/tests/filter_FASTA.csv`
- `System.AppData/Bio/tests/filter_HELM.csv`
- `System.AppData/Bio/tests/filter_MSA.csv`

The implicit matrix is **3 datasets × 1 Composition flow = 3 cells**.
Each dataset has a single Macromolecule column so the multi-column
choice dialog (atlas `bio.analyze.composition` — "Multi-column choice
dialog when more than one Macromolecule column is present") does NOT
appear in this scenario; the Composition viewer docks directly on
dispatch.

## Scenarios

### Scenario 1 — Open dataset and dispatch Composition

For each dataset in `{filter_FASTA.csv, filter_HELM.csv, filter_MSA.csv}`:

1. Open the dataset from `System.AppData/Bio/tests/`. The Macromolecule
   detector (atlas `bio.detector`) classifies the sequence column
   synchronously; the table view opens.

2. On the menu ribbon, open **Bio** > **Analyze** > **Composition**
   (atlas `bio.analyze.composition`, `compositionAnalysis` —
   `package.ts#L1041`). Verify:
   - A WebLogo viewer docks (atlas `bio.viewers.web-logo`,
     `WebLogoViewer` — `package.ts#L449`); viewer presence verified
     against `grok.shell.tv.viewers` rather than a pixel-level paint
     assertion (paint correctness is atlas `manual_only`).
   - No multi-column choice dialog appears (the dataset has a single
     Macromolecule column).

### Scenario 2 — Click letter in WebLogo selects rows in the grid

For each dataset (same matrix as Scenario 1):

3. Click a letter cell in the docked Composition WebLogo viewer. Verify:
   - At least one row becomes selected in the source grid (viewer →
     grid selection bridge — the `bio.analyze.composition` →
     `bio.viewers.web-logo` click-handler routes the canvas hit to
     a grid-row selection set on the source dataframe).
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
     property surface bound (the WebLogo `WebLogoViewer` —
     `package.ts#L449` — properties: `sequenceColumnName`,
     `positionWidth`, `startPosition`, `endPosition`, `mode`, etc.).
   - The property-editor surface is present (presence of the
     property-grid is asserted; the concrete property checklist is
     deferred — see SR-01).

5. In the Context Pane, edit at least one editable property of the
   Composition WebLogo viewer. Verify:
   - At least one property accepts a user-driven edit (the
     property-grid input commits the edited value back to the
     viewer; the specific property and its expected post-edit value
     are deferred — see SR-01).
   - No error balloon is surfaced from the property edit (the
     baseline assertion: the edit path does not crash; the property
     correctness assertion is deferred per SR-01).

## Notes

- **Source-text fixes silently applied** during migration (recorded
  in frontmatter `source_text_fixes`):
  - The original heading `#### Composition Analysis manual test` was
    normalized to `# Bio | Analyze | Composition — composition
    analysis integration` to match the sister scenarios' heading
    convention (`analyze.md`, `sequence-space.md`) and to surface
    the atlas top-menu form `Bio | Analyze | Composition` in the
    scenario title.
  - The original Step 1 listed datasets as `System.AppData/Bio/tests/
    filter_FASTA.csv,` / `filter_HELM.csv,` / `filter_MSA.csv.` with
    trailing commas and a trailing period on the last item; the
    migrated Setup uses a clean bullet list.
  - The original Step 2 said "open **Bio** > **Analyze** >
    **Composition** analysis" with a trailing "analysis" word;
    atlas registers the menu item as `Bio | Analyze | Composition`
    (no "analysis" suffix) per `bio.analyze.composition`
    (`package.ts#L1041`). The migrated scenario uses the
    atlas-canonical form.
  - The original Step 1 listed datasets without explicitly framing
    them as a per-dataset loop; surfaced here under Setup as a
    3 × 1 matrix and reflected in each Scenario as a "for each
    dataset" loop.
  - The original Step 2 said "A viewer opens" as a single trailing
    phrase; promoted here to an explicit verification step
    (WebLogo-viewer presence against `grok.shell.tv.viewers` +
    multi-column choice dialog absence).
  - The original Step 3 said "Check that you select elements from
    table if you click any letter in the viewer"; rephrased here as
    the viewer-to-grid selection bridge assertion with an explicit
    "≥1 row selected" baseline.
- **Sub-features covered:**
  - `bio.analyze.composition` (atlas L449, `package.ts#L1041`) —
    top-menu `Bio | Analyze | Composition` dispatch + viewer
    docking, with multi-column dialog when more than one
    Macromolecule column is present (not exercised here).
  - `bio.viewers.web-logo` (atlas L229, `package.ts#L449`) —
    canvas-based sequence-logo viewer; pixel-level paint
    correctness is atlas `manual_only` and is NOT asserted here.
    The user-visible behaviours that drive the viewer (docking
    presence, click-letter-selects-rows bridge, Gear-icon →
    Context-Pane wiring, property-editor edit acceptance) are
    asserted at the assertable surface.
  - Maps directly onto atlas critical path
    `bio.cp.composition-analysis` (p1,
    `derived_from: public/packages/Bio/src/package.ts#L1041`).
- **Manual-only surface touched but not asserted in pixel form:**
  `bio.viewers.web-logo`. Atlas `manual_only` declares the WebLogo
  paint (letter-height ordering, color, conservation track) is
  human-inspection only; this scenario asserts viewer-docking
  presence + click-select bridge + Gear→Context-Pane wiring +
  property-edit acceptance only. Per the chain note
  (`composition-analysis.md` notes block): "this scenario
  exercises the user-visible behaviors that drive that viewer
  without claiming pixel-level WebLogo correctness — assertable
  surface (click select, Gear panel) is automatable even though
  the underlying WebLogo paint is manual-only."
- **No related bugs.** Bug-library `bio.yaml` was consulted and no
  curated bug intersects the Composition / WebLogo-viewer-via-
  Composition surface. GROK-18474 (MSA column-header click crash
  on FASTA, `bio.rendering.column-header` + `bio.viewers.web-logo`)
  is explicitly excluded by chain
  `bug_match_attempts_skipped[GROK-18474]` — this scenario does not
  exercise the MSA column-header click path. `related_bugs: []`
  in frontmatter reflects this.
- **Cross-scenario context:** `analyze.md` exercises Composition
  alongside Sequence Space and Activity Cliffs as an umbrella
  runner (its Scenario 2 also verifies Composition viewer docking;
  its Scenario 3 also exercises the Composition Gear→Context-Pane
  wiring on `filter_FASTA.csv`). `composition-analysis.md` is the
  deeper Composition-specific scenario that additionally asserts
  the click-letter-selects-rows bridge (Scenario 2 here) — the
  unique-to-this-scenario assertion.
- **Unresolved ambiguities** (carried in frontmatter
  `unresolved_ambiguities`):
  - The Context-Pane property checklist for the Composition WebLogo
    viewer's Gear surface is not declared in atlas
    (`bio.viewers.web-logo` describes the WebLogo property surface
    as `sequenceColumnName`, `positionWidth`, `startPosition`,
    `endPosition`, `mode`, etc., but does not pin a canonical
    "which subset must be present in the Context Pane after Gear
    click" checklist). The prior run log
    (`composition-analysis-run.md`) used `Show Position Labels`
    and `Skip Empty Positions` toggles as concrete picks, but
    neither atlas nor an operator decision binds those as the
    canonical edit set; see SR-01.
  - The canvas-interactivity settle window for the Composition
    WebLogo viewer (the time between viewer-docked and
    click-letter-routes-to-grid-selection) is not declared in
    atlas. The prior run log notes a ~3s settle is required on
    the 2nd/3rd dataset before the canvas hit handlers wire
    cleanly — empirical observation, not an atlas contract.
    Operator clarification is needed before Automator can encode
    a typed assertion on the settle window.
  - The click-letter canvas hit-test coordinates are not declared
    in atlas. The prior run log probes a set of x-offsets (`18,
    30, 50, 80, 120, 160`) so narrower letter columns on HELM/MSA
    still land on a letter cell; this is an empirical robustness
    pattern, not an atlas-declared selector. Operator
    clarification is welcomed if a canonical canvas hit-test
    procedure should be encoded.
- **Candidate helpers** (carried in frontmatter `candidate_helpers`):
  - `bio.flow.composition` — registered in
    `helpers-registry.yaml` as a grok-browser procedure for the
    `Bio | Analyze | Composition` flow. Surfaced as a candidate
    Automator-side helper for the top-menu dispatch portion of
    Scenario 1 Step 2.
- **Deferred via SCOPE_REDUCTION** (see frontmatter `scope_reductions`):
  - SR-01 — the correctness assertion on the edited-property state
    in Step 5 (source Step 5, "change arbitrary properties") is
    deferred until atlas or operator supplies a concrete
    Context-Pane property checklist for the Composition WebLogo
    viewer. The Gear→Context-Pane wiring + edit-acceptance flow
    itself is preserved as a step (the Gear icon opens the
    property-grid surface and at least one property accepts a
    user-driven edit).
- **`derived_from:` provenance:** atlas entries cited above derive
  from code anchors only — `bio.analyze.composition` from
  `package.ts#L1041`, `bio.viewers.web-logo` from
  `package.ts#L449`, `bio.cp.composition-analysis` from
  `public/packages/Bio/src/package.ts#L1041`. No help-doc
  derivations were used (per the binding sourcing rule).

---
{
  "order": 10
}
