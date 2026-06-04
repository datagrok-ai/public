---
feature: peptides
sub_features_covered:
  - peptides.panels
  - peptides.panels.peptides
  - peptides.rendering
  - peptides.rendering.monomer
target_layer: playwright
coverage_type: smoke
produced_from: migrated
realized_as:
  - info-panels-spec.ts
original_path: public/packages/UsageAnalysis/files/TestTrack/peptides/info-panels.md
migration_date: 2026-05-28
source_text_fixes:
  - step-4-stray-double-asterisk-and-trailing-period-after-parenthesis
candidate_helpers: []
unresolved_ambiguities:
  - step-6-content-displayed-correctly-no-concrete-invariant
  - step-2-amino-acid-different-color-no-renderer-key-declared
scope_reductions:
  - id: SR-01
    check: D-SAN-02
    rationale: >-
      Original info-panels trace (info-panels-run.md, 2026-04-09) enumerated
      three Context Panel panes (Details, Peptides, Bioinformatics); migrated
      scenario verifies only Details + Peptides. The Bioinformatics accordion
      pane is registered by the Bio package, not Peptides; the peptides.panels.*
      atlas surface contains only peptides.panels.peptides and
      peptides.panels.manual-alignment (atlas L133-L154), so the Bioinformatics
      pane is out of the Peptides feature surface owned by this section and
      belongs to the Bio chain.
    verdict_status: null
  - id: SR-02
    check: D-STEP-01
    rationale: >-
      Original steps 5 (Expand each tab) and 6 (Verify panel content) merge into
      migrated steps 4 and 5; the original step-6 verification target
      (Bioinformatics sequence logo + WebLogo chart) is scoped out by SR-01, and
      the residual per-pane content invariant for Details and Peptides is
      surfaced as unresolved_ambiguities[step-6-content-displayed-
      correctly-no-concrete-invariant]. Bookkeeping merge only; the scenario
      body already constrains step 5 to the two enumerated panes.
    verdict_status: null
related_bugs:
  - GROK-17557
pyramid_layer: ui-smoke
ui_coverage_responsibility:
  - peptides-context-panel-details
  - peptides-context-panel-peptides-tab
  - macromolecule-column-monomer-coloring
ui_coverage_delegated_to: null
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-05-29-peptides-migrate-02
    timestamp: 2026-05-29T00:00:00Z
    review_round: 1
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-30-peptides-automate-02
    timestamp: 2026-05-30T20:42:00Z
    spec_runs:
      - spec: info-panels-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 63
        failure_keys: []
  d:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-05-29-peptides-migrate-02
    timestamp: 2026-05-29T00:00:00Z
    failure_keys: []
  f:
    verdict: PASS
    cycle_id: 2026-05-29-peptides-migrate-02
    timestamp: 2026-05-29T00:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-05-30-peptides-automate-02
    timestamp: 2026-05-30T23:55:00Z
    failure_keys: []
---

## Setup

1. Open the linked Peptides demo dataset (`System:DemoFiles/bio/peptides.csv`) — a TableView with the peptides grid and a `Macromolecule`-semtype `AlignedSequence` column should be the active view.

## Scenarios

### Context Panel surfaces Details and Peptides on Macromolecule column

Entry-point smoke for the Peptides info-panel surface. Verifies (a) per-cell amino-acid coloring on the `Macromolecule` column rendered by the Peptides cell-renderer, and (b) that selecting the peptides column reveals the expected accordion panes on the Context Panel — Details (platform-wide pane) and Peptides (the `peptidesPanel` info panel registered by the Peptides package for Macromolecule columns).

1. Verify each amino acid in the `AlignedSequence` column is rendered with a distinct color (the per-cell monomer renderer is installed on the Macromolecule column).
2. Click the `AlignedSequence` column title — the column becomes the current column and the Context Panel updates to reflect it.
3. On the Context Panel, confirm both the **Details** pane and the **Peptides** pane are present.
4. Expand the **Details** pane and the **Peptides** pane.
5. Confirm each expanded pane renders its expected content (Details — column metadata such as type, semtype, value statistics; Peptides — the SAR analyze UI / Launch SAR entry surface).

## Notes

- Pyramid role: this scenario is the chain's `ui-smoke` residual per `scenario-chains/peptides.yaml :: ui_coverage_plan.smoke_scenario` (no entity create/share/delete in the Peptides feature shape; shortest scenario by step count; lowest `order`). Sibling scenarios `peptides.md`, `sar.md`, `peptide-space.md` delegate the shared `peptides-context-panel-peptides-tab` entry-point flow to this scenario.
- Bug coverage: cross-cutting bug `GROK-17557` (SAR launched from the Peptides context panel raises "SeqHelper is not initialized") affects the `peptides.panels.peptides` surface this scenario exercises; the dedicated cross-cutting spec is proposed under chain `bug_focused_candidates` as `peptides-grok-17557-spec.ts`. This scenario verifies the panel-rendering smoke; the Launch-SAR-from-panel init-prerequisite invariant lives in the dedicated spec.
- Renderer-tag mapping: the original source text says "each amino acid is rendered with a different color" without specifying the renderer key. The sibling Playwright spec `info-panels-spec.ts` (test "Peptides — Info Panels") asserts `col.getTag('cell.renderer') === 'sequence'`; this mapping is not declared in the original scenario body and is surfaced as an unresolved ambiguity for downstream resolution.
- Panel-content verification: step 5 of the original asks for "content for each panel is displayed correctly". The migrated text constrains "each panel" to the two panels explicitly enumerated in step 4 (Details + Peptides) and lists representative content shapes per pane; the absence of a concrete per-pane invariant in the source is surfaced as an unresolved ambiguity.
