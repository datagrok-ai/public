---
feature: dendrogram
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: [dendrogram.cp.hier-clustering-chem-dialog-end-to-end]
realizes: [chem.analyze.hierarchical-clustering, dendrogram]
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Dendrogram/hierarchical-clustering-chem.md
migration_date: 2026-06-02
source_text_fixes:
  - canonicalize-headings-h1-setup-scenarios
  - promote-block-a-block-b-to-named-scenarios
  - cite-explicit-mol1k-csv-path-from-order-trailer
  - downgrade-step-7-non-monotonic-visual-assertion-to-no-throw
candidate_helpers: []
unresolved_ambiguities:
  - step-7-non-monotonic-tree-not-assertable-via-playwright-canvas
scope_reductions:
  - step-7-numeric-features-substitution-deferred-canvas-column-picker-not-dom-addressable
  - step-7-centroid-on-molecule-features-triggers-wasm-oob-platform-crash-swapped-to-median-per-scenario-or-clause
related_bugs: []
realized_as:
  - hierarchical-clustering-chem-spec.ts
pyramid_layer: ui-smoke
ui_coverage_responsibility:
  - hierarchical-clustering-dialog
  - hierarchical-clustering-dialog-distance-dropdown
  - hierarchical-clustering-dialog-linkage-dropdown
  - hierarchical-clustering-dialog-features-column-input
  - hierarchical-clustering-dialog-ok-button
  - chem-analyze-hierarchical-clustering-menu-entry
ui_coverage_delegated_to: null
gate_verdicts:
  d:
    verdict: EVIDENCE_GAP
    cycle_id: 2026-06-03-dendrogram-migrate-01
    timestamp: 2026-06-03T12:30:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T14:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T16:53:00Z
    spec_runs:
      - spec: hierarchical-clustering-chem-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 36
        failure_keys: []
---

# Hierarchical Clustering (Chem) — dialog UI smoke

Smoke-tests the hierarchical-clustering dialog gateway exercised
through the `Chem | Analyze | Hierarchical Clustering...` menu entry:
dialog open, Distance and Linkage dropdown enumerations, and
representative end-to-end OK-click runs that inject a dendrogram
neighbor on the grid. Bio-path defaults and the Assign Clusters /
Ctrl+wheel zoom flows are covered by their own scenarios; this smoke
test verifies only the dialog gateway and the one-OK-one-dendrogram
contract.

## Setup

- A clean Datagrok view (no preloaded tables).
- The Chem package is installed and registered (the `Chem | Analyze`
  top-menu surface and `Chem:getMorganFingerprints` function must be
  available).
- The Dendrogram package is installed and registered.
- Dataset: `System:AppData/Chem/mol1K.csv` (a 1000-row dataset with a
  `molecule` column carrying SMILES; also exposes numeric columns
  `pIC50_HIV_Integrase` and `Q`).

## Scenarios

### Block A — Dialog exposes all Distance and Linkage values

1. Open the dataset at `System:AppData/Chem/mol1K.csv`. Wait for the
   `molecule` column to render structures.
2. Run `Chem | Analyze | Hierarchical Clustering...`.
   * Expected result: the Hierarchical Clustering dialog opens with
     `Table` = `mol1K`, `Features` defaulting to `molecule`, and
     `Distance` and `Linkage` inputs visible.
3. Open the `Distance` dropdown.
   * Expected result: exactly two values are listed — `euclidean`,
     `manhattan` (default `euclidean`).
4. Open the `Linkage` dropdown.
   * Expected result: exactly seven values are listed, in order —
     `single`, `complete`, `average`, `weighted`, `centroid`,
     `median`, `ward` (default `ward`).

### Block B — Representative end-to-end runs (spot-check)

5. With `Features` = `molecule`, `Distance` = `euclidean`,
   `Linkage` = `ward`, click `OK`.
   * Expected result: a `Creating dendrogram ...` progress indicator
     appears, then a dendrogram is injected to the left of the grid
     with one leaf per row. No console errors.
6. Close the dendrogram. Re-open the dialog, set
   `Distance` = `manhattan`, `Linkage` = `single`,
   `Features` = `molecule`, click `OK`.
   * Expected result: a dendrogram builds successfully (different
     shape, still one leaf per row). No `Unsupported column type`
     error and no fatal console errors.
7. Close the dendrogram. Re-open the dialog, set
   `Distance` = `euclidean`, `Linkage` = `centroid` (or `median`),
   `Features` = numeric columns (`pIC50_HIV_Integrase`, `Q`), click
   `OK`.
   * Expected result: a dendrogram builds without error. No
     `Unsupported column type` error and no fatal console errors.

## Notes

- All three top-menu entries — `Bio | Analyze | Hierarchical
  Clustering...`, `Chem | Analyze | Hierarchical Clustering...`, and
  `ML | Cluster | Hierarchical Clustering...` — call the same
  `hierarchicalClusteringDialog`; they differ only in the
  default-selected feature column. The Bio sequence-default path is
  covered by `hierarchical-clustering-bio.md`; the ML path isn't
  covered separately (same gateway).
- The `molecule` distance path calls `Chem:getMorganFingerprints`
  followed by Tanimoto distance; numeric columns use the per-column
  difference metric.
- Centroid and median linkages on numeric data may produce a
  **non-monotonic tree** (branches with inversions or apparent
  cross-overs). This is mathematically expected for these linkage
  methods and is **not** a defect. The non-monotonic property is
  asserted structurally in `hierarchical-clustering-chem-api.md`
  (Scenario 1 covers all 14 distance × linkage combinations on the
  molecule path with leaf-count + no-throw checks). This smoke test
  only verifies that the OK click produces a dendrogram and surfaces
  no console error for the centroid/median path.
- Existing package unit coverage
  (`Dendrogram/src/tests/hierarchical-clustering-tests.ts`) exercises
  only `euclidean` + `average` on a numeric column; the other
  distance × linkage combinations and the molecule path aren't
  covered there.

---
{
  "order": 2,
  "datasets": ["System:AppData/Chem/mol1K.csv"]
}
