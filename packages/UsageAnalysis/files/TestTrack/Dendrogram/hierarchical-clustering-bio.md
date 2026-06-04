---
feature: dendrogram
sub_features_covered:
  - dendrogram.clustering.menu.bio
  - dendrogram.clustering.dialog
  - dendrogram.api.tree-helper.calc-distance-matrix
  - dendrogram.clustering.inject-tree-for-grid
  - dendrogram.clustering.assign-clusters-dialog
  - dendrogram.api.tree-helper.cut-tree-to-grid
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Dendrogram/hierarchical-clustering-bio.md
migration_date: 2026-06-02
source_text_fixes:
  - canonicalize-headings-h1-setup-scenarios
  - promote-block-a-b-c-to-named-scenarios
  - move-step-6-centroid-median-parenthetical-to-notes
  - cite-explicit-fasta-pt-activity-csv-path-from-order-trailer
  - mark-block-c-explicitly-as-delegated-coverage-ride-along
candidate_helpers: []
unresolved_ambiguities:
  - step-6-centroid-median-mirror-chem-dialog-step-7-pattern-not-decided
scope_reductions:
  - id: SR-01
    check: A-CONT-01
    rationale: |
      Source Step 9 ("Hold Ctrl (or Cmd) and scroll the mouse wheel over
      the dendrogram. Expected: tree zooms horizontally along the X-axis
      (clamped 1×–100×); plain wheel scrolls vertically.") exercises
      tactile Ctrl+wheel horizontal zoom and plain-wheel vertical scroll
      over the dendrogram-as-grid-neighbor. Atlas
      `manual_only[dendrogram.mo.ctrl-wheel-zoom-tactile]` (sub_feature_id:
      dendrogram.clustering.inject-tree-for-grid) declares: "Ctrl+wheel
      horizontal zoom and plain-wheel vertical scroll behaviors over the
      tree-as-grid-neighbor depend on browser-level wheel-event delivery
      and platform-specific modifier handling (Ctrl on Windows/Linux vs
      Cmd on macOS). Tactile/clamping characteristics (clamp ratio,
      zoom-out floor) are out of scope for headless automation; route to
      manual." Deferred to manual visual review per atlas. The canonical
      owner of this surface in the Dendrogram chain is
      `assign-clusters.md` (its frontmatter SR-01 already records the
      same atlas dependency); per the chain's `ui_coverage_plan`, this
      scenario delegates Ctrl+wheel coverage to that scenario rather
      than duplicating it.
    verdict_status: SCOPE_REDUCTION
related_bugs: []
realized_as:
  - hierarchical-clustering-bio-spec.ts
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-migrate-01
    timestamp: 2026-06-03T11:00:00Z
    failure_keys: []
    review_round: 1
  d:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-06-03-dendrogram-migrate-01
    timestamp: 2026-06-03T10:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T15:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T15:24:00Z
    spec_runs:
      - spec: hierarchical-clustering-bio-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 74
        failure_keys: []
---

# Hierarchical Clustering (Bio) — sequence-default dialog + bio-specific build path

Bio-path integration scenario for atlas critical_path
`dendrogram.cp.hier-clustering-bio-sequence-path` (`priority: p1`):
verify the Bio top-menu surface auto-selects a MACROMOLECULE
(sequence) column as the default Features input on the Hierarchical
Clustering dialog, then build a dendrogram from the sequence column
through the Levenshtein-on-encoded-sequences distance path. Block C
is an explicit ride-along smoke on the post-build Assign Clusters
surface; per the chain's `ui_coverage_plan`, the dialog gateway
(Distance/Linkage dropdown enumeration, OK button, menu entry) is
delegated to `hierarchical-clustering-chem.md` (section ui-smoke)
and the full Assign Clusters + Ctrl+wheel zoom surface is delegated
to `assign-clusters.md` (integration). This scenario therefore
re-exercises the gateway only as a precondition to building the
bio-path tree; its specialty coverage is everything bio-path-specific
(menu.bio entry, MACROMOLECULE auto-select, Levenshtein distance,
sequence-column leaf binding).

This scenario complements:

- `hierarchical-clustering-chem.md` (section ui-smoke) — owns the
  shared hierarchical-clustering dialog gateway (Distance/Linkage
  dropdowns, Features column input, OK button, menu entry surface
  for the Chem-default path).
- `assign-clusters.md` (integration) — owns the Assign Clusters
  dialog full surface (two-way Threshold/Clusters binding,
  column creation, replace-on-rerun, Ctrl+wheel zoom deferral).
- `hierarchical-clustering-bio-api.md` (apitest, source-matrix) —
  owns the 14-combo `{euclidean, manhattan} × {7 linkages}`
  parameter matrix on the sequence path and the MACROMOLECULE
  semType precondition assertion.

## Setup

- A clean Datagrok view (no preloaded tables).
- The Bio package is installed and registered (the `Bio | Analyze`
  top-menu surface and the MACROMOLECULE sequence renderer /
  semType detector must be available).
- The Dendrogram package is installed and registered.
- Dataset: `System:AppData/Bio/samples/FASTA_PT_activity.csv` — a
  bio sample with a `sequence` column that should be detected as
  MACROMOLECULE.
- No upstream scenario state required — this scenario opens its
  own dataset and builds its own dendrogram (chain YAML
  `depends_on: []`, `produces:
  ["dendrogram-from-fasta-pt-activity-bio-dialog",
  "cluster-column-from-bio-assign"]`).

## Scenarios

### Scenario 1 — Bio dialog opens with sequence-default Features

1. Open the dataset
   `System:AppData/Bio/samples/FASTA_PT_activity.csv`. Wait for
   the `sequence` column to be detected as MACROMOLECULE (the
   semType detector triggers the sequence renderer on the column).
   - Verification: a TableView is open on the loaded DataFrame; the
     `sequence` column is visible and rendered as a sequence (the
     MACROMOLECULE semType is applied; the column renderer shows
     residue glyphs rather than raw text). No console errors.

2. From the top menu run
   **Bio | Analyze | Hierarchical Clustering...**.
   - Verification: the **Hierarchical Clustering** dialog opens
     with **Features** auto-selected to the `sequence` MACROMOLECULE
     column (the Bio menu entry calls
     `hierarchicalClusteringSequences` which seeds the Features
     input via `bySemType(MACROMOLECULE)`; see Notes). The dialog
     also exposes **Distance** and **Linkage** inputs. No console
     errors.

### Scenario 2 — Distance and Linkage dropdowns expose the canonical value sets

Preconditions: the Hierarchical Clustering dialog is open from
Scenario 1.

3. Open the **Distance** dropdown.
   - Verification: exactly two values are listed — `euclidean`,
     `manhattan`. Default selection: `euclidean`.

4. Open the **Linkage** dropdown.
   - Verification: exactly seven values are listed, in order —
     `single`, `complete`, `average`, `weighted`, `centroid`,
     `median`, `ward`. Default selection: `ward`.

Note: dropdown enumeration is owned canonically by
`hierarchical-clustering-chem.md` (section ui-smoke); this scenario
re-verifies the enumerations on the Bio menu surface because the
auto-default Features column differs (MACROMOLECULE for Bio vs
MOLECULE for Chem) and a regression in the Bio menu entry's
dialog wiring could fail dropdown population in ways the chem path
would not detect.

### Scenario 3 — Build dendrogram from the sequence column (Levenshtein path)

Preconditions: Scenarios 1-2 have run; the dialog is open.

5. With **Features** = `sequence`, **Distance** = `euclidean`,
   **Linkage** = `ward`, click **OK**.
   - Verification: a transient progress indicator
     (`Creating dendrogram ...`) is shown while the sequence
     distance matrix is computed via the encode + Levenshtein path
     (atlas sub_feature
     `dendrogram.api.tree-helper.calc-distance-matrix` declares
     "MACROMOLECULE (Levenshtein on encoded sequences)"), then a
     dendrogram neighbor is injected to the left of the grid with
     one leaf per DataFrame row. Grid ↔ tree current /
     mouseOver / selection / filter states are synchronized
     (atlas sub_feature
     `dendrogram.clustering.inject-tree-for-grid` declares the
     bidirectional sync). No console errors and no
     `Unsupported column type` error.

6. Close the dendrogram. Re-open the dialog, set
   **Distance** = `manhattan`, **Linkage** = `complete`,
   **Features** = `sequence`, click **OK**.
   - Verification: a dendrogram builds successfully (different
     shape, still one leaf per row). No `Unsupported column type`
     error and no fatal console errors.

### Scenario 4 — Shared post-build smoke ride-along (Assign Clusters column creation)

Preconditions: a bio-built dendrogram is attached to the grid
(Scenario 3 has run). This scenario is explicitly delegated
coverage — its canonical owner is `assign-clusters.md`
(integration), which exercises the Assign Clusters surface
end-to-end on a chem-built tree. The intent here is regression
protection: confirm the same Assign Clusters surface is reachable
and functional on a tree produced through the bio Levenshtein
distance path (not only the chem Tanimoto-on-fingerprints path).

7. Click the **magic wand** icon (top-left of the dendrogram
   neighbor), or right-click the tree canvas and choose
   **Assign Clusters** from the context menu.
   - Verification: the **Assign Clusters** dialog opens with a
     **Threshold** input (slider) and a **Clusters** input
     (integer, minimum 1). The two inputs are interconnected
     per atlas sub_feature
     `dendrogram.clustering.assign-clusters-dialog`. No console
     errors.

8. Set **Clusters** to `5` and click **Assign**.
   - Verification: the dialog closes and a new categorical column
     `Cluster (<threshold>)` is added to the host DataFrame (one
     row per original sequence row; the column labels each row
     with its cluster id). Atlas sub_feature
     `dendrogram.api.tree-helper.cut-tree-to-grid` declares the
     column-creation surface. No console errors.

## Notes

- **Bio menu path provenance.** The Bio top-menu entry
  `Bio | Analyze | Hierarchical Clustering...` is registered by
  the Dendrogram package via `package.g.ts#L98` (function
  `hierarchicalClusteringSequences`), not by the Bio package
  itself — per atlas sub_feature
  `dendrogram.clustering.menu.bio` and its `derived_from:
  [SRC Dendrogram:hierarchicalClusteringSequences
  public/packages/Dendrogram/src/package.g.ts#L99]`. The
  registering function seeds the dialog's Features input via
  `bySemType(MACROMOLECULE)`, which is what produces the
  sequence-default behaviour exercised in Scenario 1. The dialog
  itself is identical across all three top-menu entries
  (`Bio | Analyze`, `Chem | Analyze`, `ML | Cluster`) — they
  differ only in the default-selected feature column. The shared
  dialog gateway is owned by `hierarchical-clustering-chem.md`
  (section ui-smoke per the chain's `ui_coverage_plan`).
- **Levenshtein distance path.** For MACROMOLECULE columns,
  `TreeHelper.calcDistanceMatrix` encodes the sequence values and
  computes pairwise distances via Levenshtein
  (`tree-helper.ts#L526`-`L545`). Numeric columns use the
  per-column difference metric; MOLECULE columns use Tanimoto on
  Morgan fingerprints (Chem path). The sequence path is exercised
  here in Scenarios 3 and 4; the corresponding parameter matrix
  (14 distance × linkage combinations on the sequence path) is
  asserted by `hierarchical-clustering-bio-api.md`
  (apitest, source-matrix layer).
- **Step 6 centroid/median parenthetical (moved from inline).**
  The original source step 6 ("Re-open the dialog, set
  Distance = manhattan, Linkage = complete, Features = sequence,
  click OK.") carried an inline parenthetical "(centroid/median
  may render non-monotonic trees — expected, not a defect.)".
  Step 6 specifies the `complete` linkage, not `centroid` or
  `median`, so the parenthetical is misplaced (chain
  `unresolved_ambiguities` flagged this as priority: low). The
  parenthetical is moved here for accuracy: centroid and median
  linkages may produce a non-monotonic tree (branches with
  inversions or apparent cross-overs), which is mathematically
  expected for these linkage methods and is not a defect. The
  non-monotonic property is a structural property of the
  merge-height sequence and is asserted in the apitest matrix
  (`hierarchical-clustering-bio-api.md`, Scenario 1 covers all
  14 distance × linkage combinations on the sequence path with
  leaf-count + no-throw checks; chem-api Scenario 3 asserts the
  positional linkage-code contract). The UI smoke here does not
  exercise centroid/median linkages on the sequence path — see
  `unresolved_ambiguities:
  step-6-centroid-median-mirror-chem-dialog-step-7-pattern-not-decided`
  for the open question of whether to add a centroid/median step
  to mirror the chem-dialog Step 7 pattern.
- **Block C deferred surface (Ctrl+wheel zoom).** The original
  source Step 9 ("Hold Ctrl or Cmd and scroll the mouse wheel
  over the dendrogram. Expected: tree zooms horizontally along
  the X-axis, clamped 1×–100×; plain wheel scrolls vertically.")
  is deferred to manual visual review per atlas `manual_only[]`
  entry `dendrogram.mo.ctrl-wheel-zoom-tactile` — see frontmatter
  `scope_reductions: [SR-01]` and the rationale field there. The
  canonical owner of this surface is `assign-clusters.md`, whose
  SR-01 records the same atlas dependency. The 1×–100× clamp
  range cited by the original step is consistent with the
  sibling-scenario citation in `assign-clusters.md` Step 12
  notes (`unresolved_ambiguities:
  step-12-max-horizontal-zoom-factor-not-specified` in that
  file).
- **Section chain context.** This scenario's chain entry declares
  `pyramid_layer: integration` (multi-subsystem within the
  feature: bio-dialog + sequence build path + delegated
  post-build smoke ride-along) with `representative source
  FASTA_PT_activity.csv` per Variant C. The chain's
  `ui_coverage_responsibility` for this scenario lists
  `bio-analyze-hierarchical-clustering-menu-entry`,
  `hierarchical-clustering-dialog-macromolecule-default`, and
  `hierarchical-clustering-dialog-sequence-path` — those flows
  are owned here. Standard dialog flows (Distance/Linkage
  dropdowns, OK button, dialog open) and post-build flows
  (Assign Clusters dialog, threshold slider, clusters input,
  assign button, Ctrl+wheel zoom) are delegated to
  `hierarchical-clustering-chem.md` and `assign-clusters.md`
  respectively (chain `ui_coverage_plan.delegated_scenarios`).
- **Sub_features cited per scenario.** Scenarios 1-2 →
  `dendrogram.clustering.menu.bio`,
  `dendrogram.clustering.dialog`; Scenario 3 →
  `dendrogram.api.tree-helper.calc-distance-matrix`,
  `dendrogram.clustering.inject-tree-for-grid`; Scenario 4 →
  `dendrogram.clustering.assign-clusters-dialog`,
  `dendrogram.api.tree-helper.cut-tree-to-grid`.
- **Bug library.** `bug-library/dendrogram.yaml` lists one entry,
  GROK-13041 ("Dendrogram needs reset after filtering"), which
  affects `dendrogram.clustering.inject-tree-for-grid` /
  `dendrogram.event.selection-changed`. That bug's repro is
  owned by a separate filter-vs-sort regression scenario (atlas
  `edge_cases[dendrogram.ec.filter-does-not-trigger-remove-revert-prompt]`
  scope), not by this scenario; `related_bugs: []` in the
  frontmatter is therefore intentional.

---
{
  "order": 4,
  "datasets": ["System:AppData/Bio/samples/FASTA_PT_activity.csv"]
}
