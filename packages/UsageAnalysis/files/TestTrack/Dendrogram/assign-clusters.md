---
feature: dendrogram
sub_features_covered:
  - dendrogram.clustering.menu.chem
  - dendrogram.clustering.dialog
  - dendrogram.clustering.inject-tree-for-grid
  - dendrogram.clustering.assign-clusters-dialog
  - dendrogram.api.tree-helper.cut-tree-to-grid
  - dendrogram.event.context-menu
  - dendrogram.viewer
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Dendrogram/assign-clusters.md
migration_date: 2026-06-02
source_text_fixes:
  - community-post-annotations-moved-to-notes
  - block-d-replace-warning-step-broken-out-as-its-own-scenario
candidate_helpers: []
unresolved_ambiguities:
  - step-7-clusters-to-threshold-tolerance-not-pinned
  - step-12-max-horizontal-zoom-factor-not-specified
scope_reductions:
  - id: SR-01
    check: A-CONT-01
    rationale: |
      Source Steps 11-13 (plain mouse-wheel vertical scroll, Ctrl+wheel
      horizontal zoom with clamping, double-click empty area resets
      horizontal zoom) exercise tactile wheel/zoom over the dendrogram
      canvas. Atlas `manual_only[dendrogram.mo.ctrl-wheel-zoom-tactile]`
      (sub_feature_id: dendrogram.clustering.inject-tree-for-grid)
      declares: "Ctrl+wheel horizontal zoom and plain-wheel vertical
      scroll behaviors over the tree-as-grid-neighbor depend on
      browser-level wheel-event delivery and platform-specific modifier
      handling (Ctrl on Windows/Linux vs Cmd on macOS). Tactile/clamping
      characteristics (clamp ratio, zoom-out floor) are out of scope for
      headless automation; route to manual." Deferred to manual visual
      review per atlas; the migrated scenario keeps Assign Clusters
      column-creation flows on the automated path.
    verdict_status: SCOPE_REDUCTION
  - id: SR-02
    check: A-CONT-01
    rationale: |
      Source Step 10 ("the dendrogram shows exactly where the cut occurs
      — the split position is rendered on the tree and updates live as
      the threshold / cluster count changes") asserts a canvas-pixel
      output. Atlas `manual_only[dendrogram.mo.tree-canvas-visual-regression]`
      (sub_feature_id: dendrogram.viewer) declares: "Tree layout on the
      canvas...is canvas-pixel output. Asserting it deterministically
      requires pixel-diff infrastructure that is out of scope for the
      deterministic UI automation layer; route to manual visual review
      until pixel-diff baseline exists." Deferred to manual visual
      review per atlas; the automated path verifies the resulting
      `Cluster (<threshold>)` column instead (the column is the
      observable, persisted output of the cut).
    verdict_status: SCOPE_REDUCTION
related_bugs: []
realized_as:
  - assign-clusters-spec.ts
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-migrate-01
    timestamp: 2026-06-03T00:00:00Z
    review_round: 1
    failure_keys: []
  d:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-migrate-01
    timestamp: 2026-06-03T00:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T13:18:00Z
    spec_runs:
      - spec: assign-clusters-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 120
        failure_keys: []
---

# Dendrogram | Assign Clusters end-to-end (column creation, two-way Threshold/Clusters binding, replace-on-rerun)

End-to-end scenario for atlas critical_path
`dendrogram.cp.assign-clusters-column-creation` (`priority: p1`) and
top-level interaction
`dendrogram.cross.cluster-assignment-roundtrip` (`coverage_type:
regression`): build a dendrogram via Chem | Analyze | Hierarchical
Clustering on the mol1K molecule column, open the Assign Clusters
dialog through both surfaces (canvas right-click context menu and the
magic-wand icon in the dendrogram toolbar), exercise the two-way
Threshold ↔ Clusters binding, click Assign, and verify that a fresh
`Cluster (<threshold>)` categorical column is appended to the host
DataFrame (and that re-running Assign with a different cluster count
appends a second auto-uniquely-named column rather than overwriting
the first). Closes with the cleanup/robustness path: re-running
Hierarchical Clustering on the same TableView warns and replaces the
existing dendrogram neighbor (no duplicate viewers).

This scenario complements:

- `hierarchical-clustering-chem.md` (section ui-smoke) — owns the
  dialog gateway surface (Distance/Linkage dropdown enumeration,
  Features column input, OK button, Chem top-menu entry).
- `hierarchical-clustering-chem-api.md` (apitest, source-matrix) —
  owns the 14-combo `{euclidean, manhattan} × {7 linkages}`
  parameter matrix on the molecule path.

This scenario therefore re-exercises the dialog gateway only as a
prerequisite to building the tree; its specialty coverage is
everything *after* the dendrogram appears: Assign Clusters,
column creation, two-way binding, second-column auto-uniqueness,
and replace-on-rerun.

## Setup

- Server: any environment with Chem, Dendrogram packages installed
  and the `System:AppData/Chem/mol1K.csv` sample dataset available.
- Datasets referenced: `System:AppData/Chem/mol1K.csv`.
- No upstream scenario state required — the scenario opens its own
  dataset and builds its own dendrogram (chain YAML
  `depends_on: []`, `produces:
  ["dendrogram-from-mol1k-chem-euclidean-ward",
  "cluster-column-from-assign"]`).
- Fresh browser session preferred (the cleanup-replace step in the
  Replace-on-rerun scenario asserts that only one dendrogram
  neighbor remains attached at the end of the run).

## Scenarios

### Scenario 1 — Build dendrogram via Chem | Analyze | Hierarchical Clustering on mol1K

1. Open the **mol1K** dataset (`System:AppData/Chem/mol1K.csv`).
   Wait for the table view to render and the `molecule` column to
   display structures.
   - Verification: a TableView is open on the loaded DataFrame; the
     `molecule` column is visible and renders structures (MOLECULE
     semType detected). No console errors.

2. From the top menu run **Chem | Analyze | Hierarchical
   Clustering...**.
   - Verification: the **Hierarchical Clustering** dialog opens with
     four inputs — **Table**, **Features**, **Distance**, **Linkage**.
     **Distance** defaults to `euclidean`, **Linkage** defaults to
     `ward`. The **Features** column-list defaults to the `molecule`
     column (MOLECULE semType auto-selection on the Chem menu entry).
     No console errors.

3. Keep **Features** = `molecule`, **Distance** = `euclidean`,
   **Linkage** = `ward`. Click **OK**.
   - Verification: a transient progress indicator (`Creating
     dendrogram ...`) is shown, then a dendrogram neighbor is
     injected to the left of the grid. Row heights, current row,
     hover, selection, and filter state are synchronized between
     the grid and the tree (Atlas sub_feature
     `dendrogram.clustering.inject-tree-for-grid` declares the
     bidirectional current/mouseOver/selection sync wiring). No
     console errors.

### Scenario 2 — Open Assign Clusters dialog via both surfaces (context menu + magic-wand icon)

Preconditions: Scenario 1 has run; the dendrogram neighbor is
attached to the grid.

4. Right-click anywhere on the dendrogram canvas and select
   **Assign Clusters** from the context menu.
   - Verification: the **Assign Clusters** dialog opens with a
     **Threshold** input (slider, defaulted to roughly half the
     tree height) and a **Clusters** input (integer, minimum 1).
     No console errors.

5. Close the dialog. Click the **magic wand** icon in the top-left
   corner of the dendrogram neighbor.
   - Verification: the same **Assign Clusters** dialog opens. The
     icon's tooltip reads `Assign Clusters`.

### Scenario 3 — Two-way Threshold ↔ Clusters binding

Preconditions: the **Assign Clusters** dialog is open from
Scenario 2 (re-open via either surface if it was closed).

6. Adjust the **Threshold** slider to a lower value.
   - Verification: the **Clusters** value automatically recalculates
     to match the new threshold; the two inputs are interconnected
     (Atlas sub_feature `dendrogram.clustering.assign-clusters-dialog`
     declares the Threshold→Clusters direction of the binding).

7. Type a specific value into the **Clusters** input (e.g. `5`).
   - Verification: the **Threshold** value automatically
     recalculates to the cut position that yields the requested
     number of clusters, or as close as possible per atlas
     edge_case "Cluster-assignment dialog Threshold/Clusters
     two-way binding — binary search runs 20 iterations to map a
     target Clusters count back to a Threshold; ambiguous (multiple
     thresholds match) falls back to closest count via minDiff"
     (inject-tree-for-grid2.ts#L394). See Notes for the tolerance
     ambiguity.

### Scenario 4 — Assign creates a `Cluster (<threshold>)` column on the host DataFrame

Preconditions: Scenarios 1-3 have run; the **Assign Clusters**
dialog is open with a known Clusters value (e.g. `5`).

8. Click **Assign**.
   - Verification: the dialog closes; a new string column named
     `Cluster (<threshold>)` (e.g. `Cluster (12.34)`, threshold
     rounded to 2 decimals) is added to the grid. Every row is
     labelled with its cluster id (`1`, `2`, …). The threshold
     value is preserved in the column name for reference. Atlas
     sub_feature `dendrogram.api.tree-helper.cut-tree-to-grid`
     declares the column-creation surface. No console errors.

9. Re-open **Assign Clusters** (via either surface), change
   **Clusters** to a different value, and click **Assign** again.
   - Verification: a second cluster column with a different
     `Cluster (<threshold>)` name is added; the previous column is
     preserved; the new column name is auto-incremented to stay
     unique (the prior cluster column is not overwritten).

### Scenario 5 — Replace-on-rerun cleanup (no duplicate viewers)

Preconditions: Scenario 1 has run; the dendrogram neighbor is
attached to the grid (Scenarios 2-4 may also have run — the
cluster columns added do not affect the rerun).

10. Close the dendrogram neighbor, then re-run **Chem | Analyze |
    Hierarchical Clustering...** on the same table with the same
    features (Features = `molecule`, Distance = `euclidean`,
    Linkage = `ward`). Click **OK**.
    - Verification: a warning notifies that the existing dendrogram
      is being closed and replaced (per atlas edge_case "Dendrogram
      already attached to a grid — re-running hierarchicalClusteringUI
      on the same TableView closes the existing GridNeighbor and
      warns before injecting a new one (DENDROGRAM_NEIGHBOR_TEMP_NAME
      tracking)", hierarchical-clustering.ts#L104). A fresh dendrogram
      neighbor is injected. No console errors. No duplicate dendrogram
      neighbors remain attached at the end of the step.

## Notes

- **Provenance.** Source Block B is annotated in the original TestTrack
  text as "community post #37"; source Block C is annotated as
  "community post #38". These are upstream community-feedback
  references for the Assign Clusters dialog (#37) and the visual cut
  indicator + horizontal-zoom surface (#38). The annotations have been
  moved here from the inline step text per source-text-fix
  `community-post-annotations-moved-to-notes`.
- **Step 7 tolerance ambiguity.** The original step asserts that
  Threshold recalculates "to the cut position that yields exactly —
  or as close as possible to — the requested number of clusters",
  but does not pin a tolerance for the automated assertion. Atlas
  edge_case on inject-tree-for-grid2.ts#L394 documents the 20-iteration
  binary search with `minDiff` fallback for ambiguous cuts; the
  automated test should reference that to constrain the accepted
  Threshold range (`unresolved_ambiguities:
  step-7-clusters-to-threshold-tolerance-not-pinned`).
- **Steps 10-13 deferred to manual.** Steps 10 (canvas-pixel
  visualisation of the cut position), 11 (plain-wheel vertical
  scroll), 12 (Ctrl/Cmd+wheel horizontal zoom with clamp), and 13
  (double-click empty area resets horizontal zoom) are deferred to
  manual visual review per atlas `manual_only[]` entries
  `dendrogram.mo.tree-canvas-visual-regression` and
  `dendrogram.mo.ctrl-wheel-zoom-tactile` — see frontmatter
  `scope_reductions: [SR-01, SR-02]` and the rationale fields there.
  The automated path verifies the column-creation observable
  (Scenarios 4-5) rather than the canvas pixel.
- **Max horizontal zoom factor ambiguity.** Original Step 12 says the
  zoom is "capped at a maximum factor" but does not specify the
  factor. The sibling `hierarchical-clustering-bio.md` (Step 9)
  cites a `1×–100×` clamp range; either harmonise with that range
  or cite the renderer constant directly when the manual review
  scenario is authored (`unresolved_ambiguities:
  step-12-max-horizontal-zoom-factor-not-specified`).
- **Section chain context.** This scenario is the integration owner
  of post-build dendrogram flows (chain
  `pyramid_layer: integration`); the dialog gateway smoke is owned
  by `hierarchical-clustering-chem.md` (section ui-smoke) and the
  parameter matrix by `hierarchical-clustering-chem-api.md`
  (source-matrix, apitest layer). The chain YAML's
  `ui_coverage_responsibility` for this scenario lists the
  Assign Clusters surface (dialog, magic-wand, threshold slider,
  clusters input, assign button), the dendrogram context menu, and
  the Ctrl+wheel/plain-wheel/double-click zoom surface; the zoom
  surface items are realized as the SR-01 deferral above.
- **Sub_features cited per scenario.** Scenario 1 →
  `dendrogram.clustering.menu.chem`,
  `dendrogram.clustering.dialog`,
  `dendrogram.clustering.inject-tree-for-grid`; Scenarios 2-3 →
  `dendrogram.event.context-menu`,
  `dendrogram.clustering.inject-tree-for-grid` (magic-wand icon),
  `dendrogram.clustering.assign-clusters-dialog`; Scenario 4 →
  `dendrogram.api.tree-helper.cut-tree-to-grid`; Scenario 5 →
  `dendrogram.clustering.menu.chem` (re-run),
  `dendrogram.clustering.inject-tree-for-grid` (replace warning
  surface, hierarchical-clustering.ts#L104).
- **Bug library.** `bug-library/dendrogram.yaml` lists one entry,
  GROK-13041 ("Dendrogram needs reset after filtering"), which
  affects `dendrogram.clustering.inject-tree-for-grid` /
  `dendrogram.event.selection-changed`. That bug's repro is owned
  by a separate filter-vs-sort regression scenario (atlas
  `edge_cases[dendrogram.ec.filter-does-not-trigger-remove-revert-prompt]`
  scope), not by this scenario; `related_bugs: []` in the
  frontmatter is therefore intentional.

---
{
  "order": 1,
  "datasets": ["System:AppData/Chem/mol1K.csv"]
}
