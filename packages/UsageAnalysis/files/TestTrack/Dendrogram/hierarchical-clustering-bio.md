---
feature: dendrogram
target_layer: playwright
coverage_type: regression
priority: p0
realizes: [dendrogram.cp.hier-clustering-bio-sequence-path]
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

Bio-path integration scenario: verify the Bio top-menu surface
auto-selects a MACROMOLECULE (sequence) column as the default
Features input on the Hierarchical Clustering dialog, then build a
dendrogram from the sequence column through the
Levenshtein-on-encoded-sequences distance path. Scenario 4 is an
explicit ride-along smoke on the post-build Assign Clusters surface;
the dialog gateway (Distance/Linkage dropdown enumeration, OK button,
menu entry) is owned by `hierarchical-clustering-chem.md`, and the
full Assign Clusters + Ctrl+wheel zoom surface is owned by
`assign-clusters.md`. This scenario re-exercises the gateway only as
a precondition to building the bio-path tree; its specialty coverage
is everything bio-path-specific (the Bio menu entry, MACROMOLECULE
auto-select, Levenshtein distance, sequence-column leaf binding).

This scenario complements:

- `hierarchical-clustering-chem.md` — owns the shared
  hierarchical-clustering dialog gateway (Distance/Linkage dropdowns,
  Features column input, OK button, menu entry surface for the
  Chem-default path).
- `assign-clusters.md` — owns the Assign Clusters dialog full surface
  (two-way Threshold/Clusters binding, column creation,
  replace-on-rerun, Ctrl+wheel zoom deferral).
- `hierarchical-clustering-bio-api.md` — owns the 14-combo
  `{euclidean, manhattan} × {7 linkages}` parameter matrix on the
  sequence path and the MACROMOLECULE semType precondition assertion.

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
     distance matrix is computed via the encode + Levenshtein path,
     then a dendrogram neighbor is injected to the left of the grid
     with one leaf per DataFrame row. Grid ↔ tree current /
     mouseOver / selection / filter states are synchronized. No
     console errors and no `Unsupported column type` error.

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
     (integer, minimum 1). The two inputs are interconnected. No
     console errors.

8. Set **Clusters** to `5` and click **Assign**.
   - Verification: the dialog closes and a new categorical column
     `Cluster (<threshold>)` is added to the host DataFrame (one
     row per original sequence row; the column labels each row
     with its cluster id). No console errors.

## Notes

- **Bio menu path provenance.** The `Bio | Analyze | Hierarchical
  Clustering...` top-menu entry is registered by the Dendrogram
  package (function `hierarchicalClusteringSequences`), not by the
  Bio package itself. The registering function seeds the dialog's
  Features input via `bySemType(MACROMOLECULE)`, which produces the
  sequence-default behavior exercised in Scenario 1. The dialog
  itself is identical across all three top-menu entries
  (`Bio | Analyze`, `Chem | Analyze`, `ML | Cluster`) — they differ
  only in the default-selected feature column. The shared dialog
  gateway is owned by `hierarchical-clustering-chem.md`.
- **Levenshtein distance path.** For MACROMOLECULE columns,
  `TreeHelper.calcDistanceMatrix` encodes the sequence values and
  computes pairwise distances via Levenshtein. Numeric columns use
  the per-column difference metric; MOLECULE columns use Tanimoto on
  Morgan fingerprints (Chem path). The sequence path is exercised
  here in Scenarios 3 and 4; the corresponding parameter matrix (14
  distance × linkage combinations on the sequence path) is asserted
  by `hierarchical-clustering-bio-api.md`.
- **Centroid/median note.** Step 6 uses `complete` linkage, not
  `centroid` or `median`. For reference: centroid and median
  linkages may produce a non-monotonic tree (branches with
  inversions or apparent cross-overs) — this is mathematically
  expected for these linkage methods, not a defect, and is asserted
  separately by the apitest matrix in
  `hierarchical-clustering-bio-api.md`. Open question: whether to
  add a dedicated centroid/median step here to mirror the
  chem-dialog pattern in `hierarchical-clustering-chem.md` Step 7.
- **Ctrl+wheel zoom deferred.** Holding Ctrl/Cmd and scrolling over
  the dendrogram to zoom horizontally (clamped 1×–100×; plain wheel
  scrolls vertically) is deferred to manual visual review —
  wheel/zoom behavior depends on browser-level event delivery.
  `assign-clusters.md` is the canonical owner of this surface.
- **Bug library.** GROK-13041 ("Dendrogram needs reset after
  filtering") affects the same dendrogram/grid sync surface, but its
  repro lives in a separate filter-vs-sort regression scenario
  (`grok-13041-filter-no-prompt.md`), not here.

---
{
  "order": 4,
  "datasets": ["System:AppData/Bio/samples/FASTA_PT_activity.csv"]
}
