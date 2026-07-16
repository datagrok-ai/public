---
feature: dendrogram
target_layer: apitest
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: []
pyramid_layer: integration
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Dendrogram/hierarchical-clustering-chem-api.md
related_bugs: []
realized_as:
  - hierarchical-clustering-chem-api-spec.ts
scope_reductions:
  - id: SR-FRONTMATTER-ATLAS-RESOLUTION
    rationale: |
      Original Test-Designer-authored sub_features_covered list used synthetic
      matrix-enumeration ids (dendrogram.hierarchical-clustering.distance.{euclidean,
      manhattan}, .linkage.{single,complete,average,weighted,centroid,median,ward},
      .molecule-path, .numeric-path) which do NOT resolve in
      feature-atlas/dendrogram.yaml :: sub_features[].id. E-STRUCT-MECH-06 failed
      Gate E with failure_keys: [E-STRUCT-MECH-06] (cycle 2026-06-03-dendrogram-
      automate-02). Remapped to the 5 atlas-resolvable ids that name the actual
      JS-API surfaces the spec exercises: the registered Dendrogram:hierarchical-
      Clustering function (dendrogram.clustering.api), the Dendrogram:getTreeHelper
      function (dendrogram.api.get-tree-helper), calcDistanceMatrix
      (dendrogram.api.tree-helper.calc-distance-matrix), parseClusterMatrix
      (dendrogram.api.tree-helper.parse-cluster-matrix), and the GridNeighbor
      mount that observability-verifies injectTreeForGrid2
      (dendrogram.clustering.inject-tree-for-grid). The distance × linkage matrix
      coverage is preserved at the SPEC level (the spec iterates all 14
      distance×linkage combos for both molecule and numeric paths); the atlas
      side does not enumerate per-distance / per-linkage as separate ids, so
      matrix breadth is asserted within the spec body rather than as separate
      atlas anchors. Upstream resolution path: Test Designer either (a) extend
      the atlas with matrix-enumeration sub_features or (b) accept the
      compute-surface atlas anchors as canonical and remove the matrix-enum
      naming convention from future scenarios.
gate_verdicts:
  e:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T00:00:00Z
    failure_keys: []
  b:
    verdict: FAIL
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T14:35:00Z
    spec_runs:
      - spec: hierarchical-clustering-chem-api-spec.ts
        result: failed
        attempts: 3
        duration_seconds: 0
        failure_keys: [B-COLLECT-ABORT, B-STAB-01]
---

# Hierarchical Clustering (chem) — Distance × Linkage matrix (JS API)

Verifies that hierarchical clustering builds a valid tree for every
supported Distance × Linkage combination on the Chem (molecule) and
numeric paths, using the JS API directly (no UI).

Dataset: **mol1K** (`System:AppData/Chem/mol1K.csv`). For runtime cost, load the table and
clone a small slice (e.g. first 60 non-null rows) so the 14 molecule combos + numeric combos
stay fast.

Matrix:
- **Distance**: `euclidean`, `manhattan` (`DistanceMetric` enum).
- **Linkage**: `single`, `complete`, `average`, `weighted`, `centroid`, `median`, `ward`
  (`LinkageMethod` enum).
- **Feature paths**: (a) `molecule` column (MOLECULE semType → Tanimoto), (b) numeric
  columns `['pIC50_HIV_Integrase', 'Q']` (FLOAT → difference).

## Setup

1. Authenticate via `spec-login.ts` (`loginToDatagrok`).
2. Load `System:AppData/Chem/mol1K.csv`, ensure semantic types are detected
   (`await df.meta.detectSemanticTypes()`) so `molecule` is recognized as MOLECULE.
3. Take a small deterministic row slice to bound runtime.

## Scenarios

### Scenario 1: molecule path — 14 combos build a valid tree

For each `(distance, linkage)` in `{euclidean, manhattan} × {single, complete, average,
weighted, centroid, median, ward}` (14 combos):

1. Compute hierarchical clustering on the `molecule` column for that combo via the
   Dendrogram compute path (`getTreeHelper()` → `calcDistanceMatrix(df, ['molecule'],
   distance)` → `getClusterMatrixWorker(matrix.data, rowCount, linkageCode)` →
   `parseClusterMatrix(clusterMatrix)`), where `linkageCode` is the enum index for `linkage`.
2. **Verify:**
   - the call resolves **without throwing** (no `Unsupported column type`);
   - the parsed Newick root is non-null and its **leaf count === sliced row count**
     (every input row appears exactly once as a leaf);
   - no fatal console error is recorded during the combo.

### Scenario 2: numeric path — 14 combos build a valid tree

Repeat Scenario 1 with `colNames = ['pIC50_HIV_Integrase', 'Q']` (numeric → difference
metric). Same three assertions per combo.

### Scenario 3: linkage-code mapping is positional and stable

1. Assert `Object.values(LinkageMethod)` equals
   `['single','complete','average','weighted','centroid','median','ward']` in this exact
   order (guards the `findIndex` contract that maps a linkage name to the worker code).
2. **Verify:** `ward` resolves to index 6 and `average` to index 2 (the index the existing
   WASM unit test pins).

## Notes

- Assertion is **structural** ("valid tree over all rows, no throw"), not a numeric cluster
  oracle — exact merge heights for each method belong in library-level unit tests
  (`@datagrok-libraries/math`, `/ml`). This scenario guards integration + the enum-order
  contract, which existing tests do not.
- centroid/median may yield non-monotonic merge heights; that is acceptable here — the
  assertion is leaf-completeness + no-throw, not height monotonicity.
- The bio (sequence/Levenshtein) matrix is covered separately in
  `hierarchical-clustering-bio-api.md`.

## Dataset metadata

```json
{
  "order": 3,
  "datasets": ["System:AppData/Chem/mol1K.csv"]
}
```
