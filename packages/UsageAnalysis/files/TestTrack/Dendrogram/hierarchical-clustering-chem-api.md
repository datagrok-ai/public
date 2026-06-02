---
feature: dendrogram
sub_features_covered:
  - dendrogram.hierarchical-clustering.distance.euclidean
  - dendrogram.hierarchical-clustering.distance.manhattan
  - dendrogram.hierarchical-clustering.linkage.single
  - dendrogram.hierarchical-clustering.linkage.complete
  - dendrogram.hierarchical-clustering.linkage.average
  - dendrogram.hierarchical-clustering.linkage.weighted
  - dendrogram.hierarchical-clustering.linkage.centroid
  - dendrogram.hierarchical-clustering.linkage.median
  - dendrogram.hierarchical-clustering.linkage.ward
  - dendrogram.hierarchical-clustering.molecule-path
  - dendrogram.hierarchical-clustering.numeric-path
target_layer: apitest
coverage_type: regression
pyramid_layer: integration
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Dendrogram/hierarchical-clustering-chem-api.md
related_bugs: []
---

# Hierarchical Clustering (chem) — Distance × Linkage matrix (JS API)

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
