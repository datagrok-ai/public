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
  - dendrogram.hierarchical-clustering.sequence-path
target_layer: apitest
coverage_type: regression
pyramid_layer: integration
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Dendrogram/hierarchical-clustering-bio-api.md
related_bugs: []
---

# Hierarchical Clustering (bio) — Distance × Linkage matrix (JS API)

Dataset: **FASTA_PT_activity** (`System:AppData/Bio/samples/FASTA_PT_activity.csv`, ~99
rows). 

Matrix:
- **Distance**: `euclidean`, `manhattan`.
- **Linkage**: `single`, `complete`, `average`, `weighted`, `centroid`, `median`, `ward`.
- **Feature path**: the `sequence` column (MACROMOLECULE semType → encode + Levenshtein).

## Setup

1. Authenticate via `spec-login.ts` (`loginToDatagrok`).
2. Load `System:AppData/Bio/samples/FASTA_PT_activity.csv` and
   `await df.meta.detectSemanticTypes()` so `sequence` is recognized as MACROMOLECULE.

## Scenarios

### Scenario 1: sequence path — 14 combos build a valid tree

For each `(distance, linkage)` in `{euclidean, manhattan} × {single, complete, average,
weighted, centroid, median, ward}` (14 combos):

1. Compute hierarchical clustering on the `sequence` column via the Dendrogram compute path
   (`getTreeHelper()` → `calcDistanceMatrix(df, ['sequence'], distance)` →
   `getClusterMatrixWorker(matrix.data, rowCount, linkageCode)` →
   `parseClusterMatrix(clusterMatrix)`), where `linkageCode` is the enum index for `linkage`.
2. **Verify:**
   - the call resolves **without throwing** (the macromolecule branch runs — no
     `Unsupported column type`);
   - the parsed Newick root is non-null and its **leaf count === row count** (every input
     sequence appears exactly once as a leaf);
   - no fatal console error during the combo.

### Scenario 2: sequence encoding precondition

1. **Verify:** `df.col('sequence').semType === DG.SEMTYPE.MACROMOLECULE` after semantic-type
   detection (the precondition that routes `calcDistanceMatrix` into the encode + Levenshtein
   branch rather than failing as unsupported).

## Notes

- Structural assertion only ("valid tree over all rows, no throw") — exact Levenshtein
  distances / merge heights belong in library-level unit tests
  (`@datagrok-libraries/ml` macromolecule distance functions).
- centroid/median may yield non-monotonic merge heights; acceptable here — the assertion is
  leaf-completeness + no-throw.
- The linkage-code positional-mapping contract (enum order → worker code) is asserted once in
  the chem `-api.md`; not duplicated here.

## Dataset metadata

```json
{
  "order": 5,
  "datasets": ["System:AppData/Bio/samples/FASTA_PT_activity.csv"]
}
```
