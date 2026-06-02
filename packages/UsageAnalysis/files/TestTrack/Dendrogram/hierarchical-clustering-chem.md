### Hierarchical Clustering (chem) — dialog UI smoke

#### Block A — Dialog exposes all Distance and Linkage values

1. Open the **mol1K** dataset.  Wait for the `molecule` column to render structures.
2. Run **Chem > Analyze > Hierarchical Clustering...**.
* Expected result: the dialog opens with **Table** = mol1K, **Features** defaulting to `molecule`, **Distance**, and **Linkage** inputs.
3. Open the **Distance** dropdown.
* Expected result: exactly two values are listed — **euclidean**, **manhattan** (default *euclidean*).
4. Open the **Linkage** dropdown.
* Expected result: exactly seven values are listed, in order — **single**, **complete**, **average**, **weighted**, **centroid**, **median**, **ward** (default *ward*).

#### Block B — Representative end-to-end runs (spot-check)

5. With **Features** = `molecule`, **Distance** = *euclidean*, **Linkage** = *ward*, click **OK**.
* Expected result: `Creating dendrogram ...` progress, then a dendrogram is injected to the left of the grid with one leaf per row. No console errors.
6. Close the dendrogram. Re-open the dialog, set **Distance** = *manhattan*, **Linkage** = *single*, **Features** = `molecule`, click **OK**.
* Expected result: a dendrogram builds successfully (different shape, still one leaf per row). No `Unsupported column type` and no fatal console errors.
7. Close the dendrogram. Re-open the dialog, set **Distance** = *euclidean*, **Linkage** = *centroid* (or *median*), **Features** = numeric columns (`pIC50_HIV_Integrase`, `Q`), click **OK**.
* Expected result: a dendrogram builds without error. **Note:** centroid/median linkage may render a **non-monotonic** tree (branches with inversions / apparent cross-overs) — this is mathematically expected for these methods and is **not** a defect.

## Notes

- All three menu entries — `Bio | Analyze`, `Chem | Analyze`, `ML | Cluster` — call the same
  `hierarchicalClusteringDialog`; they differ only in the default-selected feature column.
- The `molecule` distance path calls `Chem:getMorganFingerprints` → Tanimoto; numeric columns
  use the per-column difference metric (`tree-helper.ts:536-548`).
- Existing package unit coverage (`Dendrogram/src/tests/hierarchical-clustering-tests.ts`)
  exercises only `euclidean` + `average` on a numeric column — the other distances/linkages
  and the molecule path are not covered there.

---
{
  "order": 2,
  "datasets": ["System:AppData/Chem/mol1K.csv"]
}
