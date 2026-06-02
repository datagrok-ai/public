### Bio: Hierarchical Clustering (sequences) — dialog UI smoke

#### Block A — Dialog exposes all Distance and Linkage values

1. Open the **FASTA_PT_activity** dataset. Wait for the `sequence` column to be detected as a Macromolecule (sequence renderer).
2. Run **Bio > Analyze > Hierarchical Clustering...** from the top menu.
* Expected result: the dialog opens with **Features** defaulting to the `sequence` macromolecule column, plus **Distance** and **Linkage** inputs.
3. Open the **Distance** dropdown.
* Expected result: exactly two values — **euclidean**, **manhattan** (default *euclidean*).
4. Open the **Linkage** dropdown.
* Expected result: exactly seven values — **single**, **complete**, **average**, **weighted**, **centroid**, **median**, **ward** (default *ward*).

#### Block B — Build the dendrogram from a sequence column (bio-specific path)

5. With **Features** = `sequence`, **Distance** = *euclidean*, **Linkage** = *ward*, click **OK**.
* Expected result: `Creating dendrogram ...` progress while the sequence distance matrix is computed (encode + Levenshtein), then a dendrogram is injected to the left of the grid with one leaf per row. Grid ↔ tree hover / current / selection / filter states are synchronized. No console errors and no `Unsupported column type`.
6. Close the dendrogram. Re-open the dialog, set **Distance** = *manhattan*, **Linkage** = *complete*, **Features** = `sequence`, click **OK**.
* Expected result: a dendrogram builds successfully (different shape, one leaf per row). No fatal console errors. (centroid/median may render non-monotonic trees — expected, not a defect.)

#### Block C — Shared UI smoke on a bio-built tree (delegated coverage)

7. Click the **magic wand** icon (top-left of the dendrogram) — or right-click the tree and choose **Assign Clusters**.
* Expected result: the **Assign Clusters** dialog opens with interconnected **Threshold** (slider) and **Clusters** (integer ≥ 1) inputs.
8. Set **Clusters** to `5` and click **Assign**.
* Expected result: the dialog closes and a new string column `Cluster (<threshold>)` is added, labelling each row with its cluster id. No console errors.
9. Hold **Ctrl** (or **Cmd**) and scroll the **mouse wheel** over the dendrogram.
* Expected result: the tree zooms horizontally along the X-axis (clamped 1×–100×). Plain wheel (no modifier) scrolls vertically.

## Notes

- Bio path verified against `Dendrogram/src/package.ts:322` (`Bio | Analyze | Hierarchical Clustering...` → default column `bySemType(MACROMOLECULE)`) and `tree-helper.ts:538-542` (encodeSequences + LEVENSHTEIN).
- All three menu entries — `Bio | Analyze`, `Chem | Analyze`, `ML | Cluster` — invoke the same `hierarchicalClusteringDialog`; they differ only in the default selected feature column (MACROMOLECULE / MOLECULE / either).

---
{
  "order": 4,
  "datasets": ["System:AppData/Bio/samples/FASTA_PT_activity.csv"]
}
