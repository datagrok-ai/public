### Dendrogram: Assign Clusters and horizontal zoom

#### Block A — Build the dendrogram (Hierarchical Clustering)

1. Open the **mol1K** dataset (`System:AppData/Chem/mol1K.csv`). Wait for the table view to render and the `molecule` column to display structures.
2. Run **Chem > Analyze > Hierarchical Clustering...** from the top menu.
* Expected result: the **Hierarchical Clustering** dialog opens with **Table**, **Features**, **Distance** (default *euclidean*), and **Linkage** (default *ward*) inputs. The `molecule` column is the default selected feature. No console errors.
3. Keep **Features** = `molecule`, **Distance** = *euclidean*, **Linkage** = *ward*. Click **OK**.
* Expected result: a progress indicator (`Creating dendrogram ...`) appears, then a dendrogram viewer is injected to the left of the grid. Row heights, current row, hover, selection, and filter state are synchronized between the grid and the tree. No console errors.

#### Block B — Assign Clusters dialog (community post #37)

4. Right-click anywhere on the dendrogram and select **Assign Clusters** from the context menu.
* Expected result: the **Assign Clusters** dialog opens with a **Threshold** input (slider, defaulted to roughly half the tree height) and a **Clusters** input (integer, minimum 1). No console errors.
5. Close the dialog. Now click the **magic wand** icon in the top-left corner of the dendrogram viewer.
* Expected result: the same **Assign Clusters** dialog opens (the icon tooltip reads *Assign Clusters*).
6. Adjust the **Threshold** slider to a lower value.
* Expected result: the **Clusters** value automatically recalculates to match the new threshold (the two inputs are interconnected).
7. Type a specific value into the **Clusters** input (e.g. `5`).
* Expected result: the **Threshold** value automatically recalculates to the cut position that yields exactly — or as close as possible to — the requested number of clusters.
8. Click **Assign**.
* Expected result: the dialog closes and a new string column named `Cluster (<threshold>)` (e.g. `Cluster (12.34)`, threshold rounded to 2 decimals) is added to the grid. Each row is labelled with its cluster id (`1`, `2`, ...). The threshold value is preserved in the column name for reference. No console errors.
9. Re-open **Assign Clusters**, change **Clusters** to a different value, and click **Assign** again.
* Expected result: a second cluster column with a different `Cluster (<threshold>)` name is added (the previous column is preserved, the name is auto-incremented to stay unique).

#### Block C — Visual split and horizontal zoom (community post #38)

10. Open the **Assign Clusters** dialog again and drag the **Threshold** slider / change the **Clusters** value.
* Expected result: the dendrogram shows exactly where the cut occurs — the split position is rendered on the tree and updates live as the threshold / cluster count changes.
11. Close the dialog. Hover over the dendrogram and scroll the **mouse wheel** without any modifier key.
* Expected result: the tree scrolls vertically; the X-axis scale does not change.
12. Hold **Ctrl** (or **Cmd** on macOS) and scroll the **mouse wheel** over the dendrogram.
* Expected result: the dendrogram zooms horizontally along the X-axis, spreading out clusters that were too close together. Zoom is clamped (no zoom-out below the default fit, capped at a maximum factor). Repeated Ctrl + wheel up keeps expanding the X-axis; Ctrl + wheel down contracts it back toward the default.
13. Double-click an empty area of the dendrogram (with no node under the cursor).
* Expected result: the horizontal zoom resets to the default fit.

#### Block D — Cleanup / robustness

14. Close the dendrogram viewer, then re-run **Chem > Analyze > Hierarchical Clustering...** on the same table with the same features.
* Expected result: a warning notifies that the existing dendrogram is being closed and replaced; a fresh dendrogram is injected. No console errors and no duplicate viewers left behind.

---
{
  "order": 1,
  "datasets": ["System:AppData/Chem/mol1K.csv"]
}
