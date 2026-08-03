### Network diagram

1. Open the **demog** dataset.
2. On the **Viewers** tab of the Toolbox, click **Network diagram**. The viewer opens with
   two column combo boxes in the top-left corner (**Node 1**, **Node 2**) auto-picked by the
   viewer (first two categorical columns with the fewest categories, e.g. `sex` and `race`).
   * Expected: Nodes and edges are drawn on the canvas; no errors in the console.
3. In the **Node 1** combo box, switch the column to `RACE`; in **Node 2**, switch to `DEMOG`.
   * Expected: The graph rebuilds with new endpoints.
4. Click a node.
   * Expected: Rows backing every incident edge are selected (check the status bar /
     `grok.shell.tv.dataFrame.selection.trueCount`).
5. Shift+click another node, then Ctrl+click the first node.
   * Expected: Shift+click adds to selection; Ctrl+click toggles it.
6. Click an edge.
   * Expected: Only the rows backing that edge are selected.
7. Double-click empty canvas.
   * Expected: Selection clears on the first click; the view zooms to fit on the second.
8. Click the **Gear** icon on the viewer title bar. The **Property Pane** opens.
9. In the **Data** section, set:
   * **Edge Color Column Name** = `age`, **Edge Color Aggr Type** = `avg`.
   * **Edge Width Column Name** = `weight`, **Edge Width Aggr Type** = `avg`.
   * **Node 1 Size Column Name** = `age`, **Node 1 Color Column Name** = `sex`.
   * Expected: Edge colors follow a gradient, edge widths vary, Node 1 nodes resize and
     recolor accordingly.
10. In the **Style** section, toggle:
    * **Show Column Selectors** off → the two top-left combo boxes disappear; on → they
      reappear.
    * **Show Arrows** = `to` → directional arrows appear on edges.
    * **Suspend Simulation** = true → node positions freeze; false → physics resumes.
11. Filter the dataframe (e.g. `age > 40` in the filter panel).
    * Expected: Filtered-out nodes disappear. Toggle **Show Filtered Out Nodes** = true →
      they reappear in the filtered-out color.
12. Close the viewer via the **×** icon on its title bar.
    * Expected: The viewer is removed without errors.

---
{
  "order": 23,
  "datasets": ["System:DemoFiles/demog.csv"]
}
