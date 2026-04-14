### Pivot table

1. Open the **demog** dataset.
2. On the **Toolbox → Viewers** section, click the **Pivot table** icon.
* **Expected result**: **Pivot table** viewer opens without errors in the console. By default, **Group by** = `DIS_POP`, **Pivot** = `SEVERITY`, **Aggregate** = `AGE(avg)`.
3. Close the viewer (`×` on its title bar), then re-add it the same way.
* **Expected result**: Viewer reopens without errors; defaults are the same.
4. Directly on the viewer, modify each tag-editor row:
   * Click `+` in **Group by**, pick `SEX`; click `×` on `DIS_POP` tag to remove it.
   * Click `+` in **Pivot**, pick `RACE`; remove the previous tag.
   * Click `+` in **Aggregate**, pick `HEIGHT` with `avg`, then add `WEIGHT` with `sum`.
   * Right-click a measure tag → **Aggregation → min** to change its aggregation.
* **Expected result**: Grid updates after each change with no console errors.
5. On the viewer, click the **Gear** icon → the **Property Pane** opens.
6. In the Property Pane, toggle **Show Header**, **Show Command Bar**, and **Filtering Enabled**; switch **Row Source** between `Filtered` / `All` / `Selected`.
* **Expected result**: Each change is reflected in the viewer without errors. `Show Header` hides the three tag-editor rows; `Show Command Bar` hides the history/refresh icons.
7. Title and layout persistence:
   * Open **SPGI** dataset and add the **Pivot table** viewer.
   * Set **Title** in the Property Pane, then edit it inline in the viewer header.
   * Apply column coloring via the inner grid's column header context menu (**Color coding → Linear/Categorical**).
   * Save the layout (**View → Layout → Save to Gallery**).
   * Switch to another layout, then re-apply the saved one.
* **Expected result**: Title and coloring persist across re-apply.
8. On **demog** with a Pivot table:
   * Apply coloring to one aggregate column.
   * In the Property Pane switch **Row Source** `Filtered` → `Selected` → `Filtered`.
* **Expected result**: Coloring is preserved across row-source changes.
9. Two-way property sync:
   * Change tags directly in the viewer → the corresponding `groupByColumnNames` / `pivotColumnNames` / `aggregateColumnNames` in the Property Pane update.
   * Change those lists in the Property Pane → tags in the viewer update.
10. Push pivot result to workspace:
    * On **demog** with a Pivot table, click the **ADD** button in the counts area (top of the viewer, under the tag rows).
* **Expected result**: A new table opens with the aggregated data as its own dataframe.


---
{
  "order": 12,
  "datasets": ["System:DemoFiles/demog.csv"]
}
