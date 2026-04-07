#### Sunburst viewer

1. **Open files and initialize the viewer**
   - Open **SPGI_v2.csv** and **demog.csv**.
   - Click **Add viewer** → **Sunburst** viewer for each file.
   **Expected**: Sunburst viewer should open without errors for both files.

2. **Viewer properties**
   - Click the **Gear** icon in the Sunburst viewer.
   **Expected**: The context panel with viewer properties should open.

3. **Viewer properties functionality**

   3.1. **Table switching**
   - Switch between **SPGI_v2.csv** and **demog.csv**.
   **Expected**: The viewer updates without errors.

   3.2. **Hierarchy configuration**
   - Open the **Select Columns** dialog → Choose 2–4 columns → Click OK.
   **Expected**: The hierarchy should update based on the selected columns.
   - Reopen the dialog → Search for a column → Click Cancel.
   **Expected**: No changes should be applied.
   - Verify that the text in the viewer is fully visible (if there is enough space) or hidden (if there is not enough space), with a tooltip showing the value on hover.
   - Ensure that structural columns render at the correct size.

   3.3. **Inherit from grid**
   - In **demog.csv**, open the **Select Columns** dialog and choose the **SEX** column.
   - Enable the **Inherit from grid** property.
   - Apply categorical coloring to **SEX** column in the grid.
   **Expected**: The Sunburst viewer should reflect the grid colors.
   - Change the coloring.
   **Expected**: The viewer should update accordingly.

   3.4. **Include nulls**
   - In **SPGI_v2.csv**, select the **Core** and **R101** columns.
   - Enable **Include nulls** property.
   **Expected**: Grey segments should appear for null values.
   - Disable the **Include nulls** option.
   **Expected**: The grey segments should disappear.

4. **View reset**
   - Double-click on empty space or use the context menu → **Reset View**.
   **Expected**: The view should reset to its initial state.

5. **Multi-selection behavior**
   - Click → Single segment selected.
   - Ctrl + Click → Multi-select.
   - Ctrl + Shift + Click → Deselect segment.
   **Expected**: Grid rows should update accordingly for each action.

6. **Select/filter on empty category**
   - Open **SPGI_v2.csv**.
   - Select a column with nulls (e.g., **Sampling Time**).
   - Click on the null segment (grey).
   **Expected**: The segment should behave like any other category (selects/filters the respective rows).

7. **Projects & layouts**
   - Configure the Sunburst viewer with 3–4 columns on **SPGI_v2** → Save the project.
   - Close all → Reopen the project.
   **Expected**: The viewer should restore correctly.
   - Save the layout and apply it.
   **Expected**: The viewer should retain the saved layout.

8. **Old layout compatibility**
   - Open the layout from issue [#2979](https://github.com/datagrok-ai/public/issues/2979).
   **Expected**: The viewer should show the selected columns, and the columns in the dialog should be in sync.

9. **Collaborative filtering**
   - In **demog.csv**, configure the Sunburst viewer and apply both internal and panel filters.
   **Expected**: The filters should combine correctly, and the viewer should reflect the intersection of the filters.

---
{
  "order": 29,
  "datasets": ["System:DemoFiles/demog.csv"]
}