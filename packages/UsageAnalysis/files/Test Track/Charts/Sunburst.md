#### Sunburst viewer (Charts package)

1. Open SPGI_v2.csv. On the menu ribbon, click **Add viewer** and select **Sunburst viewer**. 
- The Sunburst viewer should open for SPGI_v2.csv.
2. Open demog.csv. On the menu ribbon, click **Add viewer** and select **Sunburst viewer**. 
- The Sunburst viewer should open for demog.csv.
3. **Open Viewer Properties:**
- In the **Sunburst viewer** viewer, click the **Gear** icon. 
- The **Context Panel** opens with viewer`s properties. 
4. **Check all properties**. The following features should work as expected:
  - a. Switching **tables**:
    - Switching between tables should be smooth and reflected in the viewer.
  - b. **Hierarchy**:
    - Open the **Select columns dialog** and choose 2–4 columns. Click **OK**. The Sunburst should update to reflect the selected columns.
    - Open the dialog again, search for a column, make changes, and click **Cancel**. No changes should be applied.
    - Text on the viewer should be fully visible or hidden, with a tooltip showing the full value on hover.
    - Structural columns should render with correct size and alignment.
  - c. **Inherit from the Grid**:
    - Open demog.csv. In the Select columns dialog, choose the SEX column.
    - Check the Inherit from the Grid option.
    - In the grid, set categorical coloring for the **SEX** column.
    - The Sunburst should reflect the same colors for SEX.
    - Update colors in the grid. The Sunburst should also update accordingly.
  - d. **Include Nulls**:
    - Open SPGI data. Select Core and R101 columns.
    - Enable the Include Nulls checkbox.
    - Two grey "null" segments should appear in the Sunburst.
    - Disable the checkbox. Grey segments should disappear.
5. **Nulls rendering**
- Open SPGI data.
- Select **Idea Author** and **Series** columns. Enable Include Nulls check-box.
- In the filter panel, filter **Series** by null values.
- Expected result: two rings should appear – green (Idea Author) and grey (Series). Viewer should still render if all values in a selected column are null.
6. **Reset view**
- Double-clicking on empty space in the Sunburst viewer should reset the view.
- Context menu > selecting "Reset View" from the context menu should reset the view.
7. **Multiple selection**
- Click a segment to select its rows in the grid. 
- Clicking multiple sectors should allow multi-selection: values for all picked sectors should be selected on the grid:
  - Click – select a single segment
  - Ctrl + Click – select multiple segments
  - Ctrl + Shift + Click – deselects the one you clicked on
8. **Select/filter on empty category**
- Open SPGI data.
- Use a column with nulls (e.g., **Sampling Time**) in Sunburst.
- Click on the empty (null) category.
- Expected result: it behaves like any other category – filters or selects the respective rows.
9. **Projects & Layouts**
- Configure the Sunburst on SPGI with 3–4 columns.
- Save the project. Close all, reopen, and check that Sunburst is restored.
- Save the layout. Switch to another layout, then back. Sunburst viewer should persist.
10. **Old Layout Compatibility**
- Open a table and apply a layout from issue [2979](https://github.com/datagrok-ai/public/issues/2979). 
- Open the Sunburst properties.
- Expected result: viewer shows selected columns, and Select columns dialog reflects the selection.
11. **Collaborative filtering**
- Close previous data. Open demog.csv and add Sunburst.
- Configure internal filtering and add filters in the filter panel.
- Expected result: Sunburst respects filter panel constraints.
- Internally, a bitset from Sunburst is combined with the one from the filter panel using a logical AND. The result should be a correct intersection view.




---
{
  "order": 29,
  "datasets": ["System:DemoFiles/demog.csv"]
}