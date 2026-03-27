1. Open spgi-100
1. Open the Filter Panel

#### 1. Filtering and resetting

1. Hover over the Structure filter - the Skretch icon appers - click it to open the Sketcher
1. In the Sketcher insert c1ccccc1 to the text input field  - the molecule appears in the sketcher
1. Click Ok - the filter should apply (32 rows are filtered) 
1. HOver over the Filter Panel and scroll down to find Stereo Category filter
1. On the Stereo Category filter, click R_ONE category - the filter should apply (~15 row is filtered)
1. Hover over the Filter Panel and scroll down to find the Average Mass filter
1. Scroll down to **Average Mass**, click filter menu icon and click Min/max item - min and max input fields apear. Click the max input field and enter 400. The filter applies (~4 rows shown).
5. Hover over the Filter Panel — icons appear at the top — uncheck the **Turn filters on/off** checkbox. All filters grey out and all rows are shown.
5. Enable the **Turn filters on/off** checkbox — the three filters reactivate and the ~298-row result is restored.
7. Hover over the **Stereo Category** filter - its icons appear - disable the **Stereo Category** filter using its **Turn on/off** checkbox — verify the row count increases to reflect only Structure + Average Mass filters.
8. Close the **Filter Panel**  - all rows should be present in the dataset
8. Reopen the Filter Panel — verify that the **Stereo Category** filter is still disabled and the row count ~1204.
9. Click the **Reset** icon in the Filter Panel header and confirm the prompt — all filters are cleared and all rows are shown.
10. Close the **Filter Panel**.
11. Reopen the **Filter Panel** — verify that no filters are present (panel is in its default state).
12. To remove the Structure filter, hover over it - its icons appear - click the “X” button
13. Hover over the Core filter and click its the “X” button, to remove the Core filter
7. Close and reopen the **Filter Panel** — **verify that the removed filters are no longer present**

#### 2. Filter Panel UI Behavior

**Precondition:** Filter Panel is open, no active filters (state after Test 1 reset).

1. Open a **Scatter plot** and a **Bar chart** viewers.
2. Go to the Filter Panel
3. Hover over a value in the Stereo Category filter — verify that the corresponding data points are highlighted in the Scatter plot and Bar chart.
4. Hover over the Competition assay Date filter — verify that the tooltip displays the date format matching the column's current format.
5. In the grid, right-click the Competition assay Date column header - the menu appears - hover over the Format item - menu expands - click the `MM/dd/yyyy` format
5. Go back to the FIlter Panel  and hover over Competition assay Date filter — verify that the tooltip in the date filter shows `MM/dd/yyyy` format.
6. Hover over the Series filter — the filter icons appear — click the **Search categories** icon. The search input field appears inside the filter.
6. Go to the Toolbox on the left
7. To save the current layout, on the Toolbox expand the Layouts pane and click SAVE - the saved layout will appear right under the Save button
8. Click the saved layout (it will be the first layout right under the Save button) and verify that the search field is still visible in the Series filter.
9. Reset filters

#### 3. Adding and Reordering Filters

**Precondition:** Filter Panel is open, no active filters (state after Test 2 reset).

1. Open the **Hamburger menu** in the Filter Panel header and select **Reorder Filters** — the dialog opens
2. In the Drag to Reorder dialog, click Structure and click the 'Move to the bottom' icon - the Structure should be removed to the bottom of the list
2. In the Drag to Reorder dialog, click Core and click the 'Move down' icon - Core should be the second in the list.
2. Click CLOSE in the bottom of the dialog
2. Verify that in the Filter Panel Core is the second filter and Structure is the last
3. On the Filter Panel, open **Hamburger menu > Select Columns** — the dialog opens
3. In the Select Columns dialog, check the checkbox next to ID - verify that ID filter card appears in the Filter Panel.
3. In the Select Columns dialog, uncheck the checkbox next to ID and click OK — verify that ID filter card is not shown in the Filter Panel.
4. On the Filter Panel, open **Hamburger menu > Remove All** — verify all filter cards are removed from the panel.
5. Add new filters one by one using each of the following methods and verify that each new filter appears **at the top** of the Filter Panel:
   - Drag the ID column header from the grid into the Filter Panel.
   - In the grid, for the CAST Idea ID column, click hamburger menu in the header and in the Filter section click Add Filter.
   - In the grid, right-click the first cell in the Structure column - the context menu opens - hover over the Current Value item and then click Use as filter 
   - In the Filter Panel, right-click the background - the context menu opens - hover over the Add filter and then click Scaffold Tree Filter.
5. verify that the order of the filters is the following:  Scaffold Tree Filter, Structure filter,   CAST Idea ID, ID
6. On the Filter Panel, open Hamburger menu and click Remove All — verify all filter cards are removed.
7. Close the **Filter Panel** — verify it closes without errors.

#### 4. Hidden Columns

**Precondition:** Filter Panel is open, no active filters (state after Test 3 reset).

1. Right-click the first grid cell - the context menu opens - click the Order Or Hide Columns... item - the Order Or Hide Columns dialog opens
1. In the Order Or Hide Columns dialog uncheck the checkbox next to Structure, Core and R1 and click OK
1. Verify the grid doesn't show Structure, Core and R1 columns.
2. Open the **Filter Panel** — verify that there are no Structure, Core and R1 filter cards on the Filter Panel.
3. Right-click a grid cell - the context menu opens - and click the Order Or Hide Columns... item - the Order Or Hide Columns dialog opens — check the checkboxes next to Structure, Core and R1 — close the dialog. Verify the grid shows all columns again.
4. Go to the Filter Panel and click the **Hamburger menu > Remove All** to clear filters, then close and reopen the **Filter Panel** — verify that filter cards for Structure, Core and R1 are now available in the panel.

---
{
"order": 1,
"datasets": ["System:DemoFiles/spgi-100.csv"]
}
