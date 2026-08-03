1. Open spgi-100
2. Open the Filter Panel

#### 1. Filtering and resetting

1. Hover over the Structure filter - the Skretch icon appers - click it to open the Sketcher
2. In the Sketcher insert c1ccccc1 to the text input field  - the molecule appears in the sketcher
3. Click Ok - the filter should apply (32 rows are filtered) 
4. HOver over the Filter Panel and scroll down to find Stereo Category filter
5. On the Stereo Category filter, click R_ONE category - the filter should apply (~15 row is filtered)
6. Hover over the Filter Panel and scroll down to find the Average Mass filter
7. Scroll down to **Average Mass**, click filter menu icon and click Min/max item - min and max input fields apear. Click the max input field and enter 400. The filter applies (~4 rows shown).
8. Hover over the Filter Panel — icons appear at the top — uncheck the **Turn filters on/off** checkbox. All filters grey out and all rows are shown.
9. Enable the **Turn filters on/off** checkbox — the three filters reactivate and the ~298-row result is restored.
10. Hover over the **Stereo Category** filter - its icons appear - disable the **Stereo Category** filter using its **Turn on/off** checkbox — verify the row count increases to reflect only Structure + Average Mass filters.
11. Close the **Filter Panel**  - all rows should be present in the dataset
12. Reopen the Filter Panel — verify that the **Stereo Category** filter is still disabled 
13. Click the **Reset** icon in the Filter Panel header and confirm the prompt — all filters are cleared and all rows are shown.
14. Close the **Filter Panel**.
15. Reopen the **Filter Panel** — verify that no filters are present (panel is in its default state).
16. To remove the Structure filter, hover over it - its icons appear - click the "X" button
17. Hover over the Core filter and click its the "X" button, to remove the Core filter
18. Close and reopen the **Filter Panel** — **verify that the removed filters are no longer present**

#### 2. Adding Filters

**Precondition:** Filter Panel is open, no active filters (state after Test 2 reset).

1. On the Filter Panel, open **Hamburger menu > Remove All** — verify all filter cards are removed from the panel.
2. Add new filters one by one using each of the following methods and verify that each new filter appears **at the top** of the Filter Panel:
   - Drag the ID column header from the grid into the Filter Panel.
   - In the grid, for the CAST Idea ID column, click hamburger menu in the header and in the Filter section click Add Filter.
   - In the grid, right-click the first cell in the Structure column - the context menu opens - hover over the Current Value item and then click Use as filter 
   - In the Filter Panel, right-click the background - the context menu opens - hover over the Add filter and then click Scaffold Tree Filter.
3. verify that the order of the filters is the following:  Scaffold Tree Filter, Structure filter,   CAST Idea ID, ID
4. On the Filter Panel, open Hamburger menu and click Remove All — verify all filter cards are removed.
5. Close the **Filter Panel** — verify it closes without errors.

#### 3. Hidden Columns

**Precondition:** Filter Panel is open, no active filters (state after Test 3 reset).

1. Right-click the first grid cell - the context menu opens - click the Order Or Hide Columns... item - the Order Or Hide Columns dialog opens
2. In the Order Or Hide Columns dialog uncheck the checkbox next to Structure, Core and R1 and click OK
3. Verify the grid doesn't show Structure, Core and R1 columns.
4. Open the **Filter Panel** — verify that there are no Structure, Core and R1 filter cards on the Filter Panel.
5. Right-click a grid cell - the context menu opens - and click the Order Or Hide Columns... item - the Order Or Hide Columns dialog opens — check the checkboxes next to Structure, Core and R1 — close the dialog. Verify the grid shows all columns again.
6. Go to the Filter Panel and click the **Hamburger menu > Remove All**
7. Close the Filter Panel.
8. Open the **Filter Panel** — verify that filter cards for Structure, Core and R1 are now available in the panel.

---
{
"order": 1,
"datasets": ["System:AppData/Chem/tests/spgi-100.csv"]
}
