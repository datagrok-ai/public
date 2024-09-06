1. Open SPGI, SPGI-linked1, SPGI-linked2.
2. Go to **Tables > SPGI.**
3. On the **Viewers tab,** click **Line chart.**
4. On the **Line chart** viewer, click the **Gear** icon. The **Property Pane** opens.
5. Go to the **Data** info panel and check all the properties, including the following:
    1. Tables switching (SPGI-linked1, SPGI-linked2, SPGI)
    2. Toggle the checkboxes
    3. Save to Layout. Check
    4. Set **Filter** to ${CAST Idea ID} <636500
    5. Set **Row Source** to Selected, select rows in the grid and check the line chart interaction
    6. Save to Layout. Check
6. Filtering
    1. Go to **Tables** and click the **Filter** icon
    2. Change the Chemical Space X tab’s range slider and check the line chart interaction.
    3. Draw a molecule in the Structure tab and check the line chart interaction.
    4. Save to Layout. Check
7. Selection:
    1. On the line chart, use **Shift+Mouse Drag** to select points. Check the selection in the grid.
    2. Go to the viewer’s Property Pane > Selection. Toggle the checkboxes. Check the line chart interaction.
    3. Save to Layout. Check
8. Scrolling & Zooming:
    1. Use the **Mouse Wheel** to zoom in and out the line chart
    2. Right-click the line chart and select **Overview > Area Chart** from the context menu
    3. Use the area chart as range sliders and scroll and zoom the line chart.
    4. Save to Layout. Check
9. Title, Description:
    1. Add a title and a description.
    2. Change their position
    3. Check the range slider functionality
    4. Save to Layout. Check
10. Formula lines
11. Line chart: Custom tooltip [#2357](https://github.com/datagrok-ai/public/issues/2357)
12. [#2623](https://github.com/datagrok-ai/public/issues/2623) Line chart: implement the one-click way to set the Split by columns
13. **Color legend**:
  * In opened line chart viewer with two or more splits set. Inspect the color legend displayed in the viewer:
  * Expected result: Visually confirm that the colors in the legend match the colors on the corresponding lines in the viewer. The color legend should be clear and distinguishable, with no two lines sharing the same color unless intended.
14. **Color-coding while setting filters**:
  * Open the SPGI dataset. Open the Line Chart viewer. Double-click on any line in the viewer to open a detailed line graph.
  * **Configure Viewer**:
    * Splits: Set "Primary Series Name" and "Series".
    * X-Axis: Set to "Stereo Category".
    * Y-Axis: Set to (avg) "Average Mass".
  * Apply Filters:
    * Open Filters Pane. Find "Primary Series Name" filter. 
    * Toggle the check-boxes for the "Primary Series Name" filter on and off.
  * **Expected Result**: Colors in the legend should match consistently with the colors of the corresponding lines in the viewer, both before and after applying filters.
15. **Split updates in properties**. 
  * Test on the SPGI dataset with opened new Line Chart viewer: 
    * Open the Line Chart's properties panel.
    * Add or remove split columns directly from the Line Chart viewer.
    * Check the number of selected columns in the properties panel.
  * **Expected results**: number of selected columns in properties should be updated immediately on adding or removing a column.
16. **Line chart with multiple axes**:
* Precondition: Open the **SPGI** dataset with **Line chart** viewer.
* Configure the Line Chart:
  * Select two Y columns and enable the "Multi Axis" option.
  * Add a split column to the chart. Ensure that a legend appears on the plot.
  * Change Y Column: Use the in-plot column selector to change one of the Y columns.
* Expected Results: The legend should be updated accordingly when the Y column is changed, reflecting the new data associated with the selected column.
17. **Y Columns Count Update in Properties Panel**:
* Precondition: Open a **Line chart** viewer with **multiple Y columns**. Properties panel should remain opened. 
* Remove Y Columns. On the plot itself, click the 'x' buttons to remove some of the Y columns.
* Expected Results: The number of Y columns in the properties panel should update immediately when columns are removed directly from the plot.
18. Use Ctrl+click to select/deselect one more category in the legend. [#2455](https://github.com/datagrok-ai/public/issues/2455) 

---  
{
  "order": 5,
  "datasets": ["System:DemoFiles/SPGI.csv","System:DemoFiles/SPGI-linked1.csv","System:DemoFiles/SPGI-linked2.csv"]
}