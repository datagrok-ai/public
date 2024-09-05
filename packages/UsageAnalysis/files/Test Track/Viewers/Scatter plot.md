1. Open SPGI, SPGI-linked1, SPGI-linked2.
2. Go to **Tables > SPGI.**
3. On the **Viewers tab,** click **Scatter plot.**
4. On the **Scatter plot** viewer, click the **Gear** icon. The **Property Pane** opens.
5. Go to the **Data** info panel and check all the properties, including the following:
    1. Tables switching (SPGI-linked1, SPGI-linked2, SPGI)
    2. **Zoom and Filter** with different options (check when zooming in and out, combine with changing range sliders)
    3. Check "pack and zoom by filter"
    4. Save to Layout. Check
    5. Set **Filter** to ${CAST Idea ID} <636500 , set color coding by Series, and set different options for **Legend Visibility** and **Legend Position**. Arbitrary change some other settings
6. Go to the **X** and **Y** info panels and check all the properties, including the following:
    1. For min and max values for the axes, use different and equal values (636000 and 638000 for the X axis)
    2. Set the Axis Type to `logarithmic` (for numerical values)
    3. Invert the axes
    4. Save to Layout. Check
7. Coloring:
    1. Go to the grid.
    2. Right-click the header of the CAST Idea ID column and select **Color coding > Linear, Edit**. In the dialog, change the color scheme and invert it.
    3. Right-click the header of the Primary Series Name column and select **Color coding > Categorical**.
    4. Right-click the header of the Chemical Space X column and select **Color coding > Conditional, Edit**. In the dialog, set the **Apply to** setting to `text`.
    5. Go to the scatter plot viewer.
    6. Check color coding when setting **Color** to Primary Series Name, and Chemical Space X, CAST Idea ID .
    7. Go to the grid.
    8. Right-click the header of the CAST Idea ID column and select **Color coding > Edit**. In the dialog, change the color scheme and invert it. The color scheme on the scatter plot should change too.
    9. Save to Layout. Check
8. Legend (Checked in 5.c.), Title, Description:
    1. Add a title and a description.
    2. Change their position
    3. Check the range sliders functionality
    4. Save to Layout. Check
9. Selection (and scatter plot-specific selection options like lasso):
    1. On the scatter plot, use **Shift+Mouse Drag** to select points. Check the selection in the grid.
    2. Check the **Misc** > **Mouse Drag** (Pan or Select) setting
    3. Go to the scatter **plot > Property Pane > Misc**. Check the **Lasso Tool** checkbox. Use it on the scatter plot. Check the selection in the grid.
    4. Save to Layout. Check
    5. Select some points. Use Ctrl+Shift+Drag for deselection.
10. Tooltip:
    1. Right-click the scatter plot and check all options on the **Tooltip** tab.
    2. Save to Layout. Check
    3. Go to Property Pane > Tooltip. Check all options.
    4. Save to Layout. Check
11. Formula lines (dataframe lines, viewer-specific lines):
    1. Right-click the scatter plot and select **Tools > Formula Lines.** A **Formula Lines** editor opens. Click the ADD NEW button and select every option from the menu.
    2. Check the Viewer and Dataframe tabs
    3. Save to Layout. Check
    4. Delete created lines.
12. Using structures as axes:
    1. Change the X axis to Structure. Check the visibility of the structures.
    2. Save to Layout. Check
13. Interaction with Filter Panel:
    1. Go to **Tables** and click the **Filter** icon
    2. Right-click a cell in the `Structure` column
    3. Select **Current value > Use as filter** from the context menu. Check the **Filter Panel** and the scatter plot.
    4. Change arbitrary settings on the **Filter Panel** and check the scatter plot interaction.
    5. Save to Layout. Check
14. Verify Table Filtering Behavior with Empty Column on Log Scale Axis in Scatter Plot
* Open a dataset that includes at least one numerical column with no values (e.g., all entries in the column are null or NaN). For example, use fruits_.csv for the test.
* Add a scatter plot viewer to the dataset. Set the scatter plot’s zoom and filter setting to "Filter by zoom." Select the empty numerical column for one of the axes (X or Y). Change the axis to a log scale.
* Check the scatter plot: Since the axis is set to log scale and the column has no values, the scatter plot should display no data points.
* Check Table Filtering: Examine the table after configuring the scatter plot. 
Expected Result: The scatter plot should display no data points, but the table should remain unfiltered, displaying all original data rows. 
15. Verify that Scatterplot ignores negatives and zero values when switching to log scale ([#2456](https://github.com/datagrok-ai/public/issues/2456)):
- Switch axes to the log scale
- Try to set categorical axis to the X


---
{
  "order": 2,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}