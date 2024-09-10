1. Open SPGI, SPGI-linked1, SPGI-linked2.
2. Go to tab with SPGI dataset
3. On the **Toolbox,** click **Histogram**. Use tooltips to find Histogram if you do not know how it looks like.
4. On the **Histogram** viewer, click the **Gear** icon. The **Context Panel** opens.
5. Click the **Toggle filters** icon. Filter rows in dataset with true in 'Has Unlableled R-groups'
6. Go to the **Data** info panel and check all the properties, including the following:
    1. Tables switching (SPGI-linked1, SPGI-linked2, SPGI)
    2. Toggle the checkbox: Show Filtered out Rows. The viewer should change its data based on the checkbox.
    3. Try to filter rows for viewer: `${TPSA} > 60`
7. Layouts:
    1. On Toolbox Click Save on Layout accordion.
    2. Right-click on sidebar > Close all
    3. Open SPGI dataset
    4. On Toolbox open Layouts accordion, wait several seconds
    5. Click on the loaded layout to apply it to opened dataset
    6. Check that everything from the configured data from previous step is applied correctly
8. Filtering
    1. Click the **Toggle filters** icon
    2. Change arbitrary settings on the **Filter Panel** and check the histogram interaction.
    3. Enable the **Filtering Enabled** option by navigating to Context Pane > Properties > Data.
    4. Change the horizontal range slider and use it to scroll and data filtering. Dataset should filter accordingly.
    5. Disable the **Filtering Enabled** option by navigating to Context Pane > Properties > Data.
    6. Dataset should not be filtered.
9. Selection
    1. Go to the grid and select the first 50 columns. The selection should be reflected on the histogram.
    2. Go to the Context Pane > Selection. Toggle the checkboxes:
        1. Show Current row: green dot on the X axis should indicate current selected row
        2. Show Mouse over Row: gray dot on the X axis should indicate current hovered row
        3. Show Mouse over Row Group:
            1. Add Bar Chart viewer
            2. Check that on hovering histogram viewer change the appearance.
    3. Click the histogram data and check the selection in grid.
    4. Change Row Source to selected and check Histogram
10. Title, Description:
    1. Context Panel > Description section
    2. Add a title and a description.
    3. Change their position
    4. Check the range slider functionality
11. Layouts:
     1. On Toolbox Click Save on Layout accordion.
     2. Right-click on sidebar > Close all
     3. Open SPGI dataset
     4. On Toolbox open Layouts accordion, wait several seconds
     5. Click on the loaded layout to apply it to opened dataset
     6. Check that everything from the configured data from previous step is applied correctly

---
{
"order": 4,
"datasets": ["System:DemoFiles/SPGI.csv"]
}
