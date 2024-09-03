1. Open SPGI, SPGI-linked1, SPGI-linked2.
2. Go to **Tables > SPGI.**
3. On the **Viewers tab,** click **Bar chart.**
4. On the **Bar chart** viewer, click the **Gear** icon. The **Property Pane** opens.
5. Split and coloring:
    1. Go to the **Category** info panel.
    2. Set Split to `Series`, `Primary Scaffold Name`
    3. Check/uncheck the checkboxes
    4. Go to the grid. Right-click the Primary Scaffold Name column and select **Color coding > Categorical, Edit.**
    5. Edit the color scheme. Check the bar chart coloring.
    6. Save to Layout. Check
6. Check aggregation functions
7. Selection
    1. Go to the grid and select the first 300 rows. The selection should be reflected on bar chart
    2. Go to the Property Pane > Selection. Toggle the **Show selected rows** checkbox.
    3. Click the bar chart and check the selection in grid.
8. Filtering
    1. Go to **Tables** and click the **Filter** icon
    2. Change arbitrary settings on the **Filter Panel** and check the bar chart interaction.
9. Go to the **Data** info panel and check all the properties, including the following:
    1. Tables switching (SPGI-linked1, SPGI-linked2, SPGI)
    2. **Row Source** with different options
    3. Set **Filter** to `${CAST Idea ID} < 636500`, set color coding by Chemical Space Y, and arbitrarily change other options.
    4. Save to Layout. Check
    5. Set **On click** to **Filter**. Check filtering (to reset filtering, simply click on any empty space around the bar chart)
10. Title, Description:
    1. Add a title and a description.
    2. Change their position
    3. Check the range slider functionality
    4. Save to Layout. Check
11. Scrolling:
    1. Change the vertical range slider and use it to scroll
    2. Set **Value** to CAST Idea ID. Scroll
12. Tooltip
    1. Right-click on the bar chart and verify that all options under the **Tooltip** tab are functioning correctly.
    2. Save to Layout. Check
13. Inverting axes:
  * Precondition: opened Bar chart viewer on the SPGI dataset.
  * Open the Bar Chart's Properties panel.
  * Navigate to Misc > Orientation. Change the viewer's orientation to Auto, Horizontal, and Vertical.
  * Expected result: The Bar Chart viewer should rotate the axes according to the applied settings. No errors or unexpected behavior should occur during the orientation changes.
---
{
  "order": 3,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}
