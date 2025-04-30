1. Open SPGI, SPGI-linked1, SPGI-linked2.
2. Go to SPGI and add a line chart
3. Click the **Gear** icon and open the **Context Panel**

#### 1. Context Panel > Data

1. Switch Tables: SPGI-linked1, SPGI-linked2, and back to SPGI
2. Toggle the checkboxes
3. Save the Layout. Check
4. Set **Filter** to `${CAST Idea ID} <636500`
5. Set **Row Source** to `Selected`, select rows in the grid and check the line chart interaction
6. Save the Layout. Check

#### 2. Split updates in properties

1. Add or remove split columns directly from the viewer - **Check the number of the Split columns on the Context Panel**
1. Line chart: **Context menu > Properties > Data > Split** and set some new Split columns  - **Check the number of the Split columns on the Context Panel and on the viewer**

#### 3. Filtering

1. Set Split to `Primary Series Name`,`Series` and `Core` - **verify that legend colors are distinct and match the colors on the chart**
1. Open the **Filter Panel**
1. Apply structure, categorical, and numeric filters — **verify correct collaborative filtering between the Filter Panel and in-viewer Filter value**
1. For the **Primary Series Name** filter, toggle the checkboxes on and off — **legend colors should remain consistent with the corresponding line colors**
1. Use **Ctrl+click** to select/deselect  categories in the legend - **check the line chart filtering **
4. Save the Layout. Check

#### 4. Selection:

1. On the line chart, use **Shift+Mouse Drag** to select points. Check the selection in the grid
2. Go to the **Context Panel** > **Selection**
1. Toggle the checkboxes. Check the line chart interaction.
3. Save to Layout. Check

#### 5. Scrolling & Zooming:

1. Use the **Mouse Wheel** to zoom in and out the line chart
1. Set **Split by** `Series` and `Core`
2. Right-click the line chart and select **Overview > Area Chart/Stacked Bar Chart**
3. Use the Overview section to scroll and zoom in/out the line chart
4. Save to Layout. Check

#### 6. Title, Description:

1. Add a title and description.
2. Change their position
3. Check the placement of the legend, title, and description
3. Check the range slider functionality
4. Save to Layout. Check

#### 7. Formula lines

1. Add a new line chart
1. Set Split to `Stereo Category`
1. Set **X** to `Chemical Space X`
1. Set **Y** to `Average Mass` and set the **log scale**
1. Add a formula line:  `${Average Mass} = 0.75* ${Chemical Space X}* ${Chemical Space X} - 4 * ${Chemical Space X} +300`
1. Change Color and Style
1. Switch **Y Axis Type** to linear/log - **check the Line on the viewer**
1. Grid: change the `Chemical Space X` name - **check the formula line caption**

####  8. Line chart: Custom tooltip [#2357](https://github.com/datagrok-ai/public/issues/2357)

1. Set X to `Chemist 521`
1. **Context Menu > Tooltip > Edit** - a dialog opens
1. Add categorical, numerical and date columns - **the available aggregations should match the column types**
1. Hover over the line chart and check the tooltip

####  9. Line chart with multiple axes

1. Add a new line chart
1. Open the Context Panel
1. Remove Y lines using the 'x' buttons on the viewer - **the number of Y columns on the Context Panel should be updated immediately**
1. Set two Y axis
1. Enable **Data > Multi Axis** 
1. Set some Split columns - **verify that legend categories' names start with the corresponding Y column **
1.  Use the in-plot column selector to change one of the Y columns - **the legend should be updated accordingly**

#### 10. [GROK-17835](https://reddata.atlassian.net/browse/GROK-17835)
1. Close All
1. Open SPGI
2. Add a line chart
1. Set **Data > Multi axis** On
1. Set Split to `Series` and `Scaffold Names`
1. Switch X axis to `Chemist 521`
1. Hover over the viewer - check, no errors occure
---  
{
  "order": 5,
  "datasets": ["System:DemoFiles/SPGI.csv","System:DemoFiles/SPGI-linked1.csv","System:DemoFiles/SPGI-linked2.csv"]
}