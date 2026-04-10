# PC plot tests (Playwright)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add PC plot

## Menu Ribbon and To Script

1. On the Menu Ribbon, click the **Add viewer** icon and select **PC plot** -- viewer should open
2. Right-click the plot > **To Script** -- a balloon with the generated script should appear
3. Close the viewer, then add PC plot via **Toolbox** -- viewer opens again

## Axis scale & normalization

1. Open Context Panel > **Value > Normalize Each Column** -- enabled by default
2. Verify each axis has its own min/max scale
3. Disable **Normalize Each Column** -- all axes should use the same global scale, a shared Y axis appears on the left
4. Re-enable -- individual scales restored, shared Y axis disappears
5. Right-click on the plot > **Y Axis** > select **Global** -- axes switch to global scale
6. Right-click > **Y Axis** > select **Normalized** -- axes return to normalized scale
7. Open Context Panel > **Value > Log Columns** > add AGE -- the AGE axis should switch to logarithmic scale
8. Add WEIGHT to log columns -- WEIGHT axis also switches to log
9. Remove both from log columns -- linear scales restored

## Selection & line display

1. Open Context Panel > **Selection**
2. Disable **Show Current Line** -- current row line should disappear
3. Re-enable **Show Current Line** -- current row line reappears
4. Disable **Show Mouse Over Line** -- hover line should disappear
5. Re-enable **Show Mouse Over Line** -- hover line reappears
6. Enable **Show Mouse Over Row Group**
7. Disable **Show All Lines** -- only current, mouse-over, and selected lines should be visible
8. Re-enable **Show All Lines** -- all lines visible again
9. Right-click > **Selection** > toggle **Show Current Line** -- same effect as Context Panel
10. Right-click > **Selection** > toggle **Show All Lines** -- same effect as Context Panel

## Style & layout

1. Open Context Panel > **Style**
2. Change **Line Width** to 3 -- lines become thicker
3. Change **Current Line Width** to 5 -- current row line becomes prominent
4. Change **Mouse Over Line Width** to 5 -- hover line becomes thicker
5. Right-click > **Lines** > adjust **Line Width** slider -- line thickness changes (same property as step 2)
6. Change **Labels Orientation** to Vert -- axis labels rotate
7. Change **Min Max Orientation** to Vert -- min/max values rotate
8. Adjust **Horz Margin** -- spacing between axes changes
9. Toggle **Auto Layout** off -- manual control over labels, min/max, and margins

## In-chart filtering & reset

1. On the first axis, drag the range slider handles to narrow the range -- rows outside range should be filtered
2. Narrow a second axis range -- only rows passing both filters should remain
3. Enable **Show Filtered Out Lines** -- filtered rows should appear as faint lines
4. Double-click white space -- all range sliders reset (Reset View)
5. Narrow range sliders on two axes again
6. Right-click > **Reset View** -- all filters should reset
7. Right-click > **Filter** > toggle **Show Filters** -- range slider visibility toggles

## Filter panel interaction

Precondition: range sliders on two axes are narrowed (from scenario 4 flow).

1. Use range sliders on PC plot axes to filter data
2. Open the Filter Panel and apply an additional filter
3. Verify filtering respects both Filter Panel and PC plot filters
4. On the PC plot, right-click > **Reset View** -- only the PC plot filtering should reset, Filter Panel filters remain
5. Use range sliders on the PC plot axes to filter again
6. On the Filter Panel click **Reset filter** -- all filters should reset and PC plot range sliders return to default

## Column management & reordering

1. Open Context Panel > **Value > Column Names** -- verify AGE, HEIGHT, WEIGHT are shown
2. Remove HEIGHT from the list -- the corresponding axis should disappear
3. Add HEIGHT back -- axis reappears
4. Reorder columns in the list -- axes should reorder on the plot
5. Hover on AGE column name on top (above column range slider), drag'n'drop it to change viewer column order
6. Verify the new order persists in Context Panel > **Value > Column Names**

## Density styles

1. Open Context Panel > **Box Plot > Density Style** -- default is **circles**
2. Switch to **box plot** -- box plot overlays should appear on each axis
3. Toggle **Show Interquartile Range** off and on
4. Toggle **Show Upper Dash** off and on
5. Toggle **Show Lower Dash** off and on
6. Toggle **Show Mean Cross** off and on
7. Toggle **Show Median** off and on
8. Enable **Show Circles** -- data point circles appear alongside the box plot
9. Switch to **violin plot** -- violin shape overlays should appear
10. Adjust **Bins** slider -- violin resolution should change
11. Adjust **Whisker Line Width** -- whisker thickness changes
12. Switch back to **circles**

## Color coding, legend & grid coloring

**Color column basics:**

1. Open Context Panel > **Color > Color Column** > set to AGE -- lines should be color-coded by age
2. Change **Color Axis Type** to **log** -- color mapping changes
3. Toggle **Invert Color Scheme** -- gradient reverses
4. Set custom **Color Min** and **Color Max** values -- color range narrows

**Color scale and legend:**

5. Set Color Column to RACE -- a categorical legend should appear
6. Change **Legend Position** to Left, Right, Top, Bottom -- legend moves accordingly
7. Set **Legend Visibility** to Never -- legend disappears
8. Set **Legend Visibility** to Auto -- legend reappears
9. Set Color Column to None -- all lines return to default color

**Color coding from grid:**

10. Open Context Panel > **Color > Color Column** > set to HEIGHT
11. In the grid, for the HEIGHT column set linear color coding type (use js-api) -- verify that PC plot legend updates accordingly
12. In the grid, for the HEIGHT column set conditional color coding type (use js-api) -- verify the PC plot legend updates accordingly

## Title and description

1. Open Context Panel > set **Title** to "My PC Plot" -- title should appear on the viewer
2. Set **Description** to "Test description" -- description text should appear
3. Change **Description Position** -- description moves to the specified position
4. Clear both fields -- title and description disappear

## Pick Up / Apply

1. Add a second PC plot
2. On the first PC plot, change the set of axes (e.g. remove HEIGHT)
3. Switch AGE axis to log scale
4. Set color column to RACE and move legend to the left
5. Set title to "Source Plot"
6. Right-click the first PC plot > **Pick up/Apply > Pick up**
7. Right-click the second PC plot > **Pick up/Apply > Apply** -- second plot should match the first (axes, log, color, legend, title)
8. Change the axes on the first PC plot -- the second plot should not be affected
9. Adjust the range slider on the second PC plot -- the first plot should update to show filtered lines, but its own range sliders should remain unchanged

## Layout and project save/restore

1. Save the layout
2. Open a scatterplot
3. Apply the saved layout -- verify only the original viewers are displayed
4. Save the project (use js-api)
5. Close All
6. Open the saved project - verify pc plots are displayed


## Table switching and transformation

1. Open spgi-100.csv, open demog
2. Add a PC plot
3. Go to Context Panel > Data and set **Table** to spgi-100
4. In the Data section click the Transformation input field and enter [{"#type":"GroupAggregation","aggType":"key","colName":"Chemist 521"},{"#type":"GroupAggregation","aggType":"pivot","colName":"Series"},{"#type":"GroupAggregation","aggType":"count","colName":"Id"}]
6. Verify the PC plot updates to show the pivoted aggregated data
7. Close All

---
{
  "order": 7,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
