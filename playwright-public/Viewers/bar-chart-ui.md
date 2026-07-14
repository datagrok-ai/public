# Bar chart tests (manual)

Manual checklist. Not included in PW automation.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Bar chart

## Selection interactions

1. Set **On Click** to **Select**
2. Click a bar — rows in that category are selected (highlighted in grid)
3. **Ctrl + click** another bar — selection adds that category
4. **Shift + drag** a rectangle over multiple bars — all covered categories are selected
5. Click empty space — selection clears

## Stacked bar filtering

1. Set Stack to SEX, **On Click** to **Filter**
2. Click on a stacked segment (e.g., Male within a RACE bar) — only matching rows remain
3. **Ctrl + click** another segment — that category is added to the filter
4. Click empty space — filter resets, all data returns

## Zoom and pan

1. **Alt + drag** vertically on the chart — chart should zoom into that range of categories
2. Use **mouse wheel** to scroll through categories when zoomed in
3. **Drag** on empty space (no modifiers) — chart should pan vertically through categories
4. **Double-click** the chart — view should reset to show all categories

## Scrolling (range sliders)

1. Change the vertical range slider and use it to scroll
2. Set **Value** to CAST Idea ID (open SPGI dataset). Scroll through categories using the range slider

## Tooltip configuration

1. Right-click on the bar chart and open the **Tooltip** tab
2. Verify that all tooltip options are functioning correctly
3. Hover over different bars — tooltip content should match configuration

## Split and coloring (color scheme editing)

1. Set Split to a categorical column (e.g., RACE)
2. Go to the grid. Right-click the RACE column and select **Color coding > Categorical, Edit**
3. Edit the color scheme — bar chart colors should update to match
4. Save to Layout. Close and reopen — colors should persist

## Context menu regions

1. Set Stack to SEX
2. Right-click on the **value axis** — axis-specific menu should appear (Axis Type, Include Nulls, Axis Font)
3. Close the menu
4. Right-click on the **chart area** — full context menu should appear (Reset View, Orientation, Data, Order, Controls, Selection)
5. Close the menu
6. Right-click on the **legend** — legend-specific menu should appear (Legend Visibility, Legend Position)

## Grid color coding sync with bar chart

Setup: Close all, open SPGI

1. Add a bar chart
2. Set Split to Series and Primary Scaffold Name
3. Go to the grid, right-click the header of Primary Scaffold Name and select Color coding > Categorical, Edit
4. Edit the color scheme — verify the bar chart coloring updates to match
5. Save the layout, close the bar chart, apply the layout — verify color coding is preserved

## Selection interaction and Show Selected Rows

1. Go to the grid and select the first 300 rows — verify the selection is reflected on the bar chart
2. Go to the Context Panel > Selection, toggle Show Selected Rows off — verify the bar chart updates
3. Toggle Show Selected Rows on — verify selected rows are highlighted again
4. Click a bar on the bar chart — verify corresponding rows are selected in the grid

## On Click filter mode

1. Go to the Context Panel > Data, set On Click to Filter
2. Click on a bar — verify the grid filters to show only that category's rows
3. Ctrl+click another bar — verify filter updates
4. Click on empty space around the bar chart — verify filter resets and all rows are visible
5. Set On Click back to Select

---
{
  "order": 101,
  "datasets": ["System:DemoFiles/demog.csv,System:DemoFiles/SPGI.csv"]
}
