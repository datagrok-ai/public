---
feature: pcplot
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes:
  - viewers.pc-plot
realized_as:
  - pc-plot-spec.ts
related_bugs: []
expected_results:
  - anchor: "Menu Ribbon and To Script"
    expectation: "Right-click > To Script > To JavaScript raises a non-empty
      `.d4-balloon` whose text carries the generated viewer-creation script
      (contains `addViewer`). Closing then re-adding the viewer via the Toolbox
      icon reopens it."
  - anchor: "Axis scale via the context menu"
    expectation: "Y Axis > Global sets normalizeEachColumn to false; Y Axis >
      Normalized sets it back to true."
  - anchor: "Selection & line display"
    expectation: "Driving the current/mouse-over/all line props raises no error and
      the viewer stays alive; the context-menu Selection > Show Current Line and
      Show All Lines items flip the matching prop and restore it on a second click."
  - anchor: "Style & layout"
    expectation: "Driving line widths, label/min-max orientation, margins and Auto
      Layout raises no error; the axis sliders are still present after the layout
      pass. The Lines > Line Width menu slider writes the same lineWidth prop
      (documented reduction — the menu slider is a canvas control)."
  - anchor: "Reset and filter visibility from the context menu"
    expectation: "Narrowing two axis range sliders drops df.filter.trueCount
      progressively; Reset View restores it to the full count. Show Filters in the
      Filter submenu round-trips the showFilters prop (true -> false -> true)."
  - anchor: "Filter panel interaction"
    expectation: "A Filter Panel histogram filter and an in-chart range slider
      combine (AND): df.filter.trueCount drops below full, then below that.
      Reset View clears ONLY the in-chart part (count returns to the Filter-Panel
      value, the panel filter survives); the Filter Panel Reset filters button
      clears everything back to the full count."
  - anchor: "Column reordering from the Context Panel list"
    expectation: "Setting columnNames reorders the axis sliders; the rendered
      slider name order matches the requested column order."
  - anchor: "Density component toggles"
    expectation: "densityStyle defaults to 'circles'; driving every box-plot /
      violin component toggle raises no error and the viewer stays alive."
  - anchor: "Color coding, legend & grid coloring"
    expectation: "Categorical coloring renders a legend listing the column
      categories; Legend Visibility Never hides it and Auto restores the same
      categories. Conditional grid color-coding surfaces its bins ('20-150',
      '150-250') in the plot legend; switching to a linear/numeric scheme drops the
      DOM legend (the gradient is canvas)."
  - anchor: "Title and description"
    expectation: "The title renders in the panel titlebar (.panel-titlebar-text) and
      the description inside the viewer element; clearing both removes them."
  - anchor: "Pick Up / Apply"
    expectation: "Pick Up on the first plot then Apply on the second copies axes,
      log-scale columns, color column, legend position and title. Changing the
      first plot's axes afterwards does not affect the second; a range slider on
      the second plot filters the shared DataFrame (df.filter.trueCount drops)."
  - anchor: "Layout round-trip"
    expectation: "A saved layout, reloaded after adding a scatterplot, restores only
      the original viewers (PC plot present, scatterplot absent)."
  - anchor: "Project save / Close All / reopen"
    expectation: "Saving the view as a project through the ribbon Save button, then
      Close All and reopening the project, restores the PC plot and its title."
  - anchor: "Table switching and transformation"
    expectation: "The PC plot bound to spgi renders the raw numeric axes; applying a
      GroupAggregation pivot on Series replaces them with the Series categories
      (axes contain 'Triazoles'); clearing the Transformation reverts the axes to
      the original set."
---

# PC plot tests (Playwright)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add PC plot

## Menu Ribbon and To Script

1. On the Menu Ribbon, click the **Add viewer** icon and select **PC plot** -- viewer should open
2. Right-click the plot > **To Script** > **To JavaScript** -- a balloon with the generated script should appear
3. Close the viewer, then add PC plot via **Toolbox** -- viewer opens again

**Actuation note:** the spec adds the viewer through the Toolbox icon
(`[name="icon-pc-plot"]`) for both adds. The Menu Ribbon **Add viewer** gallery is
a canvas-rendered dialog with no headless handles, so it is not exercised; the
Toolbox path produces the same viewer.

## Axis scale via the context menu

1. Right-click on the plot > **Y Axis** > select **Global** -- axes switch to global scale
2. Right-click > **Y Axis** > select **Normalized** -- axes return to normalized scale

## Selection & line display

1. Open Context Panel > **Selection**
2. Disable **Show Current Line** -- current row line should disappear
3. Re-enable **Show Current Line** -- current row line reappears
4. Disable **Show Mouse Over Line** -- hover line should disappear
5. Re-enable **Show Mouse Over Line** -- hover line reappears
6. Enable **Show Mouse Over Row Group**
7. Disable **Show All Lines** -- only current, mouse-over, and selected lines should be visible
8. Re-enable **Show All Lines** -- all lines visible again
9. Right-click > **Selection** > toggle **Show Current Line** -- same effect as Context Panel (the showCurrentLine prop flips)
10. Right-click > **Selection** > toggle **Show All Lines** -- same effect as Context Panel (the showAllLines prop flips)

**Actuation note:** which lines get painted is a canvas outcome with no DOM
counterpart, so steps 2-8 are a no-error floor over the prop surface. Steps 9-10
assert the menu -> prop round-trip, the same state the Context Panel toggles write.

## Style & layout

1. Open Context Panel > **Style**
2. Change **Line Width** to 3 -- lines become thicker
3. Change **Current Line Width** to 5 -- current row line becomes prominent
4. Change **Mouse Over Line Width** to 5 -- hover line becomes thicker
5. Right-click > **Lines** > adjust **Line Width** slider -- line thickness changes

   **Reduction:** the Lines menu **Line Width** slider writes the same `lineWidth`
   property exercised in step 2 (it is a canvas-drawn slider control with no
   headless handle), so it is covered by the step-2 assertion, not driven separately.
6. Change **Labels Orientation** to Vert -- axis labels rotate
7. Change **Min Max Orientation** to Vert -- min/max values rotate
8. Adjust **Horz Margin** -- spacing between axes changes
9. Toggle **Auto Layout** off -- manual control over labels, min/max, and margins

**Actuation note:** line widths, orientation and margins are pure painting, so this
is a no-error floor; the axis sliders are read afterwards to confirm the layout
pass rebuilt the plot.

## Reset and filter visibility from the context menu

1. Narrow range sliders on two axes (AGE, HEIGHT) -- df.filter.trueCount drops progressively
2. Right-click > **Reset View** -- the in-chart range filter resets and df.filter.trueCount returns to the full count
3. Right-click > **Filter** > toggle **Show Filters** -- the showFilters state round-trips (true -> false -> true)

**Actuation note:** the axis-slider DOM elements persist regardless of Show Filters
(the range-handle visuals are canvas-drawn), so the assertable signal for Show
Filters is the `showFilters` prop the menu item flips, not a DOM count.

## Filter panel interaction

Precondition: the PC plot shows AGE, HEIGHT, WEIGHT (axes carry range sliders).

1. Open the Filter Panel and apply a histogram filter (AGE 30-50) -- df.filter.trueCount drops below the full count
2. Narrow the HEIGHT range slider on the PC plot -- df.filter.trueCount drops further (the Filter Panel and in-chart filters combine with AND)
3. Right-click > **Reset View** -- only the PC plot filtering resets; df.filter.trueCount returns to the Filter-Panel value (the panel filter survives)
4. Narrow the HEIGHT range slider again -- df.filter.trueCount drops below the Filter-Panel value once more
5. On the Filter Panel click **Reset filters** -- all filters reset and df.filter.trueCount returns to the full count (the PC plot range sliders return to default)

## Column reordering from the Context Panel list

1. Open Context Panel > **Value > Column Names** -- verify AGE, HEIGHT, WEIGHT are shown
2. Reorder columns in the list -- axes should reorder on the plot (the axis-slider name order matches the requested order)

## Density component toggles

1. Open Context Panel > **Box Plot > Density Style** -- default is **circles** (Show Density defaults on)
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

**Actuation note:** every density component is a canvas overlay, so the default
(densityStyle == 'circles') is asserted and the rest is a no-error floor.

## Color coding, legend & grid coloring

**Color column basics:**

1. Open Context Panel > **Color > Color Column** > set to AGE -- lines should be color-coded by age
2. Change **Color Axis Type** to **log** -- color mapping changes
3. Toggle **Invert Color Scheme** -- gradient reverses
4. Set custom **Color Min** and **Color Max** values -- color range narrows

**Legend position and visibility:**

5. Set Color Column to RACE -- the legend lists the RACE categories
6. Change **Legend Position** to Left, Right, Top, Bottom -- legend moves accordingly
7. Set **Legend Visibility** to Never -- legend disappears
8. Set **Legend Visibility** to Auto -- legend reappears with the same categories

**Color coding from grid:**

9. Open Context Panel > **Color > Color Column** > set to HEIGHT
10. In the grid, for the HEIGHT column set conditional color coding (use js-api) -- the PC plot legend surfaces the conditional bins ('20-150', '150-250')
11. In the grid, for the HEIGHT column set linear color coding (use js-api) -- the DOM legend disappears (a linear/numeric gradient has no DOM legend; it is canvas-drawn)

## Title and description

1. Open Context Panel > set **Title** to "My PC Plot" -- the title appears in the panel titlebar
2. Set **Description** to "Test description" -- description text appears inside the viewer
3. Change **Description Position** -- description moves to the specified position
4. Clear both fields -- title and description disappear

## Pick Up / Apply

1. Add a second PC plot
2. On the first PC plot, change the set of axes (e.g. remove HEIGHT)
3. Switch AGE axis to log scale
4. Set color column to RACE and move legend to the left
5. Set title to "Source Plot"
6. Right-click the first PC plot > **Pick up/Apply > Pick up**
7. Right-click the second PC plot > **Pick up/Apply > Apply** -- second plot should match the first (axes, log scale, color, legend position, title)
8. Change the axes on the first PC plot -- the second plot should not be affected
9. Adjust the range slider on the second PC plot -- it filters the shared DataFrame (df.filter.trueCount drops), so the first plot updates to show the filtered lines while the second plot keeps its own axes

## Layout round-trip

1. Save the layout
2. Open a scatterplot
3. Apply the saved layout -- verify only the original viewers are displayed (PC plot present, scatterplot absent)

## Project save / Close All / reopen

1. Save the project through the ribbon **Save** button (only the UI Save button captures the view layout)
2. Close All
3. Open the saved project -- verify the PC plot and its title are restored

## Table switching and transformation

1. Add a PC plot on the spgi-100 table
2. Go to Context Panel > Data and set **Table** to spgi-100 -- the axes are the spgi numeric columns
3. In the Data section click the Transformation input field and enter [{"#type":"GroupAggregation","aggType":"key","colName":"Chemist 521"},{"#type":"GroupAggregation","aggType":"pivot","colName":"Series"},{"#type":"GroupAggregation","aggType":"count","colName":"Id"}]
4. Verify the PC plot updates to show the pivoted aggregated data -- the axes become the Series categories (contain 'Triazoles')
5. Clear the Transformation -- the axes revert to the original spgi columns
6. Close All

**Note:** the pivot replaces the raw numeric axes with one generated column per
Series value, so the slider names are the signal that the aggregation was applied
and then reverted.

---
{
  "order": 7,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
