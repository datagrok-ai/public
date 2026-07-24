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
      Show All Lines items flip the matching prop and restore it on a second click.
      With Show All Lines off and an empty selection the canvas ink drops far below
      the all-lines count; selecting grid rows paints only the selected lines (ink
      rises off the hidden floor yet stays below the all-lines total); re-enabling
      Show All Lines brings the ink back up (settle-gated pixel-count deltas)."
  - anchor: "Style & layout"
    expectation: "Driving line widths, label/min-max orientation, margins and Auto
      Layout raises no error; the axis sliders are still present after the layout
      pass. The Lines > Line Width menu slider writes the same lineWidth prop
      (documented reduction — the menu slider is a canvas control)."
  - anchor: "Show Filters from the context menu"
    expectation: "Show Filters in the context-menu Filter submenu round-trips the
      showFilters prop (true -> false -> true)."
  - anchor: "Title and description"
    expectation: "The title renders in the panel titlebar (.panel-titlebar-text) and
      the description inside the viewer element; clearing both removes them."
  - anchor: "Pick Up / Apply"
    expectation: "Pick Up on the first plot then Apply on the second copies axes,
      log-scale columns, color column, legend position and title. Changing the
      first plot's axes afterwards does not affect the second; a range slider on
      the second plot filters the shared DataFrame (df.filter.trueCount drops)."
  - anchor: "Table switching and transformation"
    expectation: "The PC plot bound to spgi renders the raw numeric axes; applying a
      GroupAggregation pivot on Series replaces them with the Series categories
      (axes contain 'Triazoles'); clearing the Transformation reverts the axes to
      the original set."
---

# PC plot tests

## Purpose

Verifies the PC Plot's core surface: adding the viewer and generating a script
from it, axis normalization from the context menu, selection and line-display
toggles, style and layout settings, filter visibility, title and description,
copying settings between two plots (Pick Up / Apply), and binding the plot to
another table with an aggregation transformation. Painting details that have no
readable outcome (which exact lines are drawn, how thick they look) are checked
only as far as the plot keeps working and its ink visibly changes; the
picture's correctness is not judged.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add PC plot

## Menu Ribbon and To Script

1. On the Menu Ribbon, click the **Add viewer** icon and select **PC plot** -- viewer should open
2. Right-click the plot > **To Script** > **To JavaScript** -- a balloon with the generated script should appear
3. Close the viewer, then add PC plot via **Toolbox** -- viewer opens again

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
9. Right-click > **Selection** > toggle **Show Current Line** -- same effect as the Context Panel toggle
10. Right-click > **Selection** > toggle **Show All Lines** -- same effect as the Context Panel toggle
11. Disable **Show All Lines** with an empty selection -- almost all ink disappears from the plot (only axes/labels remain)
12. Select some rows in the grid -- only the selected lines are painted (ink rises off the hidden floor but stays well below the all-lines total)
13. Re-enable **Show All Lines** -- all lines are painted again (ink returns toward the all-lines total)
14. Clear the selection -- round-trip back to the default state

## Style & layout

1. Open Context Panel > **Style**
2. Change **Line Width** to 3 -- lines become thicker
3. Change **Current Line Width** to 5 -- current row line becomes prominent
4. Change **Mouse Over Line Width** to 5 -- hover line becomes thicker
5. Right-click > **Lines** > adjust **Line Width** slider -- line thickness changes
6. Change **Labels Orientation** to Vert -- axis labels rotate
7. Change **Min Max Orientation** to Vert -- min/max values rotate
8. Adjust **Horz Margin** -- spacing between axes changes
9. Toggle **Auto Layout** off -- manual control over labels, min/max, and margins

## Show Filters from the context menu

1. Right-click > **Filter** > toggle **Show Filters** twice -- the filter visibility setting flips and then returns to its original state

## Title and description

1. Open Context Panel > set **Title** to "My PC Plot" -- the title appears in the panel titlebar
2. Set **Description** to "Test description" -- description text appears inside the viewer
3. Change **Description Position** -- description moves to the specified position
4. Clear both fields -- title and description disappear

## Pick Up / Apply

1. Add a second PC plot
2. On the first PC plot, change the set of axes in Context Panel > Value > **Column Names** (e.g. remove HEIGHT)
3. Switch the AGE axis to log scale: in Context Panel > Value, mark AGE in the **Log Columns** selector (grouped with Column Names)
4. Set **Color** (Context Panel > Color) to RACE and set **Legend Position** to **Left**
5. Set title to "Source Plot"
6. Right-click the first PC plot > **Pick up/Apply > Pick up**
7. Right-click the second PC plot > **Pick up/Apply > Apply** -- second plot should match the first (axes, log scale, color, legend position, title)
8. Change the axes on the first PC plot -- the second plot should not be affected
9. Adjust the range slider on the second PC plot -- it filters the shared table (the filtered row count drops), so the first plot updates to show the filtered lines while the second plot keeps its own axes

## Table switching and transformation

1. Add a PC plot on the spgi-100 table
2. Go to Context Panel > Data and set **Table** to spgi-100 -- the axes are the spgi numeric columns
3. In the Data section click the Transformation input field and enter [{"#type":"GroupAggregation","aggType":"key","colName":"Chemist 521"},{"#type":"GroupAggregation","aggType":"pivot","colName":"Series"},{"#type":"GroupAggregation","aggType":"count","colName":"Id"}]
4. Verify the PC plot updates to show the pivoted aggregated data -- the axes become the Series categories (contain 'Triazoles')
5. Clear the Transformation -- the axes revert to the original spgi columns
6. Close All

## Automation notes

Menu Ribbon and To Script: the spec adds the viewer through the Toolbox icon
(`[name="icon-pc-plot"]`) for both adds. The Menu Ribbon **Add viewer** gallery is
a canvas-rendered dialog with no headless handles, so it is not exercised; the
Toolbox path produces the same viewer.

Style & layout, step 5: the Lines menu **Line Width** slider writes the same
`lineWidth` property exercised in step 2 (it is a canvas-drawn slider control with
no headless handle), so automation covers it via the step-2 assertion and does not
drive the slider separately; the manual step remains valid by hand.

Selection & line display: which lines get painted is a canvas outcome with no DOM
counterpart, so steps 2-8 are a no-error floor over the prop surface. Steps 9-10
assert the menu -> prop round-trip, the same state the Context Panel toggles write.
Steps 11-14 DO carry a canvas signal — fewer painted lines means less ink — so
they are asserted as settle-gated pixel-count deltas (baseline precheck with a
>= 0 guard and a ceiling against the all-lines count, then deltas in both
directions: selection-only paints fewer lines, re-enable brings them all back).

Style & layout: line widths, orientation and margins are pure painting, so this
is a no-error floor; the axis sliders are read afterwards to confirm the layout
pass rebuilt the plot.

Show Filters from the context menu: the axis-slider DOM elements persist regardless of Show Filters
(the range-handle visuals are canvas-drawn), so the assertable signal for Show
Filters is the `showFilters` prop the menu item flips, not a DOM count. The
in-chart range-filter narrowing and Reset View restore are owned by
pcplot-setup-color-filter.md (its in-chart range-filter drop + Reset View
restore scenario).

Pick Up / Apply: "the filtered row count" is read as `df.filter.trueCount` on
the shared dataframe.

Table switching and transformation: the pivot replaces the raw numeric axes
with one generated column per Series value, so the axis-slider names are the
signal that the aggregation was applied and then reverted.

---
{
  "order": 7,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
