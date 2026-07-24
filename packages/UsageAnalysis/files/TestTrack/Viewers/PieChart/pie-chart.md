---
feature: piechart
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas:
  - piechart.in-viewer-column-selector-reconfigures-pie
  - piechart.mouseover-row-group-cross-highlight
realizes: []
realized_as:
  - pie-chart-spec.ts
related_bugs: []
expected_results:
  - anchor: "In-viewer column selector re-pick"
    expectation: "Picking SEX through the selector shown on the pie itself
      rebinds the category column to SEX and the legend lists exactly the SEX
      categories; picking RACE back restores the RACE split with exactly the
      RACE categories in the legend (round-trip through the real selector UI,
      no substitution)"
  - anchor: "Mouse-over row group cross-highlight"
    expectation: "With Show Mouse Over Row Group off, highlighting one
      category's rows leaves the pie canvas unchanged (per-color histogram
      delta below 2000 between settled frames); with the option on, the same
      highlight repaints the matching arc (delta above 20000); clearing the
      highlight returns the settled frame to the baseline (delta below 2000,
      round-trip)"
  - anchor: "Sorting"
    expectation: "With Category set to RACE, Pie Sort Type and Pie Sort Order
      read back each value in the driven sequence: by value, desc, asc,
      by category, asc, desc"
  - anchor: "Appearance"
    expectation: "Start Angle 90/180/0, Max Radius 100/150 and Shift 10/0 each
      read back the value just set"
  - anchor: "Labels"
    expectation: "Label Position Inside/Outside/Auto and the Show Label, Show
      Percentage and Show Value toggles read back the values just set
      (Inside, Outside, Auto, false, false, true)"
  - anchor: "Outline"
    expectation: "Outline Line Width 5/0/1 reads back each value in sequence"
  - anchor: "Column selector"
    expectation: "Show Column Selector reads back false then true across the
      off/on toggle"
  - anchor: "Legend"
    expectation: "Legend Visibility Always/Never/Auto and Legend Position
      LeftTop/RightBottom read back each value in the driven sequence"
  - anchor: "Row source"
    expectation: "With the first 50 rows selected, Row Source reads back
      Selected, Filtered and All in sequence"
  - anchor: "Title and description"
    expectation: "showTitle true, title 'Demographics', description 'By race',
      descriptionPosition Top and descriptionVisibilityMode Never each read
      back the value just set"
  - anchor: "Layout persistence"
    expectation: "After saving the layout, closing the pie chart and
      re-applying the saved layout, the restored viewer carries the same
      configuration — categoryColumnName RACE, segmentAngleColumnName AGE,
      startAngle 45, shift 5; the saved layout is deleted afterwards"
  - anchor: "Selection and interaction"
    expectation: "Show Selected Rows and Show Mouse Over Row Group each read
      back false then true across their off/on toggles"
  - anchor: "Auto layout"
    expectation: "Auto Layout off, Margin Left 50, Margin Top 50 and Auto
      Layout back on each read back the value just set"
  - anchor: "Table switching and row source (SPGI)"
    expectation: "Setting the Table property to SPGI rebinds the viewer's
      dataframe to SPGI and back to demog on the reverse switch; Row Source
      reads back Selected with 100 rows selected, and Filtered with
      df.filter.trueCount above zero after a RACE = Asian categorical filter"
---

# Pie chart tests

## Purpose

Verifies the Pie chart's settings surface — sorting, appearance, labels,
outline, legend, row source, title and description, auto layout — plus layout
persistence, re-picking the category column through the selector shown on the
pie itself, the cross-highlight of a hovered row group, and switching the
bound table. Most settings render to canvas only, so they are verified by
reading each setting back after the change; the two focused scenarios add
real-control and repaint checks on top.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add a Pie chart by clicking the Pie Chart icon in the Toolbox > Viewers section

## Sorting

1. Set Category to RACE
2. Go to the Context Panel > Data, set Pie Sort Type to by value — slices ordered by size
3. Set Pie Sort Order to desc — largest slice first
4. Set Pie Sort Order to asc — smallest slice first
5. Set Pie Sort Type to by category — slices ordered alphabetically
6. Set Pie Sort Order to asc — A-Z order
7. Set Pie Sort Order to desc — Z-A order

## Appearance

1. Set Start Angle to 90 — pie rotates 90 degrees
2. Set Start Angle to 180 — further rotation
3. Set Start Angle to 0
4. Set Max Radius to 100 — pie size decreases
5. Set Max Radius to 150 — pie size restored
6. Set Shift to 10 — slices separate from center (exploded view)
7. Set Shift to 0

## Labels

1. Go to the Context Panel > Style, set Label Position to Inside — labels inside slices
2. Set Label Position to Outside — labels outside with lines
3. Set Label Position to Auto — labels placed automatically
4. Toggle Show Label off — category names disappear
5. Toggle Show Percentage off — percentage values disappear
6. Toggle Show Value on — absolute count values appear
7. Re-enable Show Label and Show Percentage

## Outline

1. Set Outline Line Width to 5 — thick outlines around slices
2. Set Outline Line Width to 0 — no outlines
3. Set Outline Line Width to 1

## Column selector

1. Toggle Show Column Selector off — category dropdown disappears
2. Toggle Show Column Selector on — dropdown reappears

## In-viewer column selector re-pick

Steps:
1. Set Category to RACE, turn Show Column Selector on, and make the legend always visible.
2. Using the column selector shown on the pie itself, pick the SEX column.
3. Using the same selector, pick RACE again.
4. Restore Legend Visibility to Auto.

Expected:
- Picking SEX re-splits the pie by SEX: the category is now SEX and the legend lists exactly the SEX categories
- Picking RACE back restores the RACE split: the category is RACE and the legend lists exactly the RACE categories (round-trip)

## Legend

1. Go to the Context Panel > Legend, set Legend Visibility to Always — legend appears
2. Set Legend Position to Left Top — legend moves to top-left
3. Set Legend Position to Right Bottom — legend moves to bottom-right
4. Set Legend Visibility to Never — legend disappears
5. Set Legend Visibility to Auto — legend shown based on available space

## Row source

1. Select first 50 rows in the grid
2. Set Row Source to Selected — pie shows only selected rows
3. Set Row Source to Filtered — pie shows filtered rows
4. Set Row Source to All — pie shows all rows

## Title and description

1. Set Show Title to true
2. Set Title to "Demographics" — title appears on the viewer
3. Set Description to "By race" — description appears
4. Set Description Position to Top — description moves to top
5. Set Description Visibility Mode to Never — description disappears
6. Clear the title

## Layout persistence

1. Set Category to RACE
2. Set Segment Angle Column to AGE
3. Set Start Angle to 45
4. Set Shift to 5
5. Save the layout
6. Close the pie chart viewer by clicking the X icon on the viewer title bar
7. Apply the saved layout
8. Verify Category is RACE, Segment Angle Column is AGE, Start Angle is 45, Shift is 5
9. Delete the saved layout

## Selection and interaction

1. Set Category to RACE
2. Toggle Show Selected Rows off
3. Toggle Show Selected Rows on
4. Toggle Show Mouse Over Row Group off
5. Toggle Show Mouse Over Row Group on

## Mouse-over row group cross-highlight

Steps:
1. Set Category to RACE.
2. With Show Mouse Over Row Group off, highlight the rows of one category (for example Asian) — the pie does not change.
3. Turn Show Mouse Over Row Group on and highlight the same rows — the matching arc lights up with the highlight color.
4. Clear the highlight — the pie returns to its previous look.

Expected:
- With the option off, highlighting a category's rows repaints nothing on the pie
- With the option on, the same highlight visibly repaints the pie (the matching arc lights up)
- Clearing the highlight returns the pie to its previous look (round-trip)

## Auto layout

1. Toggle Auto Layout off
2. Set Margin Left to 50 — left margin increases
3. Set Margin Top to 50 — top margin increases
4. Toggle Auto Layout on — margins return to automatic sizing
5. Resize the viewer to be very small — take a screenshot to verify labels auto-hide
6. Resize the viewer back to normal size

## Table switching and row source (uses SPGI dataset)

Setup: Close all, open demog, then also open SPGI

1. Add a Pie chart on demog by clicking the Pie Chart icon in the Toolbox > Viewers section
2. Go to the Context Panel > Data, switch the table to SPGI — the pie chart should update to SPGI columns
3. Switch back to demog — the pie chart should restore demog columns
4. Set Row Source to Selected
5. Select some rows in the grid — the pie chart should show only selected rows
6. Set Row Source to Filtered
7. Open the filter panel and apply a filter — the pie chart should reflect the filter

## Automation notes

The settings sections (Sorting through Auto layout, Layout persistence,
Selection and interaction, and the Table/Row Source switches) are driven
through the matching viewer properties — the same values the Context Panel
controls write — and each value is read back after the set. The two focused
scenarios below are driven through the real UI as described.

In-viewer column selector re-pick: the pie's own column combo is a DOM control
(`[name="div-column-combobox-category"]`), driven through the shared
type-and-search selector helper with no JS-API substitution. The read-back is
the category property plus the legend labels compared to the column's own
category list, so a pick that did not re-split the pie fails instead of echoing
the prop.

Mouse-over row group cross-highlight: the pie highlights the mouse-over ROW
GROUP, not the single mouse-over row — hovering one grid row leaves the pie
unchanged — so "highlight the rows" in the manual steps is driven through the
dataframe row-highlight channel (`df.rows.highlight`), the same channel other
viewers write when an aggregate element is hovered. The repaint is measured as
a per-color canvas histogram delta between settled frames; the thresholds
(off below 2000, on above 20000, cleared below 2000) are conservative margins
far below the delta a real arc highlight produces, and every measured delta is
logged for threshold audits. The reverse direction — hovering a pie slice
highlights the contributing rows in the grid — has no readable headless signal
(the grid is canvas-rendered) and stays in pie-chart-ui.md.

---
{
  "order": 4,
  "datasets": ["System:DemoFiles/demog.csv", "System:DemoFiles/chem/SPGI.csv"]
}
