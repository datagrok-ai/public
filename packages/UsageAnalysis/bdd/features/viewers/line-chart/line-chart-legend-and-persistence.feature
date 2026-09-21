@journey @viewers @realizes:viewers.line-chart
Feature: Line chart legend, category colours and what survives a round-trip
  The legend that a split brings with it — where it sits, when it is hidden, what a click on one of
  its categories does to the chart — and whether the colours it shows come back from a saved layout
  and from a saved project.
  The old client-side spec read the legend by walking `.d4-legend` for an element whose computed
  display was not "none" and splitting its `innerText`; the legend publishes `data-legend-mode`,
  `data-legend-slot`, `data-legend-items` and a key per item, which is what the library's legend
  steps read. The click-filters-the-chart claim was `lc.filter.trueCount < rowCount`, an inequality
  that holds for any of the five categories: `S_ABS` has exactly 2 rows on spgi-100, and the chart
  then draws one line instead of five.
  The layout and project round-trips are the server lane of the old section (GROK-17278,
  GROK-19825), folded in here because they are the same subject and the library has both steps.
  The last scenarios translate the Legend section's line chart cases on demog-1000, opened next to
  spgi-100 (a line chart of WEIGHT by STARTED split by DIS_POP, whose six categories replace Series;
  WEIGHT, HEIGHT and AGE replace Average Mass, TPSA and NIBR logP): the default palette gives
  neighbouring categories different colors (a hue apart — two shades of one hue count as alike);
  two Y columns share one six-item legend until Multi Axis gives every Y column its own block of
  six items ("WEIGHT / RA", "HEIGHT / RA"), every item named; replacing a Y column replaces its
  block and leaves the other one with its colors; and Split, Multi Axis, the blocks and the
  replaced column come back from a layout and from a project. The Y column is replaced through
  the property, as the manual case allows.
  Multi Axis does not part company with the manual case; what it draws depends on how many Y columns
  there are, as the operator confirmed. With two Y columns the chart folds them into one box with a
  scale for each (line-chart-multi-axis-and-split.feature asserts the boxes and the scales); with
  three or more it draws one box and one scale for all of them. Both cases are claimed here, and the
  contrast is the claim: with two columns the chart reports two `y axes` and a `y2 axis` area, with
  three it reports no `y2 axis` and its `y axes` reading falls to 0, the columns sharing the one
  scale that is left. Either way every Y column keeps a legend block of its own, which is where the
  per-column split is read, together with the lines.

  Background:
    Given user is logged in
    And user opens spgi dataset
    And user colors "Stereo Category" column categorically:
      | R_ONE  | #FF0000 |
      | S_ABS  | #00FF00 |
      | S_ACHIR| #0000FF |
      | S_PART | #FFFF00 |
      | S_UNKN | #FF00FF |
    And user adds a line chart viewer with:
      | xColumnName      | Chemical Space X |
      | yColumnNames     | Chemical Space Y |
      | splitColumnNames | Stereo Category  |
    Then 100 rows should pass the filter
    And the "lines" reading of line chart viewer should be 5
    And the legend of line chart viewer should list 5 items
    And the categorical color of "R_ONE" in "Stereo Category" column should be "#FF0000"
    And line chart viewer should report no error

  Scenario: The legend comes with the split and goes when the split does
    Then the legend of line chart viewer should list 5 items
    When user sets "splitColumnNames" property of line chart viewer to ""
    Then the legend of line chart viewer should be hidden
    And the "lines" reading of line chart viewer should be 1
    When user sets "splitColumnNames" property of line chart viewer to "Stereo Category"
    Then the legend of line chart viewer should list 5 items
    And no errors should have been logged

  Scenario: Legend Position moves the legend to the side it names
    When user sets "legendPosition" property of line chart viewer to "Left"
    Then the legend of line chart viewer should be on the left
    When user sets "legendPosition" property of line chart viewer to "Top"
    Then the legend of line chart viewer should be on the top
    When user sets "legendPosition" property of line chart viewer to "Bottom"
    Then the legend of line chart viewer should be on the bottom
    When user sets "legendPosition" property of line chart viewer to "Right"
    Then the legend of line chart viewer should be on the right
    And the legend of line chart viewer should list 5 items
    When user sets "legendPosition" property of line chart viewer to "Auto"
    Then no errors should have been logged

  Scenario: Legend Visibility Never hides the legend and Auto brings it back
    When user sets "legendVisibility" property of line chart viewer to "Never"
    Then the legend of line chart viewer should be hidden
    And the "lines" reading of line chart viewer should be 5
    When user sets "legendVisibility" property of line chart viewer to "Always"
    Then the legend of line chart viewer should list 5 items
    When user sets "legendVisibility" property of line chart viewer to "Auto"
    Then the legend of line chart viewer should list 5 items
    And no errors should have been logged

  Scenario: Clicking a legend category filters the chart down to it and back
    When user clicks on "S_ABS" item in the legend of line chart viewer
    Then the "rows shown" reading of line chart viewer should be 2
    And the "lines" reading of line chart viewer should be 1
    And 100 rows should pass the filter
    And the categorical color of "R_ONE" in "Stereo Category" column should be "#FF0000"
    And line chart viewer should have repainted
    When user clicks on "S_ABS" item in the legend of line chart viewer
    Then the "rows shown" reading of line chart viewer should be 100
    And the "lines" reading of line chart viewer should be 5
    And no errors should have been logged

  Scenario: The legend items keep the colours the column was given
    Then the "R_ONE" item in the legend of line chart viewer should be colored "#FF0000"
    And the "S_ABS" item in the legend of line chart viewer should be colored "#00FF00"
    And the "R_ONE" and "S_ABS" items in the legend of line chart viewer should be colored differently
    And no errors should have been logged

  Scenario: The chart configuration and the category colour come back from a saved layout (GROK-17278)
    When user saves the layout of the current table view to the server
    And user colors "Stereo Category" column categorically:
      | R_ONE | #00AAFF |
    Then the categorical color of "R_ONE" in "Stereo Category" column should be "#00AAFF"
    When user sets "splitColumnNames" property of line chart viewer to ""
    Then the "lines" reading of line chart viewer should be 1
    When user loads the saved layout
    Then the "lines" reading of line chart viewer should be 5
    And the legend of line chart viewer should list 5 items
    And the categorical color of "R_ONE" in "Stereo Category" column should be "#FF0000"
    And the "x column" reading of line chart viewer should be "Chemical Space X"
    And no errors should have been logged

  Scenario: The chart and its colours come back from a saved project (GROK-19825)
    When user saves the current view as project "zz-linechart-legend-colors"
    And user closes all views
    And user opens the "zz-linechart-legend-colors" project
    Then line chart viewer should be added to the open tableview
    And the "lines" reading of line chart viewer should be 5
    And the legend of line chart viewer should list 5 items
    And the categorical color of "R_ONE" in "Stereo Category" column should be "#FF0000"
    And line chart viewer should be painted
    And no errors should have been logged

  Scenario: The default palette gives neighbouring split categories different colors
    Given user opens demog-1000 dataset
    And user adds a line chart viewer with:
      | xColumnName       | STARTED |
      | yColumnNames      | WEIGHT  |
      | splitColumnNames  | DIS_POP |
      | Legend Visibility | Always  |
      | Legend Position   | Right   |
    Then the legend of line chart viewer should list 6 items
    And the "AS" and "Indigestion" items in the legend of line chart viewer should be colored differently
    And the "Indigestion" and "PsA" items in the legend of line chart viewer should be colored differently
    And the "PsA" and "Psoriasis" items in the legend of line chart viewer should be colored differently
    And the "Psoriasis" and "RA" items in the legend of line chart viewer should be colored differently
    And the "RA" and "UC" items in the legend of line chart viewer should be colored differently
    And no errors should have been logged

  Scenario: Multi Axis gives every Y column its own block of legend items
    When user sets "yColumnNames" property of line chart viewer to "WEIGHT, HEIGHT"
    Then the legend of line chart viewer should list 6 items
    When user sets "multiAxis" property of line chart viewer to "true"
    Then the legend of line chart viewer should list 12 items
    And the "lines" reading of line chart viewer should be 11
    And the following elements should be visible:
      | "WEIGHT / AS" legend item in legend of line chart viewer |
      | "WEIGHT / Indigestion" legend item in legend of line chart viewer |
      | "WEIGHT / PsA" legend item in legend of line chart viewer |
      | "WEIGHT / Psoriasis" legend item in legend of line chart viewer |
      | "WEIGHT / RA" legend item in legend of line chart viewer |
      | "WEIGHT / UC" legend item in legend of line chart viewer |
      | "HEIGHT / AS" legend item in legend of line chart viewer |
      | "HEIGHT / Indigestion" legend item in legend of line chart viewer |
      | "HEIGHT / PsA" legend item in legend of line chart viewer |
      | "HEIGHT / Psoriasis" legend item in legend of line chart viewer |
      | "HEIGHT / RA" legend item in legend of line chart viewer |
      | "HEIGHT / UC" legend item in legend of line chart viewer |
    And the "WEIGHT / RA" and "HEIGHT / RA" items in the legend of line chart viewer should be colored differently
    And the "y axes" reading of line chart viewer should be 2
    And line chart viewer should have a "y2 axis" area
    And no errors should have been logged

  Scenario: Three Y columns under Multi Axis share one scale
    When user sets "yColumnNames" property of line chart viewer to "AGE, HEIGHT, WEIGHT"
    Then the "y columns" reading of line chart viewer should be "AGE, HEIGHT, WEIGHT"
    And the "charts" reading of line chart viewer should be 1
    And the "y axes" reading of line chart viewer should be 0
    And line chart viewer should not have a "y2 axis" area
    And the legend of line chart viewer should list 18 items
    And "AGE / RA" legend item in legend of line chart viewer should be visible
    And "HEIGHT / RA" legend item in legend of line chart viewer should be visible
    And "WEIGHT / RA" legend item in legend of line chart viewer should be visible
    When user sets "yColumnNames" property of line chart viewer to "WEIGHT, HEIGHT"
    Then the "y axes" reading of line chart viewer should be 2
    And line chart viewer should have a "y2 axis" area
    And the legend of line chart viewer should list 12 items
    And no errors should have been logged

  Scenario: Replacing a Y column replaces its block and leaves the other one
    Then the "WEIGHT / RA" item in the legend of line chart viewer should be colored "#9467BD"
    And the "WEIGHT / AS" item in the legend of line chart viewer should be colored "#1F77B4"
    When user sets "yColumnNames" property of line chart viewer to "WEIGHT, AGE"
    Then the "y columns" reading of line chart viewer should be "WEIGHT, AGE"
    And the legend of line chart viewer should list 12 items
    And the "lines" reading of line chart viewer should be 12
    And the following elements should be visible:
      | "WEIGHT / AS" legend item in legend of line chart viewer |
      | "WEIGHT / Indigestion" legend item in legend of line chart viewer |
      | "WEIGHT / PsA" legend item in legend of line chart viewer |
      | "WEIGHT / Psoriasis" legend item in legend of line chart viewer |
      | "WEIGHT / RA" legend item in legend of line chart viewer |
      | "WEIGHT / UC" legend item in legend of line chart viewer |
      | "AGE / AS" legend item in legend of line chart viewer |
      | "AGE / Indigestion" legend item in legend of line chart viewer |
      | "AGE / PsA" legend item in legend of line chart viewer |
      | "AGE / Psoriasis" legend item in legend of line chart viewer |
      | "AGE / RA" legend item in legend of line chart viewer |
      | "AGE / UC" legend item in legend of line chart viewer |
    And "HEIGHT / RA" legend item in legend of line chart viewer should be absent
    And the "WEIGHT / RA" item in the legend of line chart viewer should be colored "#9467BD"
    And the "WEIGHT / AS" item in the legend of line chart viewer should be colored "#1F77B4"
    And no errors should have been logged

  Scenario: Split, Multi Axis and the Y blocks come back from a saved layout
    When user saves the layout of the current table view to the server
    And user sets properties of line chart viewer:
      | multiAxis        | false  |
      | splitColumnNames |        |
      | yColumnNames     | WEIGHT |
    Then "multiAxis" property of line chart viewer should be "false"
    And the "y columns" reading of line chart viewer should be "WEIGHT"
    And the legend of line chart viewer should be hidden
    When user loads the saved layout
    Then "multiAxis" property of line chart viewer should be "true"
    And the "y columns" reading of line chart viewer should be "WEIGHT, AGE"
    And the "lines" reading of line chart viewer should be 12
    And the legend of line chart viewer should list 12 items
    And the following elements should be visible:
      | "WEIGHT / AS" legend item in legend of line chart viewer |
      | "WEIGHT / Indigestion" legend item in legend of line chart viewer |
      | "WEIGHT / PsA" legend item in legend of line chart viewer |
      | "WEIGHT / Psoriasis" legend item in legend of line chart viewer |
      | "WEIGHT / RA" legend item in legend of line chart viewer |
      | "WEIGHT / UC" legend item in legend of line chart viewer |
      | "AGE / AS" legend item in legend of line chart viewer |
      | "AGE / Indigestion" legend item in legend of line chart viewer |
      | "AGE / PsA" legend item in legend of line chart viewer |
      | "AGE / Psoriasis" legend item in legend of line chart viewer |
      | "AGE / RA" legend item in legend of line chart viewer |
      | "AGE / UC" legend item in legend of line chart viewer |
    And no errors should have been logged

  Scenario: Split, Multi Axis and the Y blocks come back from a saved project
    When user saves the current view as project "bdd-line-chart-legend-blocks"
    And user closes all views
    And user opens the "bdd-line-chart-legend-blocks" project
    Then line chart viewer should be added to the open tableview
    And "multiAxis" property of line chart viewer should be "true"
    And the "y columns" reading of line chart viewer should be "WEIGHT, AGE"
    And the "lines" reading of line chart viewer should be 12
    And the legend of line chart viewer should list 12 items
    And the following elements should be visible:
      | "WEIGHT / AS" legend item in legend of line chart viewer |
      | "WEIGHT / Indigestion" legend item in legend of line chart viewer |
      | "WEIGHT / PsA" legend item in legend of line chart viewer |
      | "WEIGHT / Psoriasis" legend item in legend of line chart viewer |
      | "WEIGHT / RA" legend item in legend of line chart viewer |
      | "WEIGHT / UC" legend item in legend of line chart viewer |
      | "AGE / AS" legend item in legend of line chart viewer |
      | "AGE / Indigestion" legend item in legend of line chart viewer |
      | "AGE / PsA" legend item in legend of line chart viewer |
      | "AGE / Psoriasis" legend item in legend of line chart viewer |
      | "AGE / RA" legend item in legend of line chart viewer |
      | "AGE / UC" legend item in legend of line chart viewer |
    And "HEIGHT / RA" legend item in legend of line chart viewer should be absent
    And no errors should have been logged
