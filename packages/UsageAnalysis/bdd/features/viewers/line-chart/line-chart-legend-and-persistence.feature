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
