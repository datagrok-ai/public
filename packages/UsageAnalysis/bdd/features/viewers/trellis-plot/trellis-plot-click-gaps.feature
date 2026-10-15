@journey @viewers @realizes:viewers.trellis-plot
Feature: Trellis plot On Click set in the context panel, and what a change keeps
  Row Source and On Click set in the context panel, the way a person sets them, correct each other as
  they do through the API (GROK-17711): On Click = Filter moves Row Source to All, Row Source =
  Filtered moves On Click back to None. Under On Click = Select a change of split column keeps the
  selection a click made; under On Click = Filter a change of the inner viewer keeps the filter a
  click made. With a filter card on, Escape takes back only the trellis's part of the filter.

  The inner viewer is changed in the trellis's own viewer selector, as a person does (MISSING.md,
  Resolved).
  Translated from TestTrack Viewers/TrellisPlot/trellis-plot.md "On Click functionality" steps 5, 8
  and 12 (the Select and Filter scenarios), trellis-plot-click-to-filter.md section 1 step 10 (the
  filter card) and section 2 steps 1-4 (the context panel). One journey on demog-1000 with SEX by RACE
  and a scatter plot inside; each scenario sets the On Click it needs.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names | SEX          |
      | Y Column Names | RACE         |
      | Viewer Type    | Scatter plot |
    Then the "cells" reading of trellis plot viewer should be 8
    And "On Click" property of trellis plot viewer should be "None"

  Scenario: Under Select, a change of split column keeps the selection
    When user sets "On Click" property of trellis plot viewer to "Select"
    And user clicks on the "cell M | Asian" area of trellis plot viewer
    Then 8 rows should be selected
    And no rows where "SEX" is "F" should be selected
    And no rows where "RACE" is "Caucasian" should be selected
    When user sets "X Column Names" property of trellis plot viewer to "CONTROL"
    Then trellis plot viewer should not have a "cell M | Asian" area
    And trellis plot viewer should have a "cell false | Asian" area
    And 8 rows should be selected
    When user sets "X Column Names" property of trellis plot viewer to "SEX"
    And user clears the row selection
    And user sets "On Click" property of trellis plot viewer to "None"
    Then no errors should have been logged

  Scenario: Under Filter, a change of the inner viewer keeps the filter, and Escape then drops it
    When user sets "On Click" property of trellis plot viewer to "Filter"
    And user clicks on the "cell F | Caucasian" area of trellis plot viewer
    Then 480 rows should pass the filter
    When user picks "Bar chart" in the viewer selector of trellis plot viewer
    Then the "inner viewer type" reading of trellis plot viewer should be "Bar chart"
    And the "current cell" reading of trellis plot viewer should be "F | Caucasian"
    And 480 rows should pass the filter
    When user presses Escape in trellis plot viewer
    Then all rows should pass the filter
    When user picks "Scatter plot" in the viewer selector of trellis plot viewer
    Then the "inner viewer type" reading of trellis plot viewer should be "Scatter plot"
    And no errors should have been logged

  Scenario: With a filter card on, Escape takes back only the trellis's part
    When user adds a categorical filter on "DIS_POP" keeping "RA"
    Then 434 rows should pass the filter
    When user clicks on the "cell M | Asian" area of trellis plot viewer
    Then 1 row should pass the filter
    And no rows where "SEX" is "F" should pass the filter
    When user presses Escape in trellis plot viewer
    Then 434 rows should pass the filter
    And all rows where "DIS_POP" is "RA" should pass the filter
    When user hovers over "DIS_POP" filter card
    And user clicks on close of "DIS_POP" filter card
    Then all rows should pass the filter
    When user sets "On Click" property of trellis plot viewer to "None"
    Then no errors should have been logged

  Scenario: Set in the context panel, On Click and Row Source correct each other
    When user clicks on settings icon of trellis plot viewer
    Then context panel should be visible
    When user selects "Filtered" in "Row Source" property in context panel
    And "Misc" category in context panel is expanded
    When user selects "Filter" in "On Click" property in context panel
    Then "On Click" property of trellis plot viewer should be "Filter"
    And "Row Source" property of trellis plot viewer should be "All"
    When user clicks on the "cell F | Caucasian" area of trellis plot viewer
    Then 480 rows should pass the filter
    When user presses Escape in trellis plot viewer
    Then all rows should pass the filter
    When user selects "Filtered" in "Row Source" property in context panel
    Then "Row Source" property of trellis plot viewer should be "Filtered"
    And "On Click" property of trellis plot viewer should be "None"
    When user sets "Row Source" property of trellis plot viewer to "All"
    Then no errors should have been logged
