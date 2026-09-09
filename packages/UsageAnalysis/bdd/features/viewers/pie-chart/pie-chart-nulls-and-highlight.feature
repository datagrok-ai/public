@journey @viewers @realizes:viewers.pie-chart
Feature: Pie chart missing values and the mouse-over row group
  A category column with gaps gets a wedge of its own, named "(no value)": it is clickable like any
  other, Include Nulls takes it away and puts it back, and the legend counts the empty category
  with the rest. And the mouse-over row group, which paints its share of the wedge it belongs to
  only while Show Mouse Over Row Group is on. One journey on demog-1000 with a calculated
  RACE_GAPS — RACE with the rows whose AGE is a multiple of ten blanked, so 112 blanks and
  Caucasian 793, Other 56, Black 25, Asian 14.
  Not translated: the reverse cross-highlight (hovering a grid row and watching the pie) from
  pie-chart-ui.md — the grid is canvas-rendered and there is no signal for the direction it
  travels; and the pixel-histogram deltas of the old spec's "Mouse-over row group cross-highlight",
  which the `mouse over <category>` overlay region replaces.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a calculated column "RACE_GAPS" with formula "if(Mod(${AGE}, 10) == 0, null, ${RACE})"
    And user adds a pie chart viewer with:
      | Category          | RACE_GAPS |
      | Legend Visibility | Always    |
    Then the "slices" reading of pie chart viewer should be 5
    And pie chart viewer should show 1000 rows

  Scenario: The blank rows get a wedge of their own
    Then pie chart viewer should have a "slice (no value)" area
    And the "angle value of (no value)" reading of pie chart viewer should be 112
    And the "share of (no value)" reading of pie chart viewer should be 11.2
    And the "angle value of Caucasian" reading of pie chart viewer should be 793
    And pie chart viewer should be painted
    And no errors should have been logged

  Scenario: Clicking the blank wedge selects exactly the blank rows
    When user clicks on the "slice (no value)" area of pie chart viewer
    Then 112 rows should be selected
    And no rows where "RACE_GAPS" is "Caucasian" should be selected
    And no rows where "RACE_GAPS" is "Asian" should be selected
    And pie chart viewer should have a "selected (no value)" area
    And pie chart viewer should show more selection highlight than before
    When user presses Escape
    Then no rows should be selected
    And no errors should have been logged

  Scenario: Include Nulls off drops the wedge and the legend item with it
    Then pie chart viewer should have a "slice (no value)" area
    And the legend of pie chart viewer should list 5 items
    When user sets "Include Nulls" property of pie chart viewer to "false"
    Then pie chart viewer should not have a "slice (no value)" area
    And the "slices" reading of pie chart viewer should be 4
    And the legend of pie chart viewer should list 4 items
    And the "angle value of Caucasian" reading of pie chart viewer should be 793
    And the "share of Caucasian" reading of pie chart viewer should be between 89.3 and 89.31
    And pie chart viewer should have repainted by at least 500 pixels
    And no errors should have been logged

  Scenario: Include Nulls on brings both back
    When user sets "Include Nulls" property of pie chart viewer to "true"
    Then pie chart viewer should have a "slice (no value)" area
    And the "slices" reading of pie chart viewer should be 5
    And the legend of pie chart viewer should list 5 items
    And the "share of (no value)" reading of pie chart viewer should be 11.2
    And pie chart viewer should have repainted by at least 500 pixels
    And no errors should have been logged

  Scenario: With Show Mouse Over Row Group off a hovered wedge paints no group overlay
    When user sets "Show Mouse Over Row Group" property of pie chart viewer to "false"
    And user moves the pointer away from pie chart viewer
    Then pie chart viewer should not have a "mouse over Caucasian" area
    When user hovers over the "slice Caucasian" area of pie chart viewer
    Then pie chart viewer should not have a "mouse over Caucasian" area
    And no errors should have been logged

  Scenario: With it on the hovered wedge lights up and goes out again
    When user moves the pointer away from pie chart viewer
    And user sets "Show Mouse Over Row Group" property of pie chart viewer to "true"
    And user hovers over the "slice Caucasian" area of pie chart viewer
    Then pie chart viewer should have a "mouse over Caucasian" area
    And pie chart viewer should not have a "mouse over Asian" area
    And pie chart viewer should have repainted
    When user hovers over the "slice Asian" area of pie chart viewer
    Then pie chart viewer should have a "mouse over Asian" area
    And pie chart viewer should not have a "mouse over Caucasian" area
    When user moves the pointer away from pie chart viewer
    Then pie chart viewer should not have a "mouse over Asian" area
    And no errors should have been logged

  Scenario: Back on the gapless column the blank wedge is gone for good
    When user sets "Category" property of pie chart viewer to "RACE"
    Then the "slices" reading of pie chart viewer should be 4
    And pie chart viewer should not have a "slice (no value)" area
    And the "angle value of Caucasian" reading of pie chart viewer should be 896
    When user removes "RACE_GAPS" column
    Then the table should not have a column "RACE_GAPS"
    And pie chart viewer should be painted
    And no errors should have been logged
