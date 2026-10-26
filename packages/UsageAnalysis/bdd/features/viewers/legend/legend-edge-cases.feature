@viewers @edge-cases @realizes:viewers.legend
Feature: Legend edge cases: the color column removed under a legend selection, every row filtered out, a constant column
  Guards for the data shapes and the changes under a scatter plot's legend. A category selected in
  the legend narrows the plot to its rows; when the color column it was made on is removed, the
  selection goes with it and the plot draws every row again (the selection belongs to the column,
  `legend-interaction.feature`; a removed column is dropped from the viewer's settings). A legend
  selection outlives a filter that removes every row: once the rows are back, the legend lists its
  four categories and the selected one still narrows the plot. A color column holding one value
  gives a legend of one item. Data: demog-1000, X = WEIGHT, Y = HEIGHT (872 rows drawn), colored by
  RACE (Black 25 rows drawn).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X                 | WEIGHT |
      | Y                 | HEIGHT |
      | Color             | RACE   |
      | Legend Position   | Right  |
    Then the legend of scatter plot viewer should list 4 items
    And scatter plot viewer should show 872 rows

  Scenario: Removing the color column under a legend selection gives every row back
    When user clicks on "Black" item in the legend of scatter plot viewer
    Then scatter plot viewer should show 25 rows
    When user removes "RACE" column
    Then "Color" property of scatter plot viewer should not be "RACE"
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  Scenario: A legend selection outlives a filter that removes every row
    When user clicks on "Black" item in the legend of scatter plot viewer
    Then scatter plot viewer should show 25 rows
    When user filters rows where "AGE" is between 500 and 600
    Then 0 rows should pass the filter
    And scatter plot viewer should show 0 rows
    When user resets the filter
    Then the legend of scatter plot viewer should list 4 items
    And "Black" legend item in legend of scatter plot viewer should be selected
    And scatter plot viewer should show 25 rows
    And no errors should have been logged

  Scenario: A color column holding one value gives a legend of one item
    When user adds a calculated column "ONE" with formula "'same'"
    And user sets "Color" property of scatter plot viewer to "ONE"
    Then the legend of scatter plot viewer should list 1 item
    And scatter plot viewer should show 872 rows
    And no errors should have been logged
