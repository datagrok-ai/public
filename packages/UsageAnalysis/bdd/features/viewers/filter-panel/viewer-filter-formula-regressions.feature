@viewers @realizes:viewers.filters
Feature: Viewer Filter formula regressions
  public#344: a viewer Filter formula that does not parse is not applied and breaks nothing. A valid
  formula filters the plot (27 Black rows on demog-1000) and an empty one gives every row back, with
  no error. An invalid formula is ignored as it should be (the scatter plot keeps all 1000 rows), but
  every time it is set the console gets a RangeError from DataFrameViewer.createCombinedRowFilter
  (core#1141). That scenario is a `@known-failure` pinned to its error claim; it passes, and asks for
  the tag to go, once the error is gone. AGE and WEIGHT have no blanks, so the plot draws all 1000
  rows. An invalid formula set after a valid one keeps the valid one's 27 rows (measured), so no
  gesture makes the ignored formula change a reading: its 1000-row claims read an unchanged plot.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | AGE    |
      | Y | WEIGHT |
    Then scatter plot viewer should show 1000 rows

  @realizes:public#344
  Scenario: A valid viewer Filter formula filters the plot and an empty one restores it (public#344)
    When user sets "Filter" property of scatter plot viewer to '${RACE} == "Black"'
    Then scatter plot viewer should show 27 rows
    And all rows should pass the filter
    When user sets "Filter" property of scatter plot viewer to ""
    Then scatter plot viewer should show 1000 rows
    And no errors should have been logged

  @known-failure:34 @realizes:public#344 @realizes:core#1141
  Scenario: An invalid viewer Filter formula is ignored and leaves the page working (public#344)
    When user sets "Filter" property of scatter plot viewer to '${RACE} = "Black"'
    Then all rows should pass the filter
    And scatter plot viewer should show 1000 rows
    And no errors should have been logged
