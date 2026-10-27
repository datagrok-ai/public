@viewers @realizes:viewers.filters
Feature: Viewer Filter formula regressions
  public#344: a viewer Filter formula that does not parse is not applied and breaks nothing. On the
  current master the formula is ignored as it should be (the scatter plot keeps all 1000 rows), but
  every time it is set the console gets a RangeError from DataFrameViewer.createCombinedRowFilter
  (core#1141). The scenario is a `@known-failure` pinned to the first error claim;
  it passes, and asks for the tag to go, once the error is gone. A valid formula then filters the
  plot (27 Black rows on demog-1000) and an empty one gives every row back. AGE and WEIGHT have no
  blanks, so the plot draws all 1000 rows.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | AGE    |
      | Y | WEIGHT |
    Then scatter plot viewer should show 1000 rows

  @known-failure:24 @realizes:public#344 @realizes:core#1141
  Scenario: An invalid viewer Filter formula leaves the page working (public#344)
    When user sets "Filter" property of scatter plot viewer to '${RACE} = "Black"'
    Then all rows should pass the filter
    And scatter plot viewer should show 1000 rows
    And no errors should have been logged
    When user sets "Filter" property of scatter plot viewer to '${RACE} == "Black"'
    Then scatter plot viewer should show 27 rows
    And all rows should pass the filter
    When user sets "Filter" property of scatter plot viewer to ""
    Then scatter plot viewer should show 1000 rows
    And no errors should have been logged
