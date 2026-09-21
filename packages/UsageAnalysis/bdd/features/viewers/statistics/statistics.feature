@journey @viewers @realizes:viewers.stats-viewer @realizes:entities.viewer.action.close-viewer
Feature: The statistics table — one row per column, one column per aggregation
  What the viewer actually computes, read cell by cell.
  The viewer *is* a column grid over a derived frame whose rows are the source table's columns, so
  it publishes the inner grid's regions renamed after what they mean here: `row <COLUMN>` is one
  source column's whole row, `cell <stat> of <COLUMN>` is one number, `header <stat>` is an
  aggregation's header, and `<stat> of <COLUMN>` is the text that cell shows. The old spec proved
  the viewer had drawn by counting more than a thousand non-white pixels and proved a statistic had
  been added by diffing five hundred of them; every claim here is a number the viewer computed.
  Only the cells the grid scrolled into view are reported, columns included — a freshly added
  viewer is narrow enough to cut `med` and `stdev` off, so the Background gives it a size that
  shows all eight default aggregations. That is the fixture, not a workaround: an area for a cell
  the grid never drew would be a lie.
  Fixture: demog-1000, 11 columns, 1000 rows, HEIGHT blank in 128 of them.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a statistics viewer
    And user resizes statistics viewer to 900 by 400
    Then 1000 rows should pass the filter
    And the "stats" reading of statistics viewer should be "values, nulls, unique, min, max, avg, med, stdev"
    And the "columns" reading of statistics viewer should be 11
    And the "columns shown" reading of statistics viewer should be 11
    And the "rows shown" reading of statistics viewer should be 1000
    And statistics viewer should report no error

  Scenario: Every column of the table gets a row, and the count statistics are filled in for all of them
    Then the "values of USUBJID" reading of statistics viewer should be "1000"
    And the "values of AGE" reading of statistics viewer should be "1000"
    And the "values of SEX" reading of statistics viewer should be "1000"
    And the "unique of USUBJID" reading of statistics viewer should be "1000"
    And the "unique of SEX" reading of statistics viewer should be "2"
    And the "unique of RACE" reading of statistics viewer should be "4"
    And the "unique of STARTED" reading of statistics viewer should be "541"
    And the "nulls of AGE" reading of statistics viewer should be "0"
    And no errors should have been logged

  Scenario: A blank is counted out of values and into nulls
    Then the "values of HEIGHT" reading of statistics viewer should be "872"
    And the "nulls of HEIGHT" reading of statistics viewer should be "128"
    And the "values of WEIGHT" reading of statistics viewer should be "1000"
    And the "nulls of WEIGHT" reading of statistics viewer should be "0"
    And no errors should have been logged

  Scenario: The numerical aggregations are computed for a numerical column
    Then the "min of AGE" reading of statistics viewer should be "18.00"
    And the "max of AGE" reading of statistics viewer should be "89.00"
    And the "avg of AGE" reading of statistics viewer should be "45.68"
    And the "med of AGE" reading of statistics viewer should be "45.00"
    And the "stdev of AGE" reading of statistics viewer should be "13.45"
    And the "min of HEIGHT" reading of statistics viewer should be "137.32"
    And the "max of HEIGHT" reading of statistics viewer should be "198.86"
    And no errors should have been logged

  Scenario: For a categorical column the numerical aggregations are empty and the counts are not
    Then the "values of SEX" reading of statistics viewer should be "1000"
    And the "nulls of SEX" reading of statistics viewer should be "0"
    And the "unique of SEX" reading of statistics viewer should be "2"
    And the "min of SEX" reading of statistics viewer should be ""
    And the "max of SEX" reading of statistics viewer should be ""
    And the "avg of SEX" reading of statistics viewer should be ""
    And the "stdev of SEX" reading of statistics viewer should be ""
    And the "avg of RACE" reading of statistics viewer should be ""
    And the "avg of DIS_POP" reading of statistics viewer should be ""
    And no errors should have been logged

  Scenario: An aggregation has a header and a cell per column; the grid's own service columns have neither
    Then statistics viewer should have a "header values" area
    And statistics viewer should have a "header avg" area
    And statistics viewer should have a "header stdev" area
    And statistics viewer should have a "cell avg of AGE" area
    And statistics viewer should have a "cell avg of STARTED" area
    And statistics viewer should have a "row AGE" area
    And statistics viewer should have a "row SEVERITY" area
    And statistics viewer should not have a "header name" area
    And statistics viewer should not have a "cell name of AGE" area
    And the "row AGE" area of statistics viewer should be painted
    And no errors should have been logged

  Scenario: The title bar closes the statistics viewer, and adding it again computes the same table
    When user clicks on close icon of statistics viewer
    Then statistics viewer should be absent
    And the open tableview should have 0 statistics viewers
    When user adds a statistics viewer
    Then the open tableview should have 1 statistics viewer
    And the "columns shown" reading of statistics viewer should be 11
    And the "rows shown" reading of statistics viewer should be 1000
    And the "avg of AGE" reading of statistics viewer should be "45.68"
    And no errors should have been logged
