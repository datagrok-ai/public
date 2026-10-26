@viewers @edge-cases @realizes:viewers.stats-viewer
Feature: Statistics edge cases: no rows counted, one row, columns added, removed and renamed
  Guards for the data shapes and the changes under an open statistics viewer. The statistics are
  computed over the rows the Row Source lets through: a filter that removes every row, and Row
  Source Selected with nothing selected, both leave zero rows counted, so every column has 0 values,
  and the viewer logs no error. A one-row table counts its one value. With Columns left empty the
  viewer shows a row for every column of the table (`stats_viewer_core.dart`, a null list means all
  columns), so a calculated column added while it is open gets a row and a removed column loses its
  row. A renamed column is followed by the viewer's settings (the viewer base renames column
  settings on the table's rename event), so its row is shown under the new name with the same
  numbers. What the numeric statistics (avg, min, stdev) should read over no values is not defined:
  see MISSING.md. Data: demog-1000 (11 columns, AGE averages 45.68 over 1000 rows).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a statistics viewer
    And user resizes statistics viewer to 900 by 400
    Then the "rows shown" reading of statistics viewer should be 1000
    And the "columns shown" reading of statistics viewer should be 11
    And the "avg of AGE" reading of statistics viewer should be "45.68"

  Scenario: A filter that removes every row leaves no values counted and no error
    When user filters rows where "AGE" is between 500 and 600
    Then 0 rows should pass the filter
    And the "rows shown" reading of statistics viewer should be 0
    And the "values of AGE" reading of statistics viewer should be "0"
    And the "values of RACE" reading of statistics viewer should be "0"
    And statistics viewer should report no error
    When user resets the filter
    Then the "rows shown" reading of statistics viewer should be 1000
    And the "avg of AGE" reading of statistics viewer should be "45.68"
    And no errors should have been logged

  Scenario: Row Source Selected with nothing selected counts no rows without an error
    Given user selects no rows
    When user sets "rowSource" property of statistics viewer to "Selected"
    Then the "rows shown" reading of statistics viewer should be 0
    And the "values of AGE" reading of statistics viewer should be "0"
    And statistics viewer should report no error
    When user sets "rowSource" property of statistics viewer to "Filtered"
    Then the "rows shown" reading of statistics viewer should be 1000
    And no errors should have been logged

  Scenario: A one-row table counts its one value
    Given user opens a table "one-row" with:
      | name | x  |
      | a    | 42 |
    And user adds a statistics viewer
    Then the "rows shown" reading of statistics viewer should be 1
    And the "values of x" reading of statistics viewer should be "1"
    And the "nulls of x" reading of statistics viewer should be "0"
    And no errors should have been logged

  Scenario: A calculated column added while the viewer is open gets a row
    When user adds a calculated column "AGE2" with formula "${AGE} * 2"
    Then the "columns shown" reading of statistics viewer should be 12
    And the "values of AGE2" reading of statistics viewer should be "1000"
    And no errors should have been logged

  Scenario: A removed column loses its row without an error
    Then statistics viewer should have a "row RACE" area
    When user removes "RACE" column
    Then the "columns shown" reading of statistics viewer should be 10
    And statistics viewer should not have a "row RACE" area
    And the "avg of AGE" reading of statistics viewer should be "45.68"
    And no errors should have been logged

  Scenario: A renamed column keeps its row under the new name
    When user renames "AGE" column to "AGE_YEARS"
    Then the "columns shown" reading of statistics viewer should be 11
    And the "avg of AGE_YEARS" reading of statistics viewer should be "45.68"
    And statistics viewer should not have a "row AGE" area
    And no errors should have been logged
