@viewers @edge-cases @realizes:viewers.heat-map
Feature: Heat map edge cases — an empty row set, every row deleted, a column removed
  The heat map on data changes its other features do not reach. A table filtered down to nothing
  leaves a map with no row and no error, and the whole table comes back with the filter; a table
  whose rows are all deleted logs no error. A column removed from the table while the map shows it
  leaves the map with its band gone and the other columns kept, as a grid shows the table's
  columns. demog-1000, whose map shows its 11 columns.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a heat map viewer
    Then the "rows shown" reading of heat map viewer should be 1000
    And the "columns shown" reading of heat map viewer should be 11
    And heat map viewer should be painted

  Scenario: With every row filtered out the map shows no row and logs no error
    When user filters rows where "AGE" is between 200 and 300
    Then 0 rows should pass the filter
    And the "rows shown" reading of heat map viewer should be 0
    And heat map viewer should report no error
    And no errors should have been logged
    When user resets the filter
    Then the "rows shown" reading of heat map viewer should be 1000

  Scenario: Once every row of the table is deleted the map logs no error
    When user selects all rows
    And user deletes the selected rows
    Then the table should have 0 rows
    And the "rows" reading of heat map viewer should be 0
    And heat map viewer should report no error
    And no errors should have been logged

  Scenario: A column removed from the table leaves the map without its band and keeps the others
    Then heat map viewer should have a "column AGE" area
    When user removes "AGE" column
    Then the table should not have a column "AGE"
    And the "columns shown" reading of heat map viewer should be 10
    And heat map viewer should not have a "column AGE" area
    And heat map viewer should have a "column WEIGHT" area
    And heat map viewer should report no error
    And no errors should have been logged
