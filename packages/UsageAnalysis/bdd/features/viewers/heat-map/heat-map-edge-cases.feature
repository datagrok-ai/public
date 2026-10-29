@viewers @edge-cases @realizes:viewers.heat-map
Feature: Heat map edge cases — every row deleted, a column removed
  The heat map on data changes its other features do not reach (a table filtered down to nothing:
  `viewer-rules-edge-cases.feature`). A table whose rows are all deleted logs no error. A column removed from the table while the map shows it
  leaves the map with its band gone and the other columns kept, as a grid shows the table's
  columns. demog-1000, whose map shows its 11 columns.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a heat map viewer
    Then the "rows shown" reading of heat map viewer should be 1000
    And the "columns shown" reading of heat map viewer should be 11
    And heat map viewer should be painted

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

  Scenario: An Alt-click without moving zooms the map to one row, an Alt-drag to the rows it spans
    Then the "y scroll span" reading of heat map viewer should be 1
    When user clicks on the "column AGE" area of heat map viewer holding Alt
    Then the "y scroll span" reading of heat map viewer should be 0
    And the "row height" reading of heat map viewer should be higher than before
    And no rows should be selected
    When user double-clicks on the "y scroll slider" area of heat map viewer
    Then the "y scroll span" reading of heat map viewer should be 1
    When user drags a zoom box over the "column AGE" area of heat map viewer
    Then the "y scroll span" reading of heat map viewer should be between 0.1 and 0.99
    And no rows should be selected
    And no errors should have been logged
