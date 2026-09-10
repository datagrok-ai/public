@journey @viewers @realizes:viewers.stats-viewer
Feature: Statistics for a date column, and the aggregations added for one from the submenu
  These scenarios come from `statistics-ui.md`, the manual-only file — "Statistics for date
  columns", "Add min and max for a date column via submenu" and "Column visibility". They were
  manual because the only thing worth claiming about a date statistic is what the cell *shows*: the
  stored value is a float, µs since the epoch, and `_formatDateTimeStatCell` turns a date-valued
  aggregation into a date, `stdev` — a spread, not a point in time — into a duration. Nothing but a
  screenshot could see that. `<stat> of <COLUMN>` is meant to be exactly that text, which is what
  makes them automatable; the last scenario here is where it is not yet.
  "Background color" stays manual: it is a colour of the viewer's own backdrop and nothing reports
  it.
  Fixture: demog-1000, STARTED running 1989-12-03 to 1991-11-30 over 541 distinct days, no blanks.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a statistics viewer
    And user resizes statistics viewer to 900 by 400
    Then the "columns shown" reading of statistics viewer should be 11
    And statistics viewer should have a "row STARTED" area
    And statistics viewer should report no error

  Scenario: The count statistics are filled in for a date column
    Then the "values of STARTED" reading of statistics viewer should be "1000"
    And the "nulls of STARTED" reading of statistics viewer should be "0"
    And the "unique of STARTED" reading of statistics viewer should be "541"
    And no errors should have been logged

  Scenario: A date column gets the value aggregations too, and they are not empty
    Then the "min of STARTED" reading of statistics viewer should not be ""
    And the "max of STARTED" reading of statistics viewer should not be ""
    And the "avg of STARTED" reading of statistics viewer should not be ""
    And the "min of USUBJID" reading of statistics viewer should be ""
    And no errors should have been logged

  Scenario: min and max are added for a date column from the Statistics submenu
    When user sets "stats" property of statistics viewer to "values, nulls, unique"
    Then the "stats" reading of statistics viewer should be "values, nulls, unique"
    And statistics viewer should not have a "header min" area
    And statistics viewer should not have a "header max" area
    When user picks "Statistics > min" from the context menu of the "row STARTED" area of statistics viewer
    Then the "stats" reading of statistics viewer should contain "min"
    And statistics viewer should have a "header min" area
    And the "min of STARTED" reading of statistics viewer should not be ""
    When user picks "Statistics > max" from the context menu of the "row STARTED" area of statistics viewer
    Then the "stats" reading of statistics viewer should contain "max"
    And statistics viewer should have a "header max" area
    And the "max of STARTED" reading of statistics viewer should not be ""
    And the "min of STARTED" and "max of STARTED" readings of statistics viewer should differ
    When user sets "stats" property of statistics viewer to "values, nulls, unique, min, max, avg, med, stdev"
    Then the "stats" reading of statistics viewer should be "values, nulls, unique, min, max, avg, med, stdev"
    And no errors should have been logged

  Scenario: Column visibility — a column dropped from Columns loses its row and gets it back
    Then statistics viewer should have a "row HEIGHT" area
    And the "values of HEIGHT" reading of statistics viewer should be "872"
    When user sets "columnNames" property of statistics viewer to "USUBJID, AGE, SEX, RACE, DIS_POP, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    Then the "columns shown" reading of statistics viewer should be 10
    And statistics viewer should not have a "row HEIGHT" area
    And statistics viewer should have a "row STARTED" area
    When user sets "columnNames" property of statistics viewer to "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    Then the "columns shown" reading of statistics viewer should be 11
    And statistics viewer should have a "row HEIGHT" area
    And the "values of HEIGHT" reading of statistics viewer should be "872"
    And the "nulls of HEIGHT" reading of statistics viewer should be "128"
    And no errors should have been logged

  Scenario: A date statistic reads as the date the cell shows
    # This scenario is why the file it came from was manual-only. A cell's displayed text comes
    # from the `GRID_CELL_PREPARE` listeners — here `_formatDateTimeStatCell`
    # (`stats_viewer_core.dart:26,42-56`) — which used to run only while the grid drew, so a status
    # reading built its own `GridCell` outside the render pass and got the raw value: the viewer
    # showed "12/3/1989" while the reading said "628646377160704.00", µs since the epoch. The grid
    # status fires the prepare event before it reads (`grid_status.dart`), so the reading is the
    # text the grid draws — for every viewer that formats through `onCellPrepare`, not just this one.
    Then the "min of STARTED" reading of statistics viewer should be "12/3/1989"
    And the "max of STARTED" reading of statistics viewer should be "11/30/1991"
