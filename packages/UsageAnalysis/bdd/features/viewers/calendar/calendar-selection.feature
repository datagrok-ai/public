@journey @viewers @realizes:viewers.calendar
Feature: Calendar selection, On Click and Show Filtered Only
  What a click on a day, a month band or a weekday header takes, what Shift and Ctrl do to it,
  what On Click = Filter does instead, and what Show Filtered Only moves.
  Every expected count is the calendar's own reading of the same target, so no scenario
  re-implements a group-by: `rows of weekday Sunday` is 153 and clicking that header selects 153
  rows. That equality is what GROK-20634 was about — Dart's `DateTime.weekday` is 1..7 with
  Sunday = 7 while the calendar draws Sunday first, and the hit test now goes through
  `x.weekday % 7` (calendar_core.dart:188), so the Sunday header no longer matches nothing.
  Show Filtered Only is read the same way: `_shownCounts` picks between the counts of every dated
  row and the counts of the filtered ones (calendar_core.dart:271), so turning it off under a
  filter puts `rows shown` back to 1000 while `dataFrame.filter` stands at 447 — the honest
  witness for GROK-20635, which the spec this replaces could only assert as "the canvas changed".

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a calendar viewer
    Then the "date column" reading of calendar viewer should be "STARTED"
    And the "rows shown" reading of calendar viewer should be 1000
    And "On Click" property of calendar viewer should be "Select"
    And calendar viewer should be painted

  Scenario: Clicking a day selects exactly the rows dated that day
    Given user clears the row selection
    When user clicks on the "day 1989-12-21" area of calendar viewer
    Then 5 rows should be selected
    And every selected row should pass the filter
    When user clears the row selection
    And user clicks on the "busiest day" area of calendar viewer
    Then some rows should be selected
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged

  Scenario: Clicking a weekday header selects every row that falls on it, Sunday included
    Given user clears the row selection
    When user clicks on the "weekday Sunday" area of calendar viewer
    Then 153 rows should be selected
    When user clears the row selection
    And user clicks on the "weekday Saturday" area of calendar viewer
    Then 117 rows should be selected
    When user clears the row selection
    And user clicks on the "weekday Thursday" area of calendar viewer
    Then 154 rows should be selected
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged

  Scenario: Shift extends the weekday selection and Ctrl toggles it back off
    Given user clears the row selection
    When user clicks on the "weekday Monday" area of calendar viewer
    Then 151 rows should be selected
    When user clicks on the "weekday Tuesday" area of calendar viewer holding Shift
    Then 287 rows should be selected
    When user clicks on the "weekday Tuesday" area of calendar viewer holding Control
    Then 151 rows should be selected
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged

  Scenario: Clicking a month band selects that month, and On Click = Filter filters it instead
    Given user clears the row selection
    When user clicks on the "month 1990-01" area of calendar viewer
    Then 32 rows should be selected
    And 1000 rows should pass the filter
    When user clears the row selection
    And user sets "onClick" property of calendar viewer to "Filter"
    And user clicks on the "month 1990-01" area of calendar viewer
    Then 32 rows should pass the filter
    And no rows should be selected
    And calendar viewer should have repainted
    When user resets the filter
    And user sets "onClick" property of calendar viewer to "Select"
    Then 1000 rows should pass the filter
    And "On Click" property of calendar viewer should be "Select"
    And the "rows shown" reading of calendar viewer should be 1000
    And no errors should have been logged

  Scenario: Show Filtered Only decides whether the discs count the filtered rows or all of them
    Then the "rows shown" reading of calendar viewer should be 1000
    When user adds a categorical filter on "SEX" keeping "M"
    Then 447 rows should pass the filter
    And the "rows shown" reading of calendar viewer should be 447
    And the "rows of weekday Sunday" reading of calendar viewer should be lower than before
    And the "days drawn" reading of calendar viewer should be lower than before
    And calendar viewer should have repainted
    When user sets "showFilteredOnly" property of calendar viewer to "false"
    Then the "rows shown" reading of calendar viewer should be 1000
    And the "rows of weekday Sunday" reading of calendar viewer should be 153
    And the "days drawn" reading of calendar viewer should be higher than before
    And 447 rows should pass the filter
    And calendar viewer should have repainted
    When user sets "showFilteredOnly" property of calendar viewer to "true"
    Then the "rows shown" reading of calendar viewer should be 447
    When user hovers over "SEX" filter card
    And user clicks on close of "SEX" filter card
    Then 1000 rows should pass the filter
    And the "rows shown" reading of calendar viewer should be 1000
    And the "rows of weekday Sunday" reading of calendar viewer should be 153
    And no errors should have been logged
