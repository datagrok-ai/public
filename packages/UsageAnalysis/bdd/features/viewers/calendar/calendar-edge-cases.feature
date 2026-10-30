@viewers @realizes:viewers.calendar
Feature: Calendar edge cases
  Guards for the calendar at the edges of its data and its lifecycle, on demog-1000 (32 rows start
  in January 1990). Every row deleted is a row of `viewer-rules-edge-cases.feature`. Closing a calendar whose On Click = Filter click is filtering the table (the open defect
  core#233) is a row of `viewer-rules-edge-cases.feature`.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a calendar viewer
    Then the "rows shown" reading of calendar viewer should be 1000

  @edge-cases @known-failure:23 @realizes:core#1221
  Scenario: The calendar's click filter stays when a category is toggled in the filter panel, and goes when On Click is switched to Select
    The category click releases the calendar's filter: 553 rows pass (every F row), not the 19 F rows
    of January 1990.
    When user sets "onClick" property of calendar viewer to "Filter"
    And user clicks on the "month 1990-01" area of calendar viewer
    Then 32 rows should pass the filter
    When user opens an empty filter panel
    And user adds a card for "SEX" to the filter panel
    And user clicks on the "category F of SEX" area of filter panel
    Then 19 rows should pass the filter
    When user sets "onClick" property of calendar viewer to "Select"
    Then 553 rows should pass the filter
    And no errors should have been logged

  @edge-cases @known-failure:41 @realizes:core#1222
  Scenario: A calendar over a date column whose values are all empty shows a "No data" message instead of a grid
    The calendar reports no message.
    Given user opens a table "empty-dates" with:
      | When       | Value |
      | 2024-01-05 | 1     |
      | 2024-01-06 | 2     |
    When user sets "When" column in row 1 to ""
    And user sets "When" column in row 2 to ""
    Then "When" column should have type "datetime"
    When user adds a calendar viewer
    Then the "date column" reading of calendar viewer should be "When"
    And the "rows shown" reading of calendar viewer should be 0
    And calendar viewer should report the error "No data"
    And calendar viewer should not have a "days" area
    And no errors should have been logged

  @edge-cases @known-failure:49 @realizes:core#1222
  Scenario: Removing the date column of an open calendar shows a "No data" message instead of a grid
    The calendar reports no message.
    When user removes "STARTED" column
    Then calendar viewer should report the error "No data"
    And calendar viewer should not have a "days" area
    And no errors should have been logged
