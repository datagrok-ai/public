@journey @viewers @realizes:viewers.calendar
Feature: The calendar filters the way the other viewers do
  With On Click = Filter a click keeps the table on the rows of that date, month or weekday, on top of
  the filters already there. The calendar re-applies it every time the table refilters, so a filter
  panel change keeps it, and lets go of it when it is closed, when On Click moves away from Filter or
  when the filter is reset (core#233). A click on the same part takes the filter off, a click on
  another part replaces it. On demog-1000 January 1990 holds 32 rows and the Sundays 153.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a calendar viewer
    Then the "date column" reading of calendar viewer should be "STARTED"
    And the "rows shown" reading of calendar viewer should be 1000
    And calendar viewer should be painted

  @realizes:core#233
  Scenario: Closing the calendar lets go of its filter
    When user sets "onClick" property of calendar viewer to "Filter"
    And user clicks on the "month 1990-01" area of calendar viewer
    Then 32 rows should pass the filter
    When user clicks on close icon of calendar viewer
    Then 1000 rows should pass the filter
    And no errors should have been logged

  @realizes:core#233
  Scenario: A second click takes the filter off, another part replaces it
    When user adds a calendar viewer
    And user sets "onClick" property of calendar viewer to "Filter"
    And user clicks on the "month 1990-01" area of calendar viewer
    Then 32 rows should pass the filter
    When user clicks on the "month 1990-01" area of calendar viewer
    Then 1000 rows should pass the filter
    When user clicks on the "month 1990-01" area of calendar viewer
    Then 32 rows should pass the filter
    When user clicks on the "weekday Sunday" area of calendar viewer
    Then 153 rows should pass the filter
    And no errors should have been logged

  @realizes:core#233
  Scenario: Leaving On Click = Filter and resetting the filter let go of it
    When user sets "onClick" property of calendar viewer to "Select"
    Then 1000 rows should pass the filter
    When user sets "onClick" property of calendar viewer to "Filter"
    And user clicks on the "month 1990-01" area of calendar viewer
    Then 32 rows should pass the filter
    When user picks "View > Reset Filter" from the top menu
    Then 1000 rows should pass the filter
    When user clicks on the "month 1990-01" area of calendar viewer
    Then 32 rows should pass the filter
    And no errors should have been logged

  @realizes:core#233
  Scenario: The calendar's filter and the filter panel apply together, whichever comes first
    When user adds a categorical filter on "SEX" keeping "M"
    Then 13 rows should pass the filter
    When user picks "View > Reset Filter" from the top menu
    And user clicks on the "month 1990-01" area of calendar viewer
    And user adds a categorical filter on "SEX" keeping "M"
    Then 13 rows should pass the filter
    And no errors should have been logged
