@viewers @realizes:viewers.calendar
Feature: Calendar known failures
  Temporary: the known failures of the calendar while their tickets are open. Each scenario is
  expected to fail on its pinned step; when its ticket is fixed, remove the tag and move the scenario
  back into calendar-selection.feature as an ordinary @edge-cases scenario.

  Each scenario is a test of its own on demog-1000 with a calendar of STARTED.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a calendar viewer
    Then the "date column" reading of calendar viewer should be "STARTED"
    And "On Click" property of calendar viewer should be "Select"
    And all rows should pass the filter

  @edge-cases @known-failure:25 @realizes:core#1221
  Scenario: The calendar's click filter stays when a category is toggled in the filter panel (core#1221)
    When user sets "onClick" property of calendar viewer to "Filter"
    And user clicks on the "month 1990-01" area of calendar viewer
    Then 32 rows should pass the filter
    When user opens an empty filter panel
    And user adds a card for "SEX" to the filter panel
    And user clicks on the "category F of SEX" area of filter panel
    Then 19 rows should pass the filter
