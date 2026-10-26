@viewers @realizes:viewers.calendar
Feature: Calendar edge cases
  Guards for the calendar at the edges of its data and its lifecycle, on demog-1000 (32 rows start
  in January 1990). A calendar over a table whose rows are all deleted counts nothing and logs no
  error. Closing a calendar whose On Click = Filter click is filtering the table is the open defect
  core#233.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a calendar viewer
    Then the "rows shown" reading of calendar viewer should be 1000

  @edge-cases
  Scenario: A calendar over a table with every row deleted counts no rows and logs no error
    When user selects all rows
    And user deletes the selected rows
    Then the table should have 0 rows
    And the "rows shown" reading of calendar viewer should be 0
    And calendar viewer should not have a "busiest day" area
    And no errors should have been logged

  @edge-cases @known-failure:33 @realizes:core#233
  Scenario: Closing a calendar releases the filter its click applied
    core#233: a click with On Click = Filter writes the table filter, and closing the calendar leaves
    the table filtered until something else asks for a refilter. Remove the tag when core#233 is
    fixed and this scenario passes.
    When user sets "onClick" property of calendar viewer to "Filter"
    And user clicks on the "month 1990-01" area of calendar viewer
    Then 32 rows should pass the filter
    When user clicks on close icon of calendar viewer
    Then calendar viewer should be absent
    And all rows should pass the filter
