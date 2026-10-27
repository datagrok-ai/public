@viewers @realizes:viewers.bar-chart
Feature: Bar chart edge cases
  Guards for the bar chart at the edges of its data and settings, on demog-1000 (RACE has four
  categories) unless a scenario builds its own table. Each one is a crash or a wrong display the
  bar chart audit core#974 found and fixed: a vertical chart over a table with no rows, a negative
  Bar Corner Radius, a category whose name holds markup characters shown in the tooltip as HTML, a
  Filter-click on a chart stacked by a date column lost by an unrelated setting change, and Show
  Filtered Rows with every row filtered out. Expected: no error, the tooltip shows the category as
  typed, the click on empty space releases the filter (core#974), and with Row Source All every bar
  stays (help bar-chart.md, Show Filtered Rows).

  Background:
    Given user is logged in

  @edge-cases @realizes:core#974
  Scenario: A vertical bar chart over a table with every row deleted draws no bars and logs no error
    Given user opens demog-1000 dataset
    And user adds a bar chart viewer with:
      | Split       | RACE     |
      | Orientation | vertical |
    Then the "bars" reading of bar chart viewer should be 4
    When user selects all rows
    And user deletes the selected rows
    Then the table should have 0 rows
    And the "bars" reading of bar chart viewer should be 0
    And no errors should have been logged

  @edge-cases @realizes:core#974
  Scenario: A negative Bar Corner Radius draws the bars and logs no error
    Given user opens demog-1000 dataset
    And user adds a bar chart viewer with:
      | Split | RACE |
      | Stack | SEX  |
    When user sets "Bar Corner Radius" property of bar chart viewer to "-5"
    Then the "bars" reading of bar chart viewer should be 4
    And the "bar Caucasian" area of bar chart viewer should be painted
    And no errors should have been logged

  @edge-cases @realizes:core#974
  Scenario: The tooltip shows a category with markup characters as typed
    Given user opens a table "markup-categories" with:
      | Category | Value |
      | a<b      | 1     |
      | a<b      | 2     |
      | x & y    | 3     |
    And user adds a bar chart viewer with:
      | Split | Category |
    Then the "bars" reading of bar chart viewer should be 2
    When user hovers over the "bar a<b" area of bar chart viewer
    Then tooltip should contain text "a<b"
    And tooltip should contain text "2 rows"
    When user hovers over the "bar x & y" area of bar chart viewer
    Then tooltip should contain text "x & y"
    And no errors should have been logged

  @edge-cases @realizes:core#974
  Scenario: A Filter click on a chart stacked by a date column is released after an unrelated setting changes
    Given user opens demog-1000 dataset
    And user adds a bar chart viewer with:
      | Split    | RACE    |
      | Stack    | STARTED |
      | On Click | Filter  |
    When user clicks on the "bar Asian" area of bar chart viewer
    Then 6 rows should pass the filter
    When user sets "Bar Corner Radius" property of bar chart viewer to "4"
    And user clicks on empty plot space of bar chart viewer
    Then all rows should pass the filter
    And no errors should have been logged

  @edge-cases @realizes:core#974
  Scenario: With Row Source All and Show Filtered Rows every bar stays when every row is filtered out
    Given user opens demog-1000 dataset
    And user adds a bar chart viewer with:
      | Split               | RACE |
      | Row Source          | All  |
      | Show Filtered Rows  | true |
    Then the "bars" reading of bar chart viewer should be 4
    When user filters rows where "AGE" is between 200 and 300
    Then 0 rows should pass the filter
    And the "bars" reading of bar chart viewer should be 4
    And no errors should have been logged

  @edge-cases
  Scenario: Closing the bar chart gives back the rows its clicked bar filtered out
    Given user opens demog-1000 dataset
    And user adds a bar chart viewer with:
      | Split    | RACE   |
      | On Click | Filter |
    When user clicks on the "bar Asian" area of bar chart viewer
    Then fewer than 1000 rows should pass the filter
    When user clicks on close icon of bar chart viewer
    Then the open tableview should have 0 bar chart viewers
    And all rows should pass the filter
    And no errors should have been logged
