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
  Scenario Outline: <caption> typed as <typed> in the context panel is not taken by the bar chart
    Given user opens demog-1000 dataset
    And user adds a bar chart viewer with:
      | Split | RACE |
    When user clicks on settings icon of bar chart viewer
    Then context panel should be visible
    Given "Style" category in context panel is expanded
    When user enters "<typed>" into "<caption>" property in context panel
    Then "<caption>" property of bar chart viewer should be "<kept>"
    And bar chart viewer should report no error
    And no errors should have been logged

    @known-failure:92 @realizes:core#1200
    Examples: The typed value is written into the property; the declarations set no lower bound
      | caption              | typed | kept |
      | Max Bar Height       | 0     | 50   |
      | Max Bar Height       | -5    | 50   |
      | Max Category Width   | 0     | 300  |
      | Category Value Width | 0     | 50   |

  @edge-cases
  Scenario Outline: "abc" typed into the numeric <caption> property is not taken and the old value stays
    Given user opens demog-1000 dataset
    And user adds a bar chart viewer with:
      | Split | RACE |
    When user clicks on settings icon of bar chart viewer
    Then context panel should be visible
    Given "<category>" category in context panel is expanded
    When user enters "abc" into "<caption>" property in context panel
    Then "<caption>" property of bar chart viewer should be "<kept>"
    And the "bars" reading of bar chart viewer should be 4
    And no errors should have been logged

    Examples:
      | category | caption            | kept |
      | Y Axis   | Min                |      |
      | Style    | Max Bar Height     | 50   |
      | Style    | Max Category Width | 300  |

  @edge-cases
  Scenario: Zoomed to two bars, then filtered to one category, the bar chart resets the zoom and shows the remaining bar
    Given user opens demog-1000 dataset
    And user adds a bar chart viewer with:
      | Split | RACE |
    Then the "bars" reading of bar chart viewer should be 4
    When user zooms into the categories from the "bar Asian" area to the "bar Black" area of bar chart viewer
    Then the "bars" reading of bar chart viewer should be 2
    When user filters rows where "RACE" is "Other"
    Then the "bars" reading of bar chart viewer should be 1
    And bar chart viewer should have a "bar Other" area
    When user resets the filter
    Then the "bars" reading of bar chart viewer should be 4
    And no errors should have been logged

  @edge-cases
  Scenario Outline: Sorted by value in <order> order, the category whose aggregate is empty has a bar only with Show Empty Bars, as the lowest value
    Given user opens a table "empty-aggregate" with:
      | Category | Value |
      | A        | 1     |
      | A        | 3     |
      | B        | 5     |
      | B        | 7     |
      | C        |       |
      | C        |       |
    And user adds a bar chart viewer with:
      | Split           | Category   |
      | Value           | Value      |
      | Value Aggr Type | avg        |
      | Orientation     | horizontal |
      | Bar Sort Type   | by value   |
      | Bar Sort Order  | <order>    |
    Then the "bars" reading of bar chart viewer should be 2
    And bar chart viewer should not have a "bar C" area
    When user sets "showEmptyBars" property of bar chart viewer to "true"
    Then the "bars" reading of bar chart viewer should be 3
    And bar chart viewer should have a "bar C" area
    And the bars of bar chart viewer should lie one under another
    And the "bar C" area of bar chart viewer should lie <where> the "bar A" area
    And the "bar C" area of bar chart viewer should lie <where> the "bar B" area
    And no errors should have been logged

    Examples:
      | order | where |
      | asc   | above |
      | desc  | below |

  @edge-cases
  Scenario: Removing the Stack column stops the stacking and logs no error
    Given user opens demog-1000 dataset
    And user adds a bar chart viewer with:
      | Split | RACE |
      | Stack | SEX  |
    Then the "stack segments" reading of bar chart viewer should be at least 1
    When user removes "SEX" column
    Then the "stack segments" reading of bar chart viewer should be 0
    And the "bars" reading of bar chart viewer should be 4
    And bar chart viewer should report no error
    And no errors should have been logged
