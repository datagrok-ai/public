@journey @viewers @realizes:viewers.box-plot
Feature: Box plot hover highlights and date categories
  What the box plot does under the pointer and with a date category: hovering a category's label
  on the category axis highlights that row group and Show Mouse Over Row Group switches it off
  (box-plot-ui "Mouse-over highlight"); a hover on a bar of another viewer keeps the hovered group's
  own categorical marker colours, and with Show Mouse Over Row Group off leaves the box plot
  unchanged (github-3065; that the other groups fade has no claim yet, MISSING.md); a date column mapped to month or quarter splits into
  one category per month (January..December) or quarter (Q1..Q4). One journey on demog-1000 with a
  box plot of AGE by RACE and 10-pixel markers; every scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a box plot viewer with:
      | Value       | AGE  |
      | Category 1  | RACE |
      | Marker Size | 10   |

  Scenario: Hovering a category's axis label highlights its row group, and Show Mouse Over Row Group switches it off
    When user takes a snapshot of box plot viewer
    And user hovers over the "category Caucasian" area of box plot viewer
    Then box plot viewer should have repainted by at least 2000 pixels
    When user moves the pointer away from grid
    And user sets "Show Mouse Over Row Group" property of box plot viewer to "false"
    And user takes a snapshot of box plot viewer
    And user hovers over the "category Caucasian" area of box plot viewer
    Then box plot viewer should not have repainted
    When user moves the pointer away from grid
    And user sets "Show Mouse Over Row Group" property of box plot viewer to "true"
    Then no errors should have been logged

  Scenario: A hover on another viewer's bar keeps the hovered group's own marker colours, and Show Mouse Over Row Group switches it off (github-3065)
    When user adds a bar chart viewer with:
      | Split | RACE |
    And user sets "Marker Color Column" property of box plot viewer to "RACE"
    Then the "Caucasian values" area of box plot viewer should contain the color "#2CA02C"
    When user takes a snapshot of box plot viewer
    And user hovers over the "bar Caucasian" area of bar chart viewer
    Then box plot viewer should have repainted by at least 2000 pixels
    And the "Caucasian values" area of box plot viewer should contain the color "#2CA02C"
    When user moves the pointer away from grid
    And user sets "Show Mouse Over Row Group" property of box plot viewer to "false"
    And user takes a snapshot of box plot viewer
    And user hovers over the "bar Caucasian" area of bar chart viewer
    Then box plot viewer should not have repainted
    When user moves the pointer away from grid
    And user sets "Show Mouse Over Row Group" property of box plot viewer to "true"
    And user clicks on close icon of bar chart viewer
    Then bar chart viewer should be absent
    When user sets "Marker Color Column" property of box plot viewer to ""
    Then no errors should have been logged

  Scenario: A date category mapped to month and to quarter splits into months and quarters
    When user sets "Category 1" property of box plot viewer to "STARTED"
    And user sets "Category 1 Map" property of box plot viewer to "month"
    Then box plot viewer should have a "category January" area
    And box plot viewer should have a "category December" area
    And box plot viewer should have a "June values" area
    And box plot viewer should not have a "category Q1" area
    When user sets "Category 1 Map" property of box plot viewer to "quarter"
    Then box plot viewer should have a "category Q1" area
    And box plot viewer should have a "category Q4" area
    And box plot viewer should have a "Q3 values" area
    And box plot viewer should not have a "category January" area
    When user sets properties of box plot viewer:
      | Category 1 Map |      |
      | Category 1     | RACE |
    Then box plot viewer should have a "category Caucasian" area
    And box plot viewer should not have a "category Q1" area
    And no errors should have been logged
