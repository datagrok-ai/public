@journey @viewers @realizes:viewers.box-plot
Feature: Box plot — hovering a category label highlights its row group
  Hovering a category's label on the category axis highlights that row group, and Show Mouse Over
  Row Group switches it off (box-plot-ui "Mouse-over highlight"). A hover on another viewer's bar is
  in `box-plot-selection.feature` and the date category mapping in `box-plot.feature`; that the
  other groups fade has no claim yet (MISSING.md). On demog-1000 with a box plot of AGE by RACE and
  10-pixel markers; the scenario puts back what it changed.

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
