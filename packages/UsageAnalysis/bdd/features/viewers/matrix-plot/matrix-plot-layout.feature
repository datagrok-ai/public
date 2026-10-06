@journey @viewers @realizes:viewers.matrix-plot
Feature: Matrix plot — the chrome around the grid and the layout that decides it
  The column labels, the two axes strips, the title and the description, and Auto Layout, which
  drops a strip that no longer fits whatever the axis flag asks for.
  Every claim here is that the plot stops REPORTING a strip, not that a canvas somewhere lost its
  width: `x axis shown` is the computed `showXAxes` getter (matrix_plot_core.dart:311), the same one
  the layout uses, and the `x axis` region is reported only while that getter and the element's own
  display agree — where the spec this replaces polled `.d4-layout-top canvas`.width until it hit 0.
  Not translated: GROK-18736, the Font property reaching the labels. `label font` is the look read
  back, and nothing else the plot reports moves with it — the label strip is a fixed 20 px and every
  label box, every cell box and every cell signature is identical at 10 px and at 24 px (measured).
  The old spec read `getComputedStyle(div).font` off whichever leaf div held the text "AGE"; an
  honest successor needs the strip's effective label height, or the labels' computed font, in the
  status.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a matrix plot viewer
    Then the "cells" reading of matrix plot viewer should be 16
    And the "auto layout" reading of matrix plot viewer should be "true"
    And the "x axis shown" reading of matrix plot viewer should be "true"
    And the "y axis shown" reading of matrix plot viewer should be "true"

  Scenario: Show X Axes and Show Y Axes take their strips away and bring them back
    Then matrix plot viewer should have an "x axis" area
    And matrix plot viewer should have a "y axis" area
    And matrix plot viewer should have an "x axis cell 1" area
    When user sets "showXAxes" property of matrix plot viewer to "false"
    Then the "x axis shown" reading of matrix plot viewer should be "false"
    And matrix plot viewer should not have an "x axis" area
    And matrix plot viewer should not have an "x axis cell 1" area
    And the "y axis shown" reading of matrix plot viewer should be "true"
    And matrix plot viewer should have a "y axis" area
    When user sets "showYAxes" property of matrix plot viewer to "false"
    Then the "y axis shown" reading of matrix plot viewer should be "false"
    And matrix plot viewer should not have a "y axis" area
    When user sets properties of matrix plot viewer:
      | showXAxes | true |
      | showYAxes | true |
    Then the "x axis shown" reading of matrix plot viewer should be "true"
    And matrix plot viewer should have an "x axis" area
    And matrix plot viewer should have a "y axis" area
    And no errors should have been logged

  Scenario: Auto Layout drops the strips that no longer fit, and turning it off puts them back (GROK-19106)
    Then the "x axis shown" reading of matrix plot viewer should be "true"
    And the "x labels shown" reading of matrix plot viewer should be "true"
    When user resizes matrix plot viewer to 260 by 240
    Then "showXAxes" property of matrix plot viewer should be "true"
    And "showYAxes" property of matrix plot viewer should be "true"
    And the "x axis shown" reading of matrix plot viewer should be "false"
    And the "y axis shown" reading of matrix plot viewer should be "false"
    And the "x labels shown" reading of matrix plot viewer should be "false"
    And the "y labels shown" reading of matrix plot viewer should be "false"
    And matrix plot viewer should not have an "x axis" area
    And matrix plot viewer should not have an "x label AGE" area
    When user sets "autoLayout" property of matrix plot viewer to "false"
    Then the "x axis shown" reading of matrix plot viewer should be "true"
    And the "y axis shown" reading of matrix plot viewer should be "true"
    And the "x labels shown" reading of matrix plot viewer should be "true"
    And matrix plot viewer should have an "x axis" area
    And matrix plot viewer should have an "x label AGE" area
    When user sets "autoLayout" property of matrix plot viewer to "true"
    And user restores the size of matrix plot viewer
    Then the "x axis shown" reading of matrix plot viewer should be "true"
    And the "x labels shown" reading of matrix plot viewer should be "true"
    And no errors should have been logged

  Scenario: The label strips carry one box per column of the viewport
    Then the "x labels shown" reading of matrix plot viewer should be "true"
    And matrix plot viewer should have an "x labels" area
    And matrix plot viewer should have an "x label AGE" area
    And matrix plot viewer should have an "x label STARTED" area
    And matrix plot viewer should have a "y labels" area
    And matrix plot viewer should have a "y label WEIGHT" area
    And matrix plot viewer should not have an "x label SEX" area
    And no errors should have been logged

  Scenario: The title shows and the description sits above the grid, then below it
    When user sets properties of matrix plot viewer:
      | showTitle   | true                 |
      | title       | Pairwise             |
      | description | Every numerical pair |
    Then title of matrix plot viewer should have text "Pairwise"
    And description of matrix plot viewer should be visible
    And description of matrix plot viewer should have text "Every numerical pair"
    And the description of matrix plot viewer should be above its content
    When user sets "descriptionPosition" property of matrix plot viewer to "Bottom"
    Then the description of matrix plot viewer should be below its content
    When user sets "descriptionVisibilityMode" property of matrix plot viewer to "Never"
    Then description of matrix plot viewer should be hidden
    When user sets properties of matrix plot viewer:
      | descriptionVisibilityMode | Auto  |
      | descriptionPosition       | Top   |
      | description               |       |
      | title                     |       |
      | showTitle                 | false |
    Then the "cells" reading of matrix plot viewer should be 16
    And no errors should have been logged

  Scenario: A saved layout brings the configured grid back
    When user sets properties of matrix plot viewer:
      | cellPlotType | Scatter plot        |
      | xColumnNames | AGE, HEIGHT         |
      | yColumnNames | AGE, HEIGHT, WEIGHT |
    Then the "cells" reading of matrix plot viewer should be 6
    And the cells of matrix plot viewer should be 2 wide and 3 tall
    And user saves the layout of the current table view
    When user sets properties of matrix plot viewer:
      | cellPlotType | Density plot                 |
      | xColumnNames | AGE, HEIGHT, WEIGHT, STARTED |
      | yColumnNames | AGE, HEIGHT, WEIGHT, STARTED |
    Then the "cells" reading of matrix plot viewer should be 16
    When user loads the saved layout
    Then matrix plot viewer should be visible
    And the "cells" reading of matrix plot viewer should be 6
    And the cells of matrix plot viewer should be 2 wide and 3 tall
    And the "cell viewer type" reading of matrix plot viewer should be "Scatter plot"
    And the "cell viewer type of HEIGHT x AGE" reading of matrix plot viewer should be "Scatter plot"
    And the "cell viewer type of AGE x AGE" reading of matrix plot viewer should be "Histogram"
    When user sets properties of matrix plot viewer:
      | cellPlotType | Density plot                 |
      | xColumnNames | AGE, HEIGHT, WEIGHT, STARTED |
      | yColumnNames | AGE, HEIGHT, WEIGHT, STARTED |
    Then the "cells" reading of matrix plot viewer should be 16
    And no errors should have been logged

  Scenario: A project round-trip brings it back too (GROK-10925)
    When user sets properties of matrix plot viewer:
      | cellPlotType | Scatter plot        |
      | xColumnNames | AGE, HEIGHT         |
      | yColumnNames | AGE, HEIGHT, WEIGHT |
    Then the "cells" reading of matrix plot viewer should be 6
    When user saves the current view as project "bdd matrix plot grid"
    And user closes all views
    And user opens the "bdd matrix plot grid" project
    Then matrix plot viewer should be visible
    And the "cells" reading of matrix plot viewer should be 6
    And the cells of matrix plot viewer should be 2 wide and 3 tall
    And the "cell viewer type" reading of matrix plot viewer should be "Scatter plot"
    And the "rows shown" reading of matrix plot viewer should be 1000
    And the "cell rows shown of HEIGHT x AGE" reading of matrix plot viewer should be 872
    When user sets properties of matrix plot viewer:
      | cellPlotType | Density plot                 |
      | xColumnNames | AGE, HEIGHT, WEIGHT, STARTED |
      | yColumnNames | AGE, HEIGHT, WEIGHT, STARTED |
    Then the "cells" reading of matrix plot viewer should be 16
    And no errors should have been logged
