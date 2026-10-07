@viewers @realizes:viewers.trellis-plot
Feature: Use in Trellis from another viewer
  General > Use in Trellis in a viewer's context menu adds a trellis plot whose inner viewer is that
  viewer, with the viewer's own settings carried into every cell, and leaves the source viewer in the
  view. Translated from the "Use in Trellis" section of TestTrack Viewers/TrellisPlot/trellis-plot.md
  and the "Create a trellis from another viewer" item of trellis-plot-ui.md, for the five viewers the
  case walks: scatter plot, bar chart, histogram, line chart and box plot, each configured first and
  its setting read back from the trellis's inner viewer. Each scenario starts on a
  fresh demog-1000 view (the pie chart path is walked by the Embedded Viewers tutorial).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  Scenario: A scatter plot in a trellis keeps its X, Y and color columns
    Given user adds a scatter plot viewer with:
      | X     | AGE    |
      | Y     | HEIGHT |
      | Color | SEX    |
    When user picks "General > Use in Trellis" from the context menu of scatter plot viewer
    Then the open tableview should have 1 trellis plot viewer
    And the open tableview should have 1 scatter plot viewer
    And the "inner viewer type" reading of trellis plot viewer should be "Scatter plot"
    And "xColumnName" inner property of trellis plot viewer should be "AGE"
    And "yColumnName" inner property of trellis plot viewer should be "HEIGHT"
    And "colorColumnName" inner property of trellis plot viewer should be "SEX"
    And the "cells drawn" reading of trellis plot viewer should be 23
    And trellis plot viewer should report no error
    And no errors should have been logged

  Scenario Outline: A <viewer> in a trellis becomes its inner viewer, its <caption> carried over
    Given user adds a <viewer> viewer
    And user sets "<caption>" property of <viewer> viewer to "<value>"
    Then <viewer> viewer should be painted
    When user picks "General > Use in Trellis" from the context menu of the "<area>" area of <viewer> viewer
    Then the open tableview should have 1 trellis plot viewer
    And the open tableview should have 1 <viewer> viewer
    And the "inner viewer type" reading of trellis plot viewer should be "<type>"
    And "<inner>" inner property of trellis plot viewer should be "<value>"
    And the "cells drawn" reading of trellis plot viewer should be <cells>
    And trellis plot viewer should report no error
    And no errors should have been logged

    Examples:
      | viewer     | type       | area | caption | value   | inner           | cells |
      | bar chart  | Bar chart  | view | Value   | WEIGHT  | valueColumnName | 28    |
      | histogram  | Histogram  | view | Value   | HEIGHT  | valueColumnName | 23    |
      | line chart | Line chart | plot | X       | STARTED | xColumnName     | 28    |
      | box plot   | Box plot   | view | Value   | WEIGHT  | valueColumnName | 28    |
