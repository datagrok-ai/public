@journey @viewers @realizes:viewers.line-chart
Feature: Line chart axis types, label orientation, line styling and a chart type per chart
  The settings of `line-chart.md` the other line chart features leave out: a logarithmic Y axis and
  back, the X axis labels turned vertical and back, the Lines section's width, transparency and
  colouring type, and — from `line-chart-ui.md` — a chart type given to one chart of a multi-axis
  chart through that chart's own menu, which leaves the other charts' types alone (in one multi-axis
  box the menu offers the submenu of the series under the pointer — HEIGHT at the plot's centre on
  this table — rather than the AGE the md right-clicks). One journey on
  demog-1000 with X = AGE and Y = WEIGHT; every scenario puts back what it changed.
  The pictures are not judged: a logarithmic axis and a colouring type are claimed by the repaint and
  the viewer's own readings, a vertical label strip by the X axis box growing taller, a line width by
  the ink it adds.
  The series' own "Chart type" is picked by the items' Dart names: its caption differs only in case
  from the chart-wide "Chart Type" group, which a caption path lands on instead.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a line chart viewer with:
      | xColumnName  | AGE    |
      | yColumnNames | WEIGHT |
    Then line chart viewer should report no error
    And the 'y axis max of "WEIGHT"' reading of line chart viewer should be a finite number

  Scenario: A logarithmic Y axis redraws the chart and a linear one draws it back
    When user sets "yAxisType" property of line chart viewer to "logarithmic"
    Then line chart viewer should have repainted
    And the 'y axis max of "WEIGHT"' reading of line chart viewer should be a finite number
    And line chart viewer should report no error
    When user sets "yAxisType" property of line chart viewer to "linear"
    Then line chart viewer should have repainted
    And line chart viewer should report no error
    And no errors should have been logged

  Scenario: Vertical X axis labels take a taller axis strip, Auto gives the room back
    When user sets "xAxisLabelOrientation" property of line chart viewer to "Vert"
    Then the "x axis" area of line chart viewer should be taller than before
    When user sets "xAxisLabelOrientation" property of line chart viewer to "Auto"
    Then the "x axis" area of line chart viewer should be shorter than before
    And no errors should have been logged

  Scenario: Line width, transparency and colouring type restyle the line
    When user sets "lineWidth" property of line chart viewer to "3"
    Then line chart viewer should have more ink than before
    When user sets "lineTransparency" property of line chart viewer to "0.5"
    Then line chart viewer should have repainted
    When user sets "lineColoringType" property of line chart viewer to "Custom"
    Then line chart viewer should have repainted
    And line chart viewer should report no error
    When user sets properties of line chart viewer:
      | lineWidth        | 1    |
      | lineTransparency | 0    |
      | lineColoringType | Auto |
    Then line chart viewer should have less ink than before
    And no errors should have been logged

  Scenario: In a multi-axis chart one chart's type changes through its own menu and the others keep theirs
    When user sets properties of line chart viewer:
      | yColumnNames | AGE, HEIGHT, WEIGHT |
      | multiAxis    | true                |
    Then the "charts" reading of line chart viewer should be 1
    And "chartTypes" property of line chart viewer should be "Line Chart, Line Chart, Line Chart"
    When user picks the item named "HEIGHT > Chart type > Area Chart" from the context menu of the "plot" area of line chart viewer
    Then "chartTypes" property of line chart viewer should be "Line Chart, Area Chart, Line Chart"
    And line chart viewer should have repainted
    When user sets properties of line chart viewer:
      | multiAxis    | false  |
      | yColumnNames | WEIGHT |
    Then the "y columns" reading of line chart viewer should be "WEIGHT"
    And no errors should have been logged
