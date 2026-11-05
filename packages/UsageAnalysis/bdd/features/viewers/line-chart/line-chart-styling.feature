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
  Closed issues on the same chart: the split legend lists no Y / split pair without data under Multi
  Axis and Y Global Scale (core#825, public#3463 — the calculated H_F holds HEIGHT for the women, W_M
  WEIGHT for the men); under Row Source = Selected the split legend follows the selected rows
  (public#1651); a split with a category per row (USUBJID) keeps the chart working (public#3210); a row
  group hovered in a bar chart repaints the chart only while Show Mouse Over Row Group is on
  (public#3709).

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
    When user picks "HEIGHT > Chart type > Area Chart" from the context menu of the "plot" area of line chart viewer
    And user closes the context menu
    Then "chartTypes" property of line chart viewer should be "Line Chart, Area Chart, Line Chart"
    And line chart viewer should have repainted
    When user sets properties of line chart viewer:
      | multiAxis    | false  |
      | yColumnNames | WEIGHT |
    Then the "y columns" reading of line chart viewer should be "WEIGHT"
    And no errors should have been logged

  @realizes:core#825 @realizes:public#3463
  Scenario: The split legend lists no Y / split pair that has no data (core#825)
    When user adds a calculated column "H_F" with formula "if(${SEX} == \"F\", ${HEIGHT}, null)"
    And user adds a calculated column "W_M" with formula "if(${SEX} == \"M\", ${WEIGHT}, null)"
    And user sets properties of line chart viewer:
      | yColumnNames     | H_F, W_M |
      | splitColumnNames | SEX      |
      | multiAxis        | true     |
      | yGlobalScale     | true     |
    Then the legend of line chart viewer should list 2 items
    And "H_F / F" legend item in legend of line chart viewer should be visible
    And "W_M / M" legend item in legend of line chart viewer should be visible
    And "H_F / M" legend item in legend of line chart viewer should be absent
    And "W_M / F" legend item in legend of line chart viewer should be absent
    When user sets properties of line chart viewer:
      | yColumnNames     | WEIGHT |
      | splitColumnNames |        |
      | multiAxis        | false  |
      | yGlobalScale     | false  |
    And user removes "H_F" column
    And user removes "W_M" column
    Then no errors should have been logged

  @realizes:public#1651
  Scenario: Under Row Source = Selected the split legend lists the categories of the selected rows (public#1651)
    When user sets properties of line chart viewer:
      | splitColumnNames | RACE     |
      | rowSource        | Selected |
    And user selects rows where "RACE" is "Asian"
    Then the legend of line chart viewer should list 1 item
    And "Asian" legend item in legend of line chart viewer should be visible
    When user selects rows where "RACE" is one of "Asian, Black"
    Then the legend of line chart viewer should list 2 items
    And "Black" legend item in legend of line chart viewer should be visible
    And "Caucasian" legend item in legend of line chart viewer should be absent
    When user clears the row selection
    And user sets properties of line chart viewer:
      | splitColumnNames |          |
      | rowSource        | Filtered |
    Then no errors should have been logged

  @realizes:public#3210
  Scenario: A split with a category per row keeps the chart working (public#3210)
    When user sets "splitColumnNames" property of line chart viewer to "USUBJID"
    Then the "categories" reading of line chart viewer should be 1000
    And line chart viewer should report no error
    When user sets "splitColumnNames" property of line chart viewer to ""
    Then the "categories" reading of line chart viewer should be 1
    And line chart viewer should be painted
    And no errors should have been logged

  @realizes:public#3709
  Scenario: A row group hovered in a bar chart repaints the line chart only while Show Mouse Over Row Group is on (public#3709)
    When user adds a bar chart viewer with:
      | Split | RACE |
    And user takes a snapshot of line chart viewer
    And user hovers over the "bar Asian" area of bar chart viewer
    Then line chart viewer should have repainted
    When user moves the pointer away from grid
    And user sets "showMouseOverRowGroup" property of line chart viewer to "false"
    And user takes a snapshot of line chart viewer
    And user hovers over the "bar Asian" area of bar chart viewer
    Then bar chart viewer should have repainted
    And line chart viewer should not have repainted
    When user moves the pointer away from grid
    And user sets "showMouseOverRowGroup" property of line chart viewer to "true"
    And user clicks on close icon of bar chart viewer
    Then no errors should have been logged
