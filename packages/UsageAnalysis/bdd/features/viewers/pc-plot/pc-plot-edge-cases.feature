@viewers @edge-cases @realizes:viewers.pc-plot
Feature: PC plot edge cases — a layout with line width 0, an unparsable font, one column, an empty source
  The PC plot on settings and data its other features do not reach. A layout saved with Line Width
  0 and density on, applied back to the view, renders and takes a hover without an error, and an
  Axis Font the browser cannot parse (`bold`) leaves the plot painting (both crashes of core#967,
  fixed). With a single column the plot reports "please select at least 2 columns" and a second
  column brings it back (`validateState`, `pc_plot_core.dart`). Filtered down to nothing or
  with every row deleted: `viewer-rules-edge-cases.feature`.
  demog-1000, AGE, HEIGHT and WEIGHT.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pc plot viewer with:
      | Column Names | AGE, HEIGHT, WEIGHT |
    Then the axes of pc plot viewer should be "AGE, HEIGHT, WEIGHT"
    And pc plot viewer should show 1000 rows
    And the "error" reading of pc plot viewer should be ""

  Scenario: A layout saved with Line Width 0 and density on is applied and hovered without an error (core#967)
    When user sets properties of pc plot viewer:
      | Show Density | true |
      | Line Width   | 0    |
    And user saves the layout of the current table view as "pc-width-0"
    And user applies the layout "pc-width-0" to the current table view
    Then "Line Width" property of pc plot viewer should be "0"
    And the "density style" reading of pc plot viewer should be "circles"
    And pc plot viewer should be painted
    When user hovers over the "view" area of pc plot viewer
    And user moves the pointer away from pc plot viewer
    Then the "error" reading of pc plot viewer should be ""
    And no errors should have been logged

  Scenario: An Axis Font the browser cannot parse leaves the plot painting without an error (core#967)
    When user sets "Axis Font" property of pc plot viewer to "bold"
    And user takes a snapshot of pc plot viewer
    And user hovers over the "view" area of pc plot viewer
    Then pc plot viewer should be painted
    And the "lines drawn" reading of pc plot viewer should be 1000
    And no errors should have been logged

  Scenario: With a single column the plot asks for at least two, and a second column brings it back
    When user sets "Column Names" property of pc plot viewer to "AGE"
    Then the "error" reading of pc plot viewer should be "please select at least 2 columns"
    And no errors should have been logged
    When user sets "Column Names" property of pc plot viewer to "AGE, WEIGHT"
    Then the "error" reading of pc plot viewer should be ""
    And the axes of pc plot viewer should be "AGE, WEIGHT"
    And the "lines drawn" reading of pc plot viewer should be 1000
    And no errors should have been logged

  @realizes:core#967
  Scenario: A transformation turned on releases the filter of the axis sliders it replaces (core#967)
    When user hovers over pc plot viewer
    And user drags the max handle of the "AGE" range slider of pc plot viewer by 120 pixels
    Then fewer than 1000 rows should pass the filter
    When user sets "Transformation" property of pc plot viewer to '[{"#type":"GroupAggregation","aggType":"key","colName":"SEX"},{"#type":"GroupAggregation","aggType":"avg","colName":"AGE"},{"#type":"GroupAggregation","aggType":"avg","colName":"WEIGHT"}]'
    Then the "transformed" reading of pc plot viewer should be "true"
    And all rows should pass the filter
    And no errors should have been logged

  @known-failure:67 @realizes:core#1203
  Scenario: Log Columns on a column holding zero says the axis cannot be logarithmic (core#1203)
    When user sets "WEIGHT" column in row 1 to "0"
    And user sets "Log Columns" property of pc plot viewer to "WEIGHT"
    Then "Log Columns" property of pc plot viewer should be "WEIGHT"
    And the "error" reading of pc plot viewer should not be ""
    And no errors should have been logged
