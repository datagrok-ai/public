@viewers @edge-cases @realizes:viewers.pc-plot
Feature: PC plot edge cases — a layout with line width 0, an unparsable font, one column, an empty source
  The PC plot on settings and data its other features do not reach. A layout saved with Line Width
  0 and density on, applied back to the view, renders and takes a hover without an error, and an
  Axis Font the browser cannot parse (`bold`) leaves the plot painting (both crashes of core#967,
  fixed). With a single column the plot reports "please select at least 2 columns" and a second
  column brings it back (`validateState`, `pc_plot_core.dart`). A table filtered down to nothing
  and a table whose rows are all deleted leave a plot that draws no line and logs no error;
  `lines drawn` and `rows shown` at 0 tell an empty source, not the ink, since the axes stay drawn.
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

  Scenario: With every row filtered out of the table the plot draws no line and logs no error
    When user filters rows where "AGE" is between 200 and 300
    Then 0 rows should pass the filter
    And pc plot viewer should show 0 rows
    And the "lines drawn" reading of pc plot viewer should be 0
    And no errors should have been logged
    When user resets the filter
    Then the "lines drawn" reading of pc plot viewer should be 1000

  Scenario: Once every row of the table is deleted the plot draws no line and logs no error
    When user selects all rows
    And user deletes the selected rows
    Then the table should have 0 rows
    And pc plot viewer should show 0 rows
    And the "lines drawn" reading of pc plot viewer should be 0
    And no errors should have been logged
