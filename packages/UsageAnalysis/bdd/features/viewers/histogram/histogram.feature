@journey @viewers @realizes:viewers.histogram
Feature: Histogram property surface
  The histogram's property surface: the value column and the bin count from the panel and from the
  in-plot controls, the counts above the bins, spline and fill, the visibility of the on-plot
  chrome, the two axes and their scales, the Y range and the clipped-bin indicators, the legend a
  split brings, title and description, the context menu as a path to properties, and the Data panel
  through a layout round-trip on the server. Bins, selection and the range filter have features of
  their own. One journey on demog-1000 with a histogram of AGE; every scenario puts back what it
  changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a histogram viewer with:
      | Value | AGE |
    And user resizes histogram viewer to 500 by 400
    Then histogram viewer should show 1000 rows
    And the "bins shown" reading of histogram viewer should be 20
    And the "axis min" reading of histogram viewer should be 18
    And histogram viewer should have a "bin 8" area
    And histogram viewer should have an "x axis" area
    And histogram viewer should have a "y axis" area

  Scenario: The value column and the bin count, from the panel and from the plot
    When user sets "Value" property of histogram viewer to "WEIGHT"
    Then Value column input in histogram viewer should contain text "WEIGHT"
    And the "axis min" reading of histogram viewer should be higher than before
    And histogram viewer should have repainted by at least 2000 pixels
    When user takes a snapshot of histogram viewer
    And user selects "HEIGHT" in Value column input in histogram viewer
    Then "Value" property of histogram viewer should be "HEIGHT"
    And histogram viewer should have repainted by at least 2000 pixels
    When user sets "Value" property of histogram viewer to "AGE"
    Then the "axis min" reading of histogram viewer should be 18
    When user sets "Bins" property of histogram viewer to "50"
    Then the "bins shown" reading of histogram viewer should be 50
    And histogram viewer should have repainted by at least 2000 pixels
    When user moves the pointer away from histogram viewer
    And user hovers over histogram viewer
    Then histogram viewer should have a "bins slider" area
    When user drags the "bins slider" area of histogram viewer by 40 pixels to the right
    Then "Bins" property of histogram viewer should not be "50"
    And the "bins shown" reading of histogram viewer should be higher than before
    And histogram viewer should have repainted by at least 2000 pixels
    When user sets "Bins" property of histogram viewer to "20"
    Then the "bins shown" reading of histogram viewer should be 20
    And histogram viewer should have a "bin 8" area
    And no errors should have been logged

  Scenario: Show Values puts the counts above the bins
    Given histogram viewer should not have a "bin labels" area
    When user sets "Show Values" property of histogram viewer to "true"
    Then histogram viewer should have a "bin labels" area
    And the "bin labels" area of histogram viewer should be painted
    And histogram viewer should have repainted by at least 500 pixels
    When user sets "Show Values" property of histogram viewer to "false"
    Then histogram viewer should not have a "bin labels" area
    And histogram viewer should have repainted by at least 500 pixels
    And no errors should have been logged

  Scenario: Spline draws a line where the bars were
    When user sets "Spline" property of histogram viewer to "true"
    Then histogram viewer should not have a "bin 8" area
    And histogram viewer should have less ink than before
    And histogram viewer should have repainted by at least 2000 pixels
    When user sets "Fill Spline" property of histogram viewer to "true"
    Then histogram viewer should have more ink than before
    When user sets "Fill Spline" property of histogram viewer to "false"
    Then histogram viewer should have less ink than before
    When user sets "Spline" property of histogram viewer to "false"
    Then histogram viewer should have a "bin 8" area
    And histogram viewer should have more ink than before
    And no errors should have been logged

  Scenario: Controls visibility
    When user moves the pointer away from histogram viewer
    And user hovers over histogram viewer
    Then Value column input in histogram viewer should be visible
    And histogram viewer should have a "bins slider" area
    And histogram viewer should have a "range min handle" area
    And histogram viewer should have a "range max handle" area
    When user sets properties of histogram viewer:
      | Show Column Selector | false |
      | Show Bin Selector    | false |
      | Show Range Slider    | false |
    And user moves the pointer away from histogram viewer
    And user hovers over histogram viewer
    Then Value column input in histogram viewer should be hidden
    And histogram viewer should not have a "bins slider" area
    And histogram viewer should not have a "range min handle" area
    When user sets properties of histogram viewer:
      | Show Column Selector | true |
      | Show Bin Selector    | true |
      | Show Range Slider    | true |
    And user moves the pointer away from histogram viewer
    And user hovers over histogram viewer
    Then Value column input in histogram viewer should be visible
    And histogram viewer should have a "bins slider" area
    When user sets "Allow Column Selection" property of histogram viewer to "false"
    Then Value column input in histogram viewer should be disabled
    When user sets "Allow Column Selection" property of histogram viewer to "true"
    Then Value column input in histogram viewer should be enabled
    And no errors should have been logged

  Scenario: The two axes and their scales
    When user sets "Show X Axis" property of histogram viewer to "false"
    Then histogram viewer should not have an "x axis" area
    And histogram viewer should have repainted by at least 500 pixels
    When user sets "Show Y Axis" property of histogram viewer to "false"
    Then histogram viewer should not have a "y axis" area
    And histogram viewer should have repainted by at least 500 pixels
    When user sets properties of histogram viewer:
      | Show X Axis | true |
      | Show Y Axis | true |
    Then histogram viewer should have an "x axis" area
    And histogram viewer should have a "y axis" area
    When user sets "X Axis Type" property of histogram viewer to "logarithmic"
    Then histogram viewer should have repainted by at least 500 pixels
    And histogram viewer should have a "bin 8" area
    When user sets "X Axis Type" property of histogram viewer to "linear"
    Then the "axis min" reading of histogram viewer should be 18
    When user sets "Y Axis Type" property of histogram viewer to "logarithmic"
    Then histogram viewer should have repainted by at least 500 pixels
    When user sets "Y Axis Type" property of histogram viewer to "linear"
    Then histogram viewer should have repainted by at least 500 pixels
    And no errors should have been logged

  Scenario: The Y range clips the tall bins and says how many
    Given the "y axis max" reading of histogram viewer should be 99
    And the "clipped bins" reading of histogram viewer should be 0
    When user sets "Y Max" property of histogram viewer to "50"
    Then the "y axis max" reading of histogram viewer should be 50
    And the "clipped bins" reading of histogram viewer should be 10
    And histogram viewer should have repainted by at least 500 pixels
    When user sets "Show Clipped Bin Indicators" property of histogram viewer to "false"
    Then the "clipped bins" reading of histogram viewer should be 0
    And histogram viewer should have repainted
    When user sets properties of histogram viewer:
      | Show Clipped Bin Indicators | true |
      | Y Max                       |      |
    Then the "y axis max" reading of histogram viewer should be 99
    And the "clipped bins" reading of histogram viewer should be 0
    And no errors should have been logged

  Scenario: A split brings a legend
    When user sets properties of histogram viewer:
      | Split             | SEX    |
      | Legend Visibility | Always |
    Then legend of histogram viewer should be visible
    And the legend of histogram viewer should list 2 items
    And legend of histogram viewer should contain text "F"
    And legend of histogram viewer should contain text "M"
    When user sets "Split" property of histogram viewer to "RACE"
    Then the legend of histogram viewer should list 4 items
    When user sets "Legend Position" property of histogram viewer to "Right"
    Then the legend of histogram viewer should be on the right
    When user sets "Legend Position" property of histogram viewer to "Left"
    Then the legend of histogram viewer should be on the left
    When user sets "Legend Visibility" property of histogram viewer to "Never"
    Then legend of histogram viewer should be hidden
    When user sets properties of histogram viewer:
      | Legend Visibility | Auto |
      | Legend Position   | Auto |
      | Split             |      |
    Then legend of histogram viewer should be hidden
    And histogram viewer should have a "bin 8" area
    And no errors should have been logged

  Scenario: The context menu as a path to properties
    When user opens the context menu of histogram viewer
    Then "Show Filtered Out Rows" menu item in context menu should be visible
    And "Selection" menu item in context menu should be visible
    And "Tools" menu item in context menu should be visible
    When user hovers over "Selection" menu item in context menu
    Then "Show Current Row" menu item in context menu should be visible
    And "Show Mouse Over Row" menu item in context menu should be visible
    And "Show Mouse Over Row Group" menu item in context menu should be visible
    When user closes the context menu
    And user picks "Selection > Show Current Row" from the context menu of histogram viewer
    Then "Show Current Row" property of histogram viewer should be "false"
    When user picks "Selection > Show Current Row" from the context menu of histogram viewer
    Then "Show Current Row" property of histogram viewer should be "true"
    When user right-clicks on the "x axis" area of histogram viewer
    Then "Show X Axis" menu item in context menu should be visible
    And "X Axis Type" menu item in context menu should be visible
    When user closes the context menu
    And user picks "Show X Axis" from the context menu of the "x axis" area of histogram viewer
    Then "Show X Axis" property of histogram viewer should be "false"
    And histogram viewer should not have an "x axis" area
    And histogram viewer should have repainted by at least 500 pixels
    When user sets "Show X Axis" property of histogram viewer to "true"
    Then histogram viewer should have an "x axis" area
    And no errors should have been logged

  Scenario: Show Filtered Out Rows drops the totals behind the filtered bins
    When user filters rows where "SEX" is "F"
    Then 553 rows should pass the filter
    And histogram viewer should show 553 rows
    And the "y axis max" reading of histogram viewer should be 99
    And the "bin 8" area of histogram viewer should contain the color "#F0F0F0"
    When user sets "Show Filtered Out Rows" property of histogram viewer to "false"
    Then the "y axis max" reading of histogram viewer should be 53
    And the "bin 8" area of histogram viewer should not contain the color "#F0F0F0"
    And histogram viewer should have repainted by at least 2000 pixels
    When user sets "Show Filtered Out Rows" property of histogram viewer to "true"
    Then the "y axis max" reading of histogram viewer should be 99
    And the "bin 8" area of histogram viewer should contain the color "#F0F0F0"
    When user resets the filter
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: The Data panel through a layout round-trip on the server
    Given user opens spgi dataset
    And user switches to the "demog-1000" table view
    When user sets "Table" property of histogram viewer to "spgi-100"
    Then histogram viewer should be bound to table "spgi-100"
    And "Value" property of histogram viewer should be "CAST Idea ID"
    And histogram viewer should show 100 rows
    When user sets "Filter" property of histogram viewer to "${CAST Idea ID} < 634835"
    Then histogram viewer should show 50 rows
    And 100 rows of table "spgi-100" should pass the filter
    When user saves the layout of the current table view to the server
    And user clicks on close icon of histogram viewer
    Then histogram viewer should be absent
    When user loads the saved layout
    Then histogram viewer should be visible
    And properties of histogram viewer should be:
      | Table  | spgi-100                 |
      | Filter | ${CAST Idea ID} < 634835 |
      | Value  | CAST Idea ID             |
    And histogram viewer should show 50 rows
    And no errors should have been logged
    When user sets properties of histogram viewer:
      | Filter |            |
      | Table  | demog-1000 |
    And user sets "Value" property of histogram viewer to "AGE"
    Then histogram viewer should be bound to table "demog-1000"
    And histogram viewer should show 1000 rows
    And the "axis min" reading of histogram viewer should be 18
