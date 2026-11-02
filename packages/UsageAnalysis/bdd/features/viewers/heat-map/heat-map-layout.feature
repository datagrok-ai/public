@viewers @realizes:viewers.heat-map
Feature: Heat map layouts that used to load broken
  Regressions of what a person saw after loading a layout, not of what a layout stores (that is an
  API test). public#2927: a heat map stacked with another viewer as a tab and zoomed on its
  horizontal slider came back collapsed to a single column; it now comes back with as many
  columns on screen as it had. public#2901: a layout saved after Is Heatmap was turned off and on
  again broke on Max Heatmap Columns; it loads as a heat map with every column, and one saved in
  grid mode with Max Heatmap Columns 200 loads as a grid that turns back into the whole heat map.
  GROK-16239: a heat map layout saved on the thousand rows of demog-1000 broke the map when loaded
  on a 100-row table; it now loads with all 100 rows shown, down to the last row's cells.
  Each scenario saves the layout through the server, closes the table and opens it again before
  loading, as the reports did; the layout is deleted at feature end and read back gone, and what
  a killed run left is swept on the first save (`user saves the layout of the current table view
  to the server`). The scroll window a layout restores (public#2926) is not claimed on its own:
  it is what the layout stores, and public#2927's claim covers what a person sees of it.
  core#931: the layout brings the stacked heat map back on the tab that was active at the save, the
  Heat map tab, with no switch, and the histogram of the other tab is not on the page. The grid's
  status publishes no viewer error, so no scenario claims "should report no error" here (MISSING.md).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a heat map viewer
    Then the "is heatmap" reading of heat map viewer should be "true"
    And the "columns shown" reading of heat map viewer should be 11
    And heat map viewer should be painted

  @realizes:core#931
  Scenario: A zoomed heat map stacked as a tab is not collapsed after the layout loads, on the tab it was saved on (public#2927, core#931)
    Given user adds a histogram viewer
    When user stacks heat map viewer onto histogram viewer as a tab
    And user switches the tabbed panel to the "Heat map" tab
    And user drags the "x scroll max handle" area of heat map viewer by 200 pixels to the left
    Then the "x scroll span" reading of heat map viewer should be between 0.4 and 0.9
    And the "columns shown" reading of heat map viewer should be between 4 and 10
    When user remembers the "columns shown" reading of heat map viewer
    And user saves the layout of the current table view to the server
    And user closes all views
    And user opens demog-1000 dataset
    And user loads the saved layout
    Then the viewers "Histogram, Heat map" should be the tabs of one panel
    And histogram viewer should be absent
    And the "columns shown" reading of heat map viewer should be as remembered
    And the "x scroll span" reading of heat map viewer should be between 0.4 and 0.9
    And heat map viewer should have a "column AGE" area
    And heat map viewer should be painted
    And no errors should have been logged

  Scenario: A layout saved after Is Heatmap off and on loads as the whole heat map (public#2901)
    When user sets "isHeatmap" property of heat map viewer to "false"
    Then the "row height" reading of heat map viewer should be between 20 and 40
    When user sets "isHeatmap" property of heat map viewer to "true"
    And user saves the layout of the current table view to the server
    And user closes all views
    And user opens demog-1000 dataset
    And user loads the saved layout
    Then the "is heatmap" reading of heat map viewer should be "true"
    And the "columns shown" reading of heat map viewer should be 11
    And the "row height" reading of heat map viewer should be between 0 and 8
    And the "y scroll span" reading of heat map viewer should be 1
    And heat map viewer should be painted
    And no errors should have been logged

  Scenario: A layout saved in grid mode with Max Heatmap Columns 200 turns back into the heat map (public#2901)
    When user sets "maxHeatmapColumns" property of heat map viewer to "200"
    And user sets "isHeatmap" property of heat map viewer to "false"
    And user saves the layout of the current table view to the server
    And user closes all views
    And user opens demog-1000 dataset
    And user loads the saved layout
    Then the "is heatmap" reading of heat map viewer should be "false"
    And the "row height" reading of heat map viewer should be between 20 and 40
    And heat map viewer should have a "cell 1 of AGE" area
    When user sets "isHeatmap" property of heat map viewer to "true"
    Then the "is heatmap" reading of heat map viewer should be "true"
    And the "columns shown" reading of heat map viewer should be 11
    And the "row height" reading of heat map viewer should be between 0 and 8
    And the "y scroll span" reading of heat map viewer should be 1
    And heat map viewer should have a "column AGE" area
    And no errors should have been logged

  Scenario: A layout saved on 1000 rows loads on a 100-row table with all 100 rows shown (GROK-16239)
    When user saves the layout of the current table view to the server
    And user closes all views
    And user opens demog-1000 dataset keeping the first 100 rows
    And user loads the saved layout
    Then the "is heatmap" reading of heat map viewer should be "true"
    And the "rows" reading of heat map viewer should be 100
    And the "rows shown" reading of heat map viewer should be 100
    And the "columns shown" reading of heat map viewer should be 11
    And the "y scroll span" reading of heat map viewer should be 1
    And heat map viewer should have a "cell 1 of AGE" area
    And heat map viewer should have a "cell 100 of AGE" area
    And heat map viewer should be painted
    And no errors should have been logged
