@viewers @realizes:viewers.heat-map
Feature: Heat map regressions of closed issues
  Closed heat map issues, each claimed by what a person sees. core#928: with Max Heatmap Columns 0
  the map is empty, and it now says why inside the viewer ("No columns to show: Max Heatmap Columns
  is 0"); 100 brings the columns back and takes the hint away. The hint is read as text shown in the
  viewer: the grid's status does not publish it as its `error` (see MISSING-regressions.md). core#931: a heat map
  stacked with a histogram as tabs comes back from a layout on the tab that was active at the save,
  the Heat map tab, with no switch; the histogram of the other tab is not on the page. The layout is
  deleted at feature end and read back gone.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a heat map viewer
    Then the "columns shown" reading of heat map viewer should be 11
    And heat map viewer should be painted

  @realizes:core#928
  Scenario: Max Heatmap Columns 0 says why the map is empty, and 100 takes the hint away (core#928)
    When user sets "maxHeatmapColumns" property of heat map viewer to "0"
    Then heat map viewer should show the text "No columns to show: Max Heatmap Columns is 0"
    And the "columns shown" reading of heat map viewer should be 0
    And heat map viewer should not have a "column AGE" area
    When user sets "maxHeatmapColumns" property of heat map viewer to "100"
    Then the "columns shown" reading of heat map viewer should be 11
    And heat map viewer should have a "column AGE" area
    And heat map viewer should not show the text "No columns to show"
    And no errors should have been logged

  @realizes:core#931
  Scenario: A layout brings the stacked heat map back on its own tab (core#931)
    Given user adds a histogram viewer
    When user stacks heat map viewer onto histogram viewer as a tab
    And user switches the tabbed panel to the "Histogram" tab
    And user switches the tabbed panel to the "Heat map" tab
    Then heat map viewer should have a "column AGE" area
    When user saves the layout of the current table view to the server
    And user closes all views
    And user opens demog-1000 dataset
    And user loads the saved layout
    Then the viewers "Histogram, Heat map" should be the tabs of one panel
    And heat map viewer should have a "column AGE" area
    And heat map viewer should be painted
    And histogram viewer should be absent
    And no errors should have been logged
