@viewers @realizes:viewers.heat-map
Feature: Heat map known failures
  Temporary: the known failures of the heat map while their tickets are open. Each scenario is
  expected to fail on its pinned step; when its ticket is fixed, remove the tag and move the scenario
  back into heat-map.feature as an ordinary @edge-cases scenario.

  Each scenario is a test of its own on demog-1000 with a heat map of every column.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a heat map viewer
    Then the "is heatmap" reading of heat map viewer should be "true"
    And the "rows shown" reading of heat map viewer should be 1000

  @edge-cases @known-failure:32 @realizes:core#1202
  Scenario: A column hidden under Max Heatmap Columns stays hidden when the cap is raised (core#1202)
    When user sets "maxHeatmapColumns" property of heat map viewer to "5"
    Then the "columns shown" reading of heat map viewer should be 5
    And heat map viewer should have a "column AGE" area
    When user clicks on settings icon of heat map viewer
    Given "Columns" category in context panel is expanded
    When user clicks on first button in "Columns" property in context panel
    Then Order or Hide Columns dialog should be visible
    When user toggles the "AGE" column in the column list of Order or Hide Columns dialog
    And user clicks on CLOSE button in Order or Hide Columns dialog
    Then the "columns shown" reading of heat map viewer should be 4
    And heat map viewer should not have a "column AGE" area
    When user remembers the "columns shown" reading of heat map viewer
    And user enters "6" into "Max Heatmap Columns" property in context panel
    Then the "columns shown" reading of heat map viewer should be higher than remembered
    And heat map viewer should not have a "column AGE" area
