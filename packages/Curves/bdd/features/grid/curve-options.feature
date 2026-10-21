@journey
Feature: Curve options and the level they are set at
  Translated from TestTrack Curves/curves-options.md, scenarios 1 and 2. An option set at Column level
  in the Options pane of one curve cell is what the other cells of that column show, and another column
  keeps its own; an option set at Cell level stays in its cell (GROK-12239). The option is read where a
  person reads it, in the Options pane of another cell: since #4306 a click on a second curve cell keeps
  the curve panes in the Context Panel.

  In the package's tests instead: a Column level option kept by a layout applied to the table opened
  again (panel-renderer-tests "a column-level option survives a layout applied to a fresh table") and by
  a project (project-tests "column-level option survives a project"), scenario 6. Not claimed: what an
  option changes on the canvas - connect dots, markers, the confidence band, candlesticks (scenarios
  3-5) are the renderer's.

  Background:
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And the "Curves" package is installed
    And user opens curves dataset
    And the context panel is open

  Scenario: A Column level fit function reaches the other cells of its column, not another column
    When user scrolls the grid to the "styled proprtional with IC50" column
    And user clicks on the "cell 1 of styled proprtional with IC50" area of grid
    Then the current column should be "styled proprtional with IC50"
    Given "Options" pane in context panel is expanded
    Then "Level" input in context panel should have value "Column"
    And "Fit Function" input in context panel should have value "sigmoid"
    When user selects "linear" in "Fit Function" input in context panel
    And user clicks on the "cell 5 of styled proprtional with IC50" area of grid
    Then the current column should be "styled proprtional with IC50"
    And "Fit Function" input in context panel should have value "linear"
    When user scrolls the grid to the "custom fit function" column
    And user clicks on the "cell 1 of custom fit function" area of grid
    Then the current column should be "custom fit function"
    And "Fit Function" input in context panel should have value "Polynomial"
    And no errors should have been logged

  Scenario: A Cell level title stays in its cell
    When user scrolls the grid to the "styled proprtional with IC50" column
    And user clicks on the "cell 1 of styled proprtional with IC50" area of grid
    Then the current column should be "styled proprtional with IC50"
    And "Title" input in context panel should have value "Dose-Response curves"
    When user selects "Cell" in "Level" input in context panel
    And user enters "Only this one" into "Title" input in context panel
    And user clicks on the "cell 3 of styled proprtional with IC50" area of grid
    Then the current column should be "styled proprtional with IC50"
    And "Title" input in context panel should have value "Dose-Response curves"
    When user clicks on the "cell 1 of styled proprtional with IC50" area of grid
    Then "Title" input in context panel should have value "Only this one"
    And no errors should have been logged
