@journey
Feature: Curve options and the level they are set at
  Translated from TestTrack Curves/curves-options.md, scenarios 1 and 6. An option set at Column level
  in the Options pane of one curve cell is what the other cells of that column show, and it is kept
  by a layout applied to the table opened again, and by a project.

  The option is read where a person reads it, in the Options pane of another cell. Not translated:
  that another column keeps its own fit function and the Cell level of scenario 2 — after a click on
  a curve cell of another column, or after Level is set to Cell, the Context Panel often holds the
  plain table cell instead of the curve panes, and the panes cannot be read there. Not claimed: what
  an option changes on the canvas — connect dots, markers, the confidence band, candlesticks
  (scenarios 3-5 of the md) are pictures. The layout and the project are saved and opened through the
  library's steps, the way the rest of the suite sets a scene.

  Background:
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And the "Curves" package is installed
    And no project named "CurvesOptions-{run}" is on the server
    And the layouts named "curves-options-{run}" are deleted when the feature ends
    And user opens curves dataset
    And the context panel is open

  Scenario: A Column level fit function reaches the other cells of its column
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
    And no errors should have been logged

  Scenario: A Column level fit function is kept by a layout applied to the table opened again
    When user saves the layout of the current table view as "curves-options-{run}"
    And user closes all views
    Given user opens curves dataset
    And the context panel is open
    When user applies the layout "curves-options-{run}" to the current table view
    And user scrolls the grid to the "styled proprtional with IC50" column
    And user clicks on the "cell 4 of styled proprtional with IC50" area of grid
    Then the current column should be "styled proprtional with IC50"
    Given "Options" pane in context panel is expanded
    Then "Fit Function" input in context panel should have value "linear"
    And no errors should have been logged

  Scenario: A Column level fit function is kept by a project
    When user saves the current view as project "CurvesOptions-{run}"
    And user closes all views
    And user opens the "CurvesOptions-{run}" project and waits for its table
    And user scrolls the grid to the "styled proprtional with IC50" column
    And user clicks on the "cell 4 of styled proprtional with IC50" area of grid
    Then the current column should be "styled proprtional with IC50"
    Given "Options" pane in context panel is expanded
    Then "Fit Function" input in context panel should have value "linear"
    And no errors should have been logged
