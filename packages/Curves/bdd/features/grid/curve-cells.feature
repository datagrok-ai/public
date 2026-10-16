@journey
Feature: Curve cells in the grid and their Context Panel
  Translated from TestTrack Curves/curves-grid.md. The curves demo table (System:DemoFiles/curves.csv,
  30 rows) holds five curve columns the package detects as semantic type "fit". A click on a curve
  cell shows its Options, Chart and Fit panes in the Context Panel; a double-click opens the chart in
  the Edit chart dialog. Every click is preceded by scrolling the grid to the column, and followed by
  the current column it made, so a pane is read for the cell that was meant.

  Not claimed: what a cell draws on its canvas (points, fit line, droplines, statistics, the red
  cross of a cell that cannot be drawn) — the renderer reports no readings.

  Background:
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And the "Curves" package is installed
    And user opens curves dataset
    And the context panel is open

  Scenario: The curve columns are detected as fitted curves
    Then the table should have 30 rows
    And "multiple prefit" column should have semantic type "fit"
    And "multiple styled series" column should have semantic type "fit"
    And "styled proprtional with IC50" column should have semantic type "fit"
    And "log min/max candlesticks" column should have semantic type "fit"
    And "custom fit function" column should have semantic type "fit"
    And no errors should have been logged

  Scenario: A cell with two curves shows its panes, the fit functions and the summary inputs
    When user scrolls the grid to the "multiple styled series" column
    And user clicks on the "cell 2 of multiple styled series" area of grid
    Then the current column should be "multiple styled series"
    And "Options" pane in context panel should be visible
    And "Chart" pane in context panel should be visible
    And "Fit" pane in context panel should be visible
    Given "Options" pane in context panel is expanded
    Then "Level" input in context panel should have value "Column"
    And "Level" input in context panel should offer the choice "Dataframe"
    And "Level" input in context panel should offer the choice "Cell"
    And "Fit Function" input in context panel should offer the choice "linear"
    And "Fit Function" input in context panel should offer the choice "sigmoid"
    And "Fit Function" input in context panel should offer the choice "log-linear"
    And "Fit Function" input in context panel should offer the choice "exponential"
    And "Fit Function" input in context panel should offer the choice "4pl-regression"
    And "Fit Function" input in context panel should offer the choice "4pl-dose-response"
    And "Fit Function" input in context panel should not offer the choice ""
    And "Droplines" input in context panel should contain text "IC10"
    And "Droplines" input in context panel should contain text "IC90"
    And "Statistics mode" input in context panel should be visible
    And "Aggregation" input in context panel should be visible
    And no errors should have been logged

  Scenario: A cell with one curve has nothing to summarise
    When user scrolls the grid to the "styled proprtional with IC50" column
    And user clicks on the "cell 1 of styled proprtional with IC50" area of grid
    Then the current column should be "styled proprtional with IC50"
    Given "Options" pane in context panel is expanded
    Then "Fit Function" input in context panel should be visible
    And "Statistics mode" input in context panel should not be visible
    And "Aggregation" input in context panel should not be visible
    And no errors should have been logged

  Scenario: A double-click opens the chart in the Edit chart dialog, and Close closes it
    When user scrolls the grid to the "multiple prefit" column
    And user double-clicks on the "cell 3 of multiple prefit" area of grid
    Then "Edit chart" dialog should be visible
    When user clicks on "Close" button in "Edit chart" dialog
    Then "Edit chart" dialog should be absent
    And the table should have 30 rows
    And no errors should have been logged

  Scenario: The filter panel and a pinned curve column take the curves without errors
    When user clicks on filter icon in toolbar
    Then filter panel should be visible
    And no errors should have been logged
    When user scrolls the grid to the "multiple prefit" column
    And user picks "Pin > Pin Column" from the context menu of the "header multiple prefit" area of grid
    Then the grid should pin the columns "multiple prefit"
    And no errors should have been logged
