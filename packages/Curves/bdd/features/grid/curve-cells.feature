@journey
Feature: Curve cells in the grid and their Context Panel
  Translated from TestTrack Curves/curves-grid.md. The curves demo table (System:DemoFiles/curves.csv,
  30 rows) holds five curve columns. A click on a curve cell shows its Options, Chart and Fit panes in
  the Context Panel, its Options set at the Column level, which offers Dataframe and Cell as well
  (scenario 1, step 4); a double-click opens the chart in the Edit chart dialog. Every click is preceded by
  scrolling the grid to the column, and followed by the current column it made, so a pane is read for
  the cell that was meant.

  In the package's tests instead: the five columns detected as curves (detector-tests), the fit
  functions and droplines the Options pane offers ("a custom fit function is offered once…", "…droplines
  only for a fit that has asymptotes") and the summary inputs a one-curve cell does without ("a
  single-series cell has nothing to aggregate", panel-renderer-tests). Not claimed: what a cell draws on
  its canvas (points, fit line, droplines, statistics, the red cross of a cell that cannot be drawn) -
  the renderer's, tested by the package.

  Background:
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And the "Curves" package is installed
    And user opens curves dataset
    And the context panel is open

  Scenario: A cell with two curves shows its panes and the level its options are set at
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
    And "Statistics mode" input in context panel should be visible
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
