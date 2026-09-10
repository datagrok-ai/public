@journey @viewers @realizes:viewers.correlation-plot
Feature: Correlation plot — a cell as a click, hover and menu target
  What a cell answers to: a click that records the pair and the value, the context panel that names
  it, a double-click that opens it as a real scatter plot, the tooltip with the pair's own scatter
  inside, the pinned name column that hovers the column's statistics instead, the cell menu, and
  "Open as table", which hands the whole matrix over as a four-row table.
  Every gesture here lands on a region the plot reports — `cell HEIGHT x AGE`, `row header HEIGHT`,
  `pinned band` — where the spec this replaces ran a six-attempt calibration loop that moved its own
  `pinnedW` by a cell width and its `headerH` by four pixels until a probe click produced the pair it
  wanted, and then reused those constants in a second file. The click's event arguments are not
  something the runtime can read, so the pair and the value are read back as `last clicked cell` and
  `last clicked value` instead of counting an event: that value is the cell's own float32
  (-0.23482920229434967) while `correlation of` recomputes in double (-0.2348292062338333), which is
  why both are bounded rather than compared for equality.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a correlation plot viewer
    Then the "cells" reading of correlation plot viewer should be 16
    And the "show pearson r" reading of correlation plot viewer should be "true"
    And the "last clicked cell" reading of correlation plot viewer should be ""
    And correlation plot viewer should be painted

  Scenario: A click records which cell it hit and the value that cell holds
    When user clicks on the "cell HEIGHT x AGE" area of correlation plot viewer
    Then the "last clicked cell" reading of correlation plot viewer should be "HEIGHT x AGE"
    And the "last clicked value" reading of correlation plot viewer should be between -0.2349 and -0.2348
    And the "correlation of HEIGHT and AGE" reading of correlation plot viewer should be between -0.2349 and -0.2348
    And the "current column" reading of correlation plot viewer should be "HEIGHT"
    When user clicks on the "cell WEIGHT x AGE" area of correlation plot viewer
    Then the "last clicked cell" reading of correlation plot viewer should be "WEIGHT x AGE"
    And the "last clicked value" reading of correlation plot viewer should be between 0.0647 and 0.0649
    And no errors should have been logged

  Scenario: The click puts that pair in the context panel with a scatter plot of its own
    When user clicks on the "cell HEIGHT x AGE" area of correlation plot viewer
    Then context panel should contain the text "HEIGHT vs AGE"
    And "Scatter plot" section in context panel should be present
    And no errors should have been logged

  Scenario: A double-click opens the pair as a real scatter plot, and closing it puts the view back
    Then the open tableview should have 0 scatter plot viewers
    When user double-clicks on the "cell WEIGHT x AGE" area of correlation plot viewer
    Then the open tableview should have 1 scatter plot viewer
    And "xColumnName" property of scatter plot viewer should be "WEIGHT"
    And "yColumnName" property of scatter plot viewer should be "AGE"
    When user clicks on close icon of scatter plot viewer
    Then the open tableview should have 0 scatter plot viewers
    And no errors should have been logged

  Scenario: Ignore Double Click suppresses the same gesture on the same cell
    When user sets "ignoreDoubleClick" property of correlation plot viewer to "true"
    And user double-clicks on the "cell WEIGHT x AGE" area of correlation plot viewer
    Then the open tableview should have 0 scatter plot viewers
    When user sets "ignoreDoubleClick" property of correlation plot viewer to "false"
    And user double-clicks on the "cell WEIGHT x AGE" area of correlation plot viewer
    Then the open tableview should have 1 scatter plot viewer
    When user clicks on close icon of scatter plot viewer
    Then the open tableview should have 0 scatter plot viewers
    And no errors should have been logged

  Scenario: The cell menu is the grid's, with the plot's own items on it
    # The whole top level, in order: the plot adds "Open as table", "Show Pearson R" and "Columns"
    # to what the inner grid puts on a cell. Correlation Type, Show Tooltip and Ignore Double Click
    # are NOT here — they live only in the mirror under "Properties...", where every property of a
    # Dart viewer is repeated in a zero-size submenu; a group label that occurs twice cannot be
    # opened by name (the library's opener waits on the first item carrying it, the hidden one), so
    # those three are claimed through their properties instead, in correlation-plot.feature.
    When user right-clicks on the "cell HEIGHT x AGE" area of correlation plot viewer
    Then the open menu should list "Show Pearson R"
    And the open menu should list "Open as table"
    And the open menu should list "Columns"
    And the open menu should list "Tooltip"
    And the open menu should list "Grid"
    And the open menu should list "'AGE' column"
    When user closes the context menu
    Then no errors should have been logged

  Scenario: Show Pearson R from the cell menu does what the property does
    When user picks "Show Pearson R" from the context menu of the "cell HEIGHT x AGE" area of correlation plot viewer
    Then the "show pearson r" reading of correlation plot viewer should be "false"
    And the "text of cell HEIGHT x AGE" reading of correlation plot viewer should be ""
    And the "cell width" reading of correlation plot viewer should be 20
    When user picks "Show Pearson R" from the context menu of the "cell HEIGHT x AGE" area of correlation plot viewer
    Then the "show pearson r" reading of correlation plot viewer should be "true"
    And the "text of cell HEIGHT x AGE" reading of correlation plot viewer should be "-0.23"
    And no errors should have been logged

  Scenario: Open as table hands over the matrix, coefficients and all (GROK-19053)
    When user picks "Open as table" from the context menu of the "cell HEIGHT x AGE" area of correlation plot viewer
    Then table "corr" should be open
    And table "corr" should have 4 rows
    When user switches to the "corr" table view
    Then the table should have a column "__name"
    And the table should have a column "AGE"
    And the value of "__name" column in row 2 should be "HEIGHT"
    And every value of "AGE" column should lie between -1 and 1
    And "AGE" column should have its maximum in row 3
    When user closes the current view
    And user switches to the "demog-1000" table view
    Then the "cells" reading of correlation plot viewer should be 16
    And no errors should have been logged

  Scenario: Narrowed to a fraction of its width, the matrix scrolls and keeps the name column
    Then correlation plot viewer should not have a "horz scroll" area
    And the "columns shown" reading of correlation plot viewer should be 6
    When user resizes correlation plot viewer to 190 by 400
    Then correlation plot viewer should have a "horz scroll" area
    And the "columns shown" reading of correlation plot viewer should be lower than before
    And correlation plot viewer should have a "row header HEIGHT" area
    And correlation plot viewer should have a "pinned band" area
    When user restores the size of correlation plot viewer
    Then correlation plot viewer should not have a "horz scroll" area
    And the "columns shown" reading of correlation plot viewer should be 6
    And no errors should have been logged

  Scenario: The cell tooltip names the pair, its coefficient and the plot drawn inside it (GROK-20125)
    When user hovers over the "cell HEIGHT x AGE" area of correlation plot viewer
    Then exactly one tooltip should be shown
    And tooltip should contain the text "Pearson R: -0.235"
    And the "tooltip plot x column" reading of correlation plot viewer should be "HEIGHT"
    And the "tooltip plot y column" reading of correlation plot viewer should be "AGE"
    When user hovers over the "cell WEIGHT x HEIGHT" area of correlation plot viewer
    Then exactly one tooltip should be shown
    And tooltip should contain the text "Pearson R: 0.412"
    And the "tooltip plot x column" reading of correlation plot viewer should be "WEIGHT"
    And the "tooltip plot y column" reading of correlation plot viewer should be "HEIGHT"
    When user moves the pointer away from correlation plot viewer
    Then no errors should have been logged

  Scenario: The pinned name column hovers the column's statistics, not a coefficient
    When user hovers over the "row header HEIGHT" area of correlation plot viewer
    Then exactly one tooltip should be shown
    And tooltip should contain the text "min:"
    And tooltip should contain the text "nulls: 128"
    And tooltip should not contain the text "Pearson R"
    When user moves the pointer away from correlation plot viewer
    Then no errors should have been logged

  @known-failure
  Scenario: Show Tooltip off still shows the cell tooltip (correlation_plot_core.dart:57)
    # `showTooltip` is declared on the look (correlation_plot_look.dart:30, "Shows the tooltip with
    # the corresponding scatter plot inside") and the context menu's Tooltip > Visible toggles it
    # (correlation_plot_core.dart:129), but the `onCellTooltip` handler at :57 never reads it: it
    # builds the tooltip and calls `tooltip.showElement` unconditionally. Turning the setting off
    # changes nothing a user can see. Left last in the journey because a known failure aborts before
    # its own restore step.
    When user sets "showTooltip" property of correlation plot viewer to "false"
    And user moves the pointer away from correlation plot viewer
    And user hovers over the "cell HEIGHT x AGE" area of correlation plot viewer
    Then tooltip should not contain the text "Pearson R"
    When user sets "showTooltip" property of correlation plot viewer to "true"
    Then no errors should have been logged
