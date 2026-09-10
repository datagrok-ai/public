@journey @viewers @realizes:viewers.heat-map
Feature: Heat map navigation and grid mode
  Which row a click lands on, how an Alt-drag and a double-click move the two sliders, and what
  the viewer becomes when Is Heatmap is turned off.
  The spec this replaces read the sliders' geometry off the `cx`/`cy` attributes of the SVG
  circles that draw their handles; the grid now reports `x scroll span` and `y scroll span` — the
  share of the track the window covers, 1 when everything fits — so a zoom and a reset are two
  numbers. A click used to assert only that `currentRowIdx` had changed; `current column` says
  which column it landed in, which is the half that could have been wrong. And "Is Heatmap off
  redraws as a plain grid" used to be a repaint of at least 1000 pixels; it is `row height` going
  from under a pixel to 28 and the per-cell areas replacing the column bands.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a heat map viewer
    Then the "is heatmap" reading of heat map viewer should be "true"
    And the "row height" reading of heat map viewer should be between 0 and 8
    And heat map viewer should be painted

  Scenario: A click on a column band makes a row of that column current
    Then the "current row" reading of heat map viewer should be 1
    And the "current column" reading of heat map viewer should be "USUBJID"
    When user clicks on the "column AGE" area of heat map viewer
    Then the "current column" reading of heat map viewer should be "AGE"
    And the "current row" reading of heat map viewer should differ from before
    And the table should have a current row
    When user makes row 1 current
    Then the "current row" reading of heat map viewer should be 1
    And no errors should have been logged

  Scenario: An Alt-drag zooms both sliders in and a double-click on each resets its own axis
    Then the "x scroll span" reading of heat map viewer should be 1
    And the "y scroll span" reading of heat map viewer should be 1
    When user drags a zoom box over the "column AGE" area of heat map viewer
    Then the "x scroll span" reading of heat map viewer should be lower than before
    And the "y scroll span" reading of heat map viewer should be lower than before
    And the "row height" reading of heat map viewer should be higher than before
    And heat map viewer should have repainted
    When user double-clicks on the "y scroll slider" area of heat map viewer
    Then the "y scroll span" reading of heat map viewer should be 1
    And the "x scroll span" reading of heat map viewer should be between 0 and 0.99
    When user double-clicks on the "x scroll slider" area of heat map viewer
    Then the "x scroll span" reading of heat map viewer should be 1
    And the "y scroll span" reading of heat map viewer should be 1
    And no errors should have been logged

  Scenario: Is Heatmap off draws the same table as a spreadsheet
    Then heat map viewer should have a "column AGE" area
    And heat map viewer should not have a "cell 1 of AGE" area
    When user sets "isHeatmap" property of heat map viewer to "false"
    Then the "is heatmap" reading of heat map viewer should be "false"
    And the "row height" reading of heat map viewer should be between 20 and 40
    And the "y scroll span" reading of heat map viewer should be lower than before
    And heat map viewer should have a "cell 1 of AGE" area
    And the "text of cell 1 of AGE" reading of heat map viewer should be "26"
    And heat map viewer should not have a "column AGE" area
    And heat map viewer should have repainted
    When user sets "isHeatmap" property of heat map viewer to "true"
    Then the "is heatmap" reading of heat map viewer should be "true"
    And no errors should have been logged

  @known-failure
  Scenario: Is Heatmap on again brings the whole table back on screen (grid_look.dart:435)
    # `isHeatmap = true` flips the mode flag and `refreshGrid()` rebuilds the columns, but it
    # passes `updateVertScroll: false, keepVisualRange: false`, so the vertical scroll window the
    # grid mode left behind (about 34 rows of 1000) is kept. `_rowHeight` in heat map mode is
    # `contentBox.height / (maxRow - minRow + 1)` (grid_core.dart:1686), so it stays at the grid's
    # 28 px, the per-cell areas stay instead of the column bands, and `y scroll span` stays at
    # 0.034 — the heat map does not come back. The spec this replaces asserted a repaint of at
    # least 1000 pixels in each direction, which the chrome change alone produced.
    # Left last: a known failure aborts before its restore step.
    Given user sets "isHeatmap" property of heat map viewer to "false"
    Then the "row height" reading of heat map viewer should be between 20 and 40
    When user sets "isHeatmap" property of heat map viewer to "true"
    Then the "is heatmap" reading of heat map viewer should be "true"
    And the "row height" reading of heat map viewer should be between 0 and 8
    And the "y scroll span" reading of heat map viewer should be 1
    And heat map viewer should have a "column AGE" area
