@journey @viewers @realizes:viewers.trellis-plot
Feature: Trellis plot viewer filter, menus and a second undo cycle
  The viewer's own Filter formula narrows the rows its cells draw without touching the table's
  filter, and clearing it gives them all back. A right-click on a cell opens the trellis's menu with
  a group named after the inner viewer, holding that viewer's own items; To Script > To JavaScript
  prints the call that rebuilds the trellis. Closing the viewer, undoing and redoing twice in a row
  leaves no error. Translated from the "Viewer filter formula", "Context menu", "To Script" and
  "Undo/redo" sections of TestTrack Viewers/TrellisPlot/trellis-plot.md. One journey on demog-1000
  with SEX by RACE and a scatter plot inside; every scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names | SEX          |
      | Y Column Names | RACE         |
      | Viewer Type    | Scatter plot |
    Then the "cells" reading of trellis plot viewer should be 8
    And the "rows shown" reading of trellis plot viewer should be 1000

  Scenario: The viewer's Filter formula narrows its rows and leaves the table's filter alone
    When user sets "Filter" property of trellis plot viewer to "${AGE} > 40"
    Then the "rows shown" reading of trellis plot viewer should be 635
    And all rows should pass the filter
    When user sets "Filter" property of trellis plot viewer to ""
    Then the "rows shown" reading of trellis plot viewer should be 1000
    And no errors should have been logged

  Scenario: A cell's context menu holds the inner viewer's group
    When user right-clicks on the "cell body F | Caucasian" area of trellis plot viewer
    Then the open menu should list "Scatter plot > Lasso Tool"
    And the open menu should list "Scatter plot > Markers"
    And the open menu should list "Scatter plot > Selection"
    And the open menu should list "General > Clone"
    And the open menu should list "Properties..."
    When user closes the context menu
    Then no errors should have been logged

  Scenario: To Script > To JavaScript prints the call that rebuilds the trellis
    Given the package autostarts have completed
    When user picks "To Script > To JavaScript" from the context menu of trellis plot viewer
    Then balloon should contain text "addViewer"
    And balloon should contain text "Trellis"
    And no errors should have been logged

  Scenario: The inner viewer's tab of the context panel changes every cell
    When user clicks on settings icon of trellis plot viewer
    Then context panel should be visible
    When user clicks on "Scatter plot" tab in context panel
    Given "X Axis" category in context panel is expanded
    When user remembers the "cell signature F | Caucasian" reading of trellis plot viewer
    And user remembers the "cell signature M | Asian" reading of trellis plot viewer
    And user selects "WEIGHT" in "X" property in context panel
    Then "xColumnName" inner property of trellis plot viewer should be "WEIGHT"
    And the "cell signature F | Caucasian" reading of trellis plot viewer should not be as remembered
    And the "cell signature M | Asian" reading of trellis plot viewer should not be as remembered
    And no errors should have been logged

  Scenario: Two cycles of undo and redo after the title-bar close leave no error
    When user clicks on close icon of trellis plot viewer
    Then the open tableview should have 0 trellis plot viewers
    When user presses Control+Z
    Then the open tableview should have 1 trellis plot viewer
    And the "cells drawn" reading of trellis plot viewer should be 8
    When user presses Control+Shift+Z
    Then the open tableview should have 0 trellis plot viewers
    When user presses Control+Z
    Then the open tableview should have 1 trellis plot viewer
    And the "cells drawn" reading of trellis plot viewer should be 8
    When user presses Control+Shift+Z
    Then the open tableview should have 0 trellis plot viewers
    And no error or warning balloon should have been shown
    And no errors should have been logged
