@journey @viewers @realizes:viewers.trellis-plot
Feature: Trellis plot legend and lifecycle
  The legend the inner viewer's colour column gives the trellis — its visibility, its four slots and
  the categories it lists — and what survives the viewer leaving and coming back: undo after the
  title-bar close, a layout round-trip through the server, and a project round-trip. One journey on
  demog-1000 with SEX by RACE and a scatter plot inside.

  A box plot inside, with Show All Categories on, keeps that setting through the legend's label and
  cross clicks (GROK-20432); the legend's dimming is read as the subset it marks current
  (`aria-selected`): nothing before the first click, the label-clicked entry alone, a second one
  added with Control, one taken away by its cross, a plain label click replacing the subset with its
  own entry, and nothing again once the last one's cross is clicked. Show All Categories is read
  from the inner look the trellis keeps and builds its cells from, where it was written — not from
  each inner box plot or from the Box plot tab of the context panel the md reads. The md's first click is a cross on an entry nothing has chosen yet; the current legend
  shows a cross only on a chosen entry (`d4.css`, `.d4-legend-item-current > .d4-legend-cross`), so
  the subset is built with the label and Control first.
  Not translated here: "Pick Up / Apply" between two trellises and "Use in Trellis" from another
  viewer (`viewer-chrome.feature` and a later round). The md's persistence tail with two trellises
  and a Selected trellis through a project is `trellis-plot-persistence.feature`.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names | SEX          |
      | Y Column Names | RACE         |
      | Viewer Type    | Scatter plot |
    Then the "cells" reading of trellis plot viewer should be 8

  Scenario: A colour column gives the trellis a legend
    Then legend of trellis plot viewer should be hidden
    When user sets "colorColumnName" inner property of trellis plot viewer to "SEX"
    And user sets "Legend Visibility" property of trellis plot viewer to "Always"
    Then legend of trellis plot viewer should be visible
    And the legend of trellis plot viewer should list 2 items
    When user sets "colorColumnName" inner property of trellis plot viewer to "RACE"
    Then the legend of trellis plot viewer should list 4 items
    And legend of trellis plot viewer should contain text "Caucasian"
    And no errors should have been logged

  Scenario: The legend takes each of the four slots
    When user sets "Legend Position" property of trellis plot viewer to "Left"
    Then the legend of trellis plot viewer should be in the "left" slot
    When user sets "Legend Position" property of trellis plot viewer to "Right"
    Then the legend of trellis plot viewer should be in the "right" slot
    When user sets "Legend Position" property of trellis plot viewer to "Top"
    Then the legend of trellis plot viewer should be in the "top" slot
    When user sets "Legend Position" property of trellis plot viewer to "Bottom"
    Then the legend of trellis plot viewer should be in the "bottom" slot
    And the legend of trellis plot viewer should list 4 items
    When user sets "Legend Visibility" property of trellis plot viewer to "Never"
    Then legend of trellis plot viewer should be hidden
    When user sets properties of trellis plot viewer:
      | Legend Visibility | Always |
      | Legend Position   | Auto   |
    Then legend of trellis plot viewer should be visible
    And no errors should have been logged

  Scenario: A box plot inside keeps Show All Categories through the legend's clicks
    When user sets "Viewer Type" property of trellis plot viewer to "Box plot"
    And user sets "showAllCategories" inner property of trellis plot viewer to "true"
    Then "showAllCategories" inner property of trellis plot viewer should be "true"
    And legend of trellis plot viewer should be visible
    And the legend of trellis plot viewer should list 6 items
    And "AS" legend item in legend of trellis plot viewer should not be selected
    And "Indigestion" legend item in legend of trellis plot viewer should not be selected
    And "PsA" legend item in legend of trellis plot viewer should not be selected
    When user clicks on "AS" item in the legend of trellis plot viewer
    Then "AS" legend item in legend of trellis plot viewer should be selected
    And "Indigestion" legend item in legend of trellis plot viewer should not be selected
    And "PsA" legend item in legend of trellis plot viewer should not be selected
    And "showAllCategories" inner property of trellis plot viewer should be "true"
    When user clicks on "Indigestion" item in the legend of trellis plot viewer holding Control
    Then "AS" legend item in legend of trellis plot viewer should be selected
    And "Indigestion" legend item in legend of trellis plot viewer should be selected
    And "PsA" legend item in legend of trellis plot viewer should not be selected
    And "showAllCategories" inner property of trellis plot viewer should be "true"
    When user clicks on the cross of "AS" item in the legend of trellis plot viewer
    Then "AS" legend item in legend of trellis plot viewer should not be selected
    And "Indigestion" legend item in legend of trellis plot viewer should be selected
    And "PsA" legend item in legend of trellis plot viewer should not be selected
    And "showAllCategories" inner property of trellis plot viewer should be "true"
    When user clicks on "PsA" item in the legend of trellis plot viewer
    Then "PsA" legend item in legend of trellis plot viewer should be selected
    And "Indigestion" legend item in legend of trellis plot viewer should not be selected
    And "AS" legend item in legend of trellis plot viewer should not be selected
    And "showAllCategories" inner property of trellis plot viewer should be "true"
    When user clicks on the cross of "PsA" item in the legend of trellis plot viewer
    Then "PsA" legend item in legend of trellis plot viewer should not be selected
    And "Indigestion" legend item in legend of trellis plot viewer should not be selected
    And "AS" legend item in legend of trellis plot viewer should not be selected
    And "PsA" legend item in legend of trellis plot viewer should not be selected
    And "UC" legend item in legend of trellis plot viewer should not be selected
    And the legend of trellis plot viewer should list 6 items
    And "showAllCategories" inner property of trellis plot viewer should be "true"
    And no errors should have been logged
    When user sets "Viewer Type" property of trellis plot viewer to "Scatter plot"
    Then the "inner viewer type" reading of trellis plot viewer should be "Scatter plot"

  Scenario: Undo brings the closed viewer back and redo closes it again
    Then the open tableview should have 1 trellis plot viewer
    When user clicks on close icon of trellis plot viewer
    Then the open tableview should have 0 trellis plot viewers
    When user presses Control+Z
    Then the open tableview should have 1 trellis plot viewer
    And the "cells" reading of trellis plot viewer should be 8
    When user presses Control+Shift+Z
    Then the open tableview should have 0 trellis plot viewers
    When user presses Control+Z
    Then the open tableview should have 1 trellis plot viewer
    And the "cells" reading of trellis plot viewer should be 8
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A layout round-trip through the server restores the split and the type
    When user sets properties of trellis plot viewer:
      | X Column Names | RACE      |
      | Y Column Names | SEX       |
      | Viewer Type    | Bar chart |
    Then the cells of trellis plot viewer should be 4 wide and 2 tall
    When user saves the layout of the current table view to the server
    And user clicks on close icon of trellis plot viewer
    Then the open tableview should have 0 trellis plot viewers
    When user loads the saved layout
    Then the open tableview should have 1 trellis plot viewer
    And properties of trellis plot viewer should be:
      | X Column Names | RACE      |
      | Y Column Names | SEX       |
      | Viewer Type    | Bar chart |
    And the cells of trellis plot viewer should be 4 wide and 2 tall
    And the "cells drawn" reading of trellis plot viewer should be 8
    And no errors should have been logged

  Scenario: A project round-trip restores the split, the type and the row source
    When user sets properties of trellis plot viewer:
      | X Column Names | SEX          |
      | Y Column Names | RACE         |
      | Viewer Type    | Scatter plot |
      | Row Source     | All          |
    Then the cells of trellis plot viewer should be 2 wide and 4 tall
    When user saves the current view as project "bdd-trellis-plot"
    And user closes all views
    And user opens the "bdd-trellis-plot" project
    Then the open tableview should have 1 trellis plot viewer
    And properties of trellis plot viewer should be:
      | X Column Names | SEX          |
      | Y Column Names | RACE         |
      | Viewer Type    | Scatter plot |
      | Row Source     | All          |
    And the cells of trellis plot viewer should be 2 wide and 4 tall
    And the "cells drawn" reading of trellis plot viewer should be 8
    And trellis plot viewer should show 1000 rows
    And no errors should have been logged

