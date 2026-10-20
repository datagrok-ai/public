@journey @viewers @realizes:viewers.line-chart
Feature: Line chart through a saved layout and a saved project
  What `line-chart.md` saves with the view and expects back, each time with the viewer closed before
  the layout is applied, so that nothing on the screen is left over from before: the configured
  chart (X = STARTED, Y = AGE and HEIGHT, split by SEX, Multi Axis, Line Width 3, spline
  interpolation), the four Selection checkboxes switched away from their defaults, and the Data
  panel's Pack Categories with Multi Axis. Then `legend-color-and-persistence.md`'s second
  scenario: a category colour picked in the legend itself comes back from a layout applied after
  the viewer was closed, and from a project. The checks read the viewer's own look and readings
  after the reload; the legend colour is taken off the column between the save and the load, so the
  layout has to bring it back. The layout and the project go through the server by the API (the md
  uses View > Layout and the ribbon's SAVE).
  One journey on demog-1000, whose SEX has F and M; the project and the layouts are deleted at the
  end.
  Not translated: the md's right-click paths to Split Columns and Multi Axis — the same properties,
  set directly; the gesture is not what these scenarios are about.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a line chart viewer with:
      | xColumnName  | AGE    |
      | yColumnNames | HEIGHT |
    Then line chart viewer should report no error

  Scenario: The configured chart comes back from a layout applied after the viewer was closed
    When user sets properties of line chart viewer:
      | xColumnName      | STARTED     |
      | yColumnNames     | AGE, HEIGHT |
      | splitColumnNames | SEX         |
      | multiAxis        | true        |
      | lineWidth        | 3           |
      | interpolation    | Spline      |
    Then the "x column" reading of line chart viewer should be "STARTED"
    And the "y columns" reading of line chart viewer should be "AGE, HEIGHT"
    And the "split columns" reading of line chart viewer should be 1
    And the "multi axis" reading of line chart viewer should be "true"
    When user saves the layout of the current table view to the server
    And user clicks on close icon of line chart viewer
    Then line chart viewer should be absent
    When user loads the saved layout
    Then line chart viewer should be visible
    And properties of line chart viewer should be:
      | xColumnName      | STARTED     |
      | yColumnNames     | AGE, HEIGHT |
      | splitColumnNames | SEX         |
      | multiAxis        | true        |
      | lineWidth        | 3           |
      | interpolation    | Spline      |
    And the "x column" reading of line chart viewer should be "STARTED"
    And the "y columns" reading of line chart viewer should be "AGE, HEIGHT"
    And the "split columns" reading of line chart viewer should be 1
    And the "categories" reading of line chart viewer should be 2
    And the "multi axis" reading of line chart viewer should be "true"
    And the "charts" reading of line chart viewer should be 1
    And line chart viewer should report no error
    And line chart viewer should be painted
    And no errors should have been logged

  Scenario: The Selection checkboxes keep their states through a layout
    When user sets properties of line chart viewer:
      | xColumnName  | AGE    |
      | yColumnNames | HEIGHT |
      | splitColumnNames |    |
      | multiAxis    | false  |
    And user drags a selection box over the "plot" area of line chart viewer
    Then some rows should be selected
    When user sets "showSelectedRows" property of line chart viewer to "false"
    Then line chart viewer should show less selection highlight than before
    When user sets "showSelectedRows" property of line chart viewer to "true"
    Then line chart viewer should show more selection highlight than before
    When user sets properties of line chart viewer:
      | showSelectedRows      | false |
      | showCurrentRowLine    | true  |
      | showMouseOverCategory | false |
      | showMouseOverRowLine  | false |
    And user saves the layout of the current table view to the server
    And user clicks on close icon of line chart viewer
    Then line chart viewer should be absent
    When user loads the saved layout
    Then line chart viewer should be visible
    And properties of line chart viewer should be:
      | showSelectedRows      | false |
      | showCurrentRowLine    | true  |
      | showMouseOverCategory | false |
      | showMouseOverRowLine  | false |
    When user sets properties of line chart viewer:
      | showSelectedRows      | true  |
      | showCurrentRowLine    | false |
      | showMouseOverCategory | true  |
      | showMouseOverRowLine  | true  |
    And user clears the row selection
    Then no errors should have been logged

  Scenario: Pack Categories off and Multi Axis on come back from a layout
    When user sets properties of line chart viewer:
      | yColumnNames   | AGE, HEIGHT |
      | packCategories | false       |
      | multiAxis      | true        |
    Then the "multi axis" reading of line chart viewer should be "true"
    When user saves the layout of the current table view to the server
    And user clicks on close icon of line chart viewer
    Then line chart viewer should be absent
    When user loads the saved layout
    Then line chart viewer should be visible
    And properties of line chart viewer should be:
      | packCategories | false |
      | multiAxis      | true  |
    And the "multi axis" reading of line chart viewer should be "true"
    And the "charts" reading of line chart viewer should be 1
    When user sets properties of line chart viewer:
      | packCategories | true   |
      | multiAxis      | false  |
      | yColumnNames   | HEIGHT |
    Then no errors should have been logged

  Scenario: A colour picked in the legend comes back from a layout and from a project
    When user sets "splitColumnNames" property of line chart viewer to "SEX"
    Then legend of line chart viewer should be visible
    When user right-clicks on "F" legend item in legend of line chart viewer
    Then "F" dialog should be visible
    When user picks the color "#9467BD" in the color picker dialog
    And user clicks on OK button in "F" dialog
    Then the categorical color of "F" in "SEX" column should be "#9467BD"
    And the "F" item in the legend of line chart viewer should be colored "#9467BD"
    When user saves the layout of the current table view to the server
    And user clicks on close icon of line chart viewer
    Then line chart viewer should be absent
    When user removes the coloring of "SEX" column
    Then "SEX" column should have no color coding
    When user loads the saved layout
    Then line chart viewer should be visible
    And the categorical color of "F" in "SEX" column should be "#9467BD"
    And the "F" item in the legend of line chart viewer should be colored "#9467BD"
    When user saves the current view as project "bdd-line-chart-layout"
    And user closes all views
    And user opens the "bdd-line-chart-layout" project
    Then line chart viewer should be visible
    And legend of line chart viewer should be visible
    And the categorical color of "F" in "SEX" column should be "#9467BD"
    And the "F" item in the legend of line chart viewer should be colored "#9467BD"
    And no errors should have been logged
