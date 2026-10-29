@viewers @edge-cases @realizes:viewers.trellis-plot
Feature: Trellis plot edge cases: Tiles Per Row 0, keys after a table swap, a label click under a filter
  Guards for the property bounds and the changes under an open trellis plot (no rows, a renamed
  split column and Tiles Per Row above its maximum: `viewer-rules-edge-cases.feature`).
  Closing the trellis gives back the rows its clicked cell filtered out (a filter that outlives its
  viewer is a defect, core#233). Three defects of core#1057 are pinned as known failures. Tiles Per
  Row declares 1..10 (`trellis_plot_look.dart`), yet 0 typed in the context panel is taken, and on a
  one-column trellis every render then throws a truncating division by zero. After the trellis is
  set to another table and back, the arrow keys no longer move the current cell. A click on a
  category label under a filter selects rows the filter removed, while a cell click selects only the
  rows that pass it. Each tag goes when core#1057 is fixed. Data: demog-1000, SEX by RACE, scatter
  plots inside.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names | SEX          |
      | Y Column Names | RACE         |
      | Viewer Type    | Scatter plot |
    Then the cells of trellis plot viewer should be 2 wide and 4 tall
    And the "rows shown" reading of trellis plot viewer should be 1000

  @known-failure:35 @realizes:core#1057
  Scenario: Tiles Per Row 0 typed in the context panel on a one-column trellis logs no error
    core#1057: Tiles Per Row 0 on a one-column trellis throws a truncating division by zero on every
    render. Remove the tag when core#1057 is fixed.
    When user sets properties of trellis plot viewer:
      | Y Column Names |  |
    Then the "y categories" reading of trellis plot viewer should be 1
    When user clicks on settings icon of trellis plot viewer
    Then context panel should be visible
    Given "Misc" category in context panel is expanded
    When user enters "0" into "Tiles Per Row" property in context panel
    Then no errors should have been logged
    And trellis plot viewer should report no error

  @known-failure:56 @realizes:core#1057
  Scenario: The arrow keys still move the current cell after the trellis is set to another table and back
    core#1057: the trellis' key subscription is cancelled on the first detach and never added again,
    so after a table swap the arrows no longer move the current cell. Remove the tag when core#1057
    is fixed.
    Given user opens cars dataset
    And user switches to the "demog-1000" table view
    When user sets "table" property of trellis plot viewer to "cars"
    Then trellis plot viewer should be bound to table "cars"
    When user sets "table" property of trellis plot viewer to "demog-1000"
    And user sets properties of trellis plot viewer:
      | X Column Names | SEX          |
      | Y Column Names | RACE         |
      | Viewer Type    | Scatter plot |
    Then trellis plot viewer should be bound to table "demog-1000"
    When user clicks on the "cell F | Caucasian" area of trellis plot viewer
    Then the "current cell" reading of trellis plot viewer should be "F | Caucasian"
    When user presses ArrowRight in trellis plot viewer
    Then the "current cell" reading of trellis plot viewer should be "M | Caucasian"

  @known-failure:66 @realizes:core#1057
  Scenario: A category label click under a filter selects only the rows that pass the filter
    core#1057: a category label click selects the rows of the category over every row, the
    filtered-out ones included. Remove the tag when core#1057 is fixed.
    Given user clears the row selection
    When user filters rows where "AGE" is between 20 and 40
    And user clicks on the "x label F" area of trellis plot viewer
    Then some rows should be selected
    And every selected row should pass the filter

  Scenario: Closing the trellis gives back the rows its clicked cell filtered out
    When user sets "On Click" property of trellis plot viewer to "Filter"
    And user clicks on the "cell F | Caucasian" area of trellis plot viewer
    Then 480 rows should pass the filter
    When user clicks on close icon of trellis plot viewer
    Then the open tableview should have 0 trellis plot viewers
    And all rows should pass the filter
    And no errors should have been logged
