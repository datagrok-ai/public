@viewers @realizes:viewers.tooltips
Feature: Tooltips — the known failures
  Temporary: the known failures of the tooltips while their tickets are open. Each scenario is
  expected to fail on its pinned step; when its ticket is fixed, remove the tag and move the scenario
  back into edit-tooltip.feature as an ordinary @edge-cases scenario.

  A tooltip column renamed should show under its new name (core#1215). Today the tooltip keeps the
  old name in the table tag (`tooltip.dart`), nothing renames it, and the column drops out of the
  tooltip. The column is renamed in its Column Properties. On demog-1000 with the grid on "inherit
  from table" and Show Column Names Always on the grid and the box plot, as in edit-tooltip.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user sets properties of grid:
      | Show Tooltip                    | inherit from table |
      | Show Column Names               | Always             |
      | Show Visible Columns In Tooltip | true               |
    And user adds a box plot viewer
    And user sets properties of box plot viewer:
      | showLabels | Always |

  @edge-cases @known-failure:42 @realizes:core#1215
  Scenario: A tooltip column renamed is shown under its new name
    When user picks "Tooltip > Edit..." from the context menu of box plot viewer
    Then "Edit Tooltip" dialog should be visible
    When user checks checkbox in "Edit Tooltip" dialog
    And user unchecks checkbox in "Edit Tooltip" dialog
    Then the "WEIGHT" column should not be checked in the column list of "Edit Tooltip" dialog
    When user toggles the "AGE" column in the column list of "Edit Tooltip" dialog
    And user toggles the "WEIGHT" column in the column list of "Edit Tooltip" dialog
    And user clicks on OK button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    When user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns in the order "AGE, WEIGHT"
    When user moves the pointer away from grid
    And user scrolls the grid to the "WEIGHT" column
    And user picks "Column Properties..." from the context menu of the "header WEIGHT" area of grid
    And user enters "MASS" into "New name:" input in "WEIGHT" dialog
    And user clicks on OK button in "WEIGHT" dialog
    And user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns in the order "AGE, MASS"
