@viewers @edge-cases @realizes:viewers.tooltips
Feature: Tooltip edge cases: Edit Tooltip cancelled, every column unticked, a tooltip column removed
  Guards for the table's tooltip as the Edit Tooltip dialog sets it. The dialog writes the picked
  columns on OK only (`edit_tooltip.dart`), so a change made in it and then cancelled leaves the
  tooltip as it was. Every column unticked and confirmed is a tooltip explicitly set to no column
  (`tooltip.dart`, "Columns explicitly set to none"), which is not the same as no tooltip set: the
  box plot then shows no tooltip on its marker, where one ticked column brings the tooltip back. A
  column of the tooltip removed from the table is skipped and the remaining ones still show, without
  an error; a renamed one should show under its new name, a known failure with no ticket yet. Data: demog-1000 with a box plot that inherits the table's tooltip (it lists exactly the
  table's tooltip columns), Show Column Names Always, so the tooltip's columns are read by name; the
  tooltip is first set to AGE, SEX and WEIGHT in the dialog.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a box plot viewer
    And user sets properties of box plot viewer:
      | showLabels | Always |
    When user picks "Tooltip > Edit..." from the context menu of box plot viewer
    Then "Edit Tooltip" dialog should be visible
    When user unchecks checkbox in "Edit Tooltip" dialog
    And user toggles the "AGE" column in the column list of "Edit Tooltip" dialog
    And user toggles the "SEX" column in the column list of "Edit Tooltip" dialog
    And user toggles the "WEIGHT" column in the column list of "Edit Tooltip" dialog
    And user clicks on OK button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    When user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns in the order "AGE, SEX, WEIGHT"
    When user moves the pointer away from grid

  Scenario: A change made in Edit Tooltip and cancelled leaves the tooltip as it was
    When user picks "Tooltip > Edit..." from the context menu of box plot viewer
    Then "Edit Tooltip" dialog should be visible
    When user toggles the "AGE" column in the column list of "Edit Tooltip" dialog
    And user toggles the "RACE" column in the column list of "Edit Tooltip" dialog
    Then the "AGE" column should not be checked in the column list of "Edit Tooltip" dialog
    And the "RACE" column should be checked in the column list of "Edit Tooltip" dialog
    When user clicks on CANCEL button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    When user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns in the order "AGE, SEX, WEIGHT"
    And the tooltip should not show columns "RACE"
    And no errors should have been logged

  Scenario: Every column unticked and confirmed leaves no column in the tooltip
    When user picks "Tooltip > Edit..." from the context menu of box plot viewer
    Then "Edit Tooltip" dialog should be visible
    When user checks checkbox in "Edit Tooltip" dialog
    Then the "RACE" column should be checked in the column list of "Edit Tooltip" dialog
    When user unchecks checkbox in "Edit Tooltip" dialog
    Then the "AGE" column should not be checked in the column list of "Edit Tooltip" dialog
    And the "RACE" column should not be checked in the column list of "Edit Tooltip" dialog
    When user clicks on OK button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    When user hovers over the "marker" area of box plot viewer
    Then tooltip should be hidden
    When user moves the pointer away from grid
    And user picks "Tooltip > Edit..." from the context menu of box plot viewer
    Then "Edit Tooltip" dialog should be visible
    When user toggles the "AGE" column in the column list of "Edit Tooltip" dialog
    And user clicks on OK button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    When user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns in the order "AGE"
    And no errors should have been logged

  Scenario: A tooltip column removed from the table is skipped and the others still show
    When user removes "SEX" column
    And user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns in the order "AGE, WEIGHT"
    And no errors should have been logged

  @known-failure:79 @realizes:core#1215
  Scenario: A tooltip column renamed is shown under its new name
    A renamed column is followed (`viewer-rules-edge-cases.feature`); the tooltip keeps the old name in the
    table tag (`tooltip.dart`), nothing renames it, and the column drops out of the tooltip.
    When user renames "WEIGHT" column to "MASS"
    And user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns in the order "AGE, SEX, MASS"
    And no errors should have been logged

  Scenario: An aggregation row added in Edit Aggregated Tooltip takes a column by default
    When user adds a line chart viewer with:
      | xColumnName  | AGE    |
      | yColumnNames | WEIGHT |
    Then the "aggregated" reading of line chart viewer should be "true"
    And "aggTooltipColumns" property of line chart viewer should be ""
    When user picks "Tooltip > Edit..." from the context menu of line chart viewer
    Then "Edit Aggregated Tooltip" dialog should be visible
    When user clicks on "Add aggregation" button in "Edit Aggregated Tooltip" dialog
    And user clicks on OK button in "Edit Aggregated Tooltip" dialog
    Then the "Edit Aggregated Tooltip" dialog should close
    And "aggTooltipColumns" property of line chart viewer should be "avg(AGE)"
    And no errors should have been logged
