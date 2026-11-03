@journey @viewers @realizes:viewers.tooltips
Feature: Editing the table's tooltip from a viewer
  Tooltip > Edit... on a viewer opens the Edit Tooltip dialog: a check box that ticks every column,
  a searchable list of the table's columns with a box each, the Reset group tooltip and Design
  custom tooltip... actions under it, and OK, CANCEL and the history icon in its footer. The search
  ignores case. The columns picked there and confirmed with OK become the tooltip every viewer that
  inherits the table's tooltip shows, the grid included, and stay the same in a cloned view and in
  the original (public#1208). Edit Tooltip on a viewer bound to another table lists that table's
  columns (public#373). Translated from the TestTrack case
  Tooltips/edit-tooltip, on demog-1000 rather than SPGI: the dialog and the tooltip do not depend on the
  table, and demog-1000 has no molecules to render.

  The grid is switched to "inherit from table" with Show Visible Columns In Tooltip on (the case
  says), and Show Column Names is Always on the grid and both plots, so the tooltip's
  columns are read by name; see default-tooltip-visibility for why. Between two hovers the pointer
  rests above the grid. The scatter plot adds its own axis columns to the picked ones (its Data
  Values property is Merge by default), so on demog-1000 it also lists HEIGHT, its X; the box plot and
  the grid list exactly the picked columns.

  The search is claimed with a lower-case text, a mixed-case one and a part from the middle of two
  names. Above the list the dialog also offers Use tooltip (Table here) and Show column names.

  Not translated: that the list's grid has three columns (type, name, box) — the list steps read the
  names and the boxes, not the grid's columns; and that every viewer lists the picked columns in the
  same order — the tooltip steps claim the set, not the order (MISSING.md).

  Edge cases: all on the box plot, after the scenarios above. Edit Tooltip writes the picked columns
  on OK only (`edit_tooltip.dart`), so a change made in it and cancelled leaves the tooltip as it was.
  Every column unticked and confirmed is a tooltip set to no column (`tooltip.dart`, "Columns
  explicitly set to none"): the dialog opened again has every column unticked, and one ticked column
  brings the tooltip back. A tooltip column removed from the table (a calculated copy of SEX, removed
  from the grid header's menu) is skipped and the others still show, without an error. That a
  renamed one should show under its new name is a known failure (core#1215), claimed in
  tooltips-known-failures.feature while the ticket is open.
  The check box above the list reads unchecked while only some columns are ticked, so unticking it
  then changes nothing; a scenario that starts from a partial tooltip ticks it first, then unticks it.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user sets properties of grid:
      | Show Tooltip                    | inherit from table |
      | Show Column Names               | Always             |
      | Show Visible Columns In Tooltip | true               |
    And user adds a scatter plot viewer
    And user adds a box plot viewer
    And user sets properties of scatter plot viewer:
      | showLabels | Always |
    And user sets properties of box plot viewer:
      | showLabels | Always |

  Scenario: Tooltip > Edit... opens the Edit Tooltip dialog with every column ticked
    When user picks "Tooltip > Edit..." from the context menu of scatter plot viewer
    Then "Edit Tooltip" dialog should be visible
    And "Use tooltip" input in "Edit Tooltip" dialog should have value "Table"
    And "Show column names" input in "Edit Tooltip" dialog should be visible
    And checkbox in "Edit Tooltip" dialog should be checked
    And the column list of "Edit Tooltip" dialog should be exactly "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    And the "AGE" column should be checked in the column list of "Edit Tooltip" dialog
    And the "SEVERITY" column should be checked in the column list of "Edit Tooltip" dialog
    And "Reset group tooltip" label in "Edit Tooltip" dialog should be visible
    And "Design custom tooltip..." label in "Edit Tooltip" dialog should be visible
    And OK button in "Edit Tooltip" dialog should be visible
    And CANCEL button in "Edit Tooltip" dialog should be visible
    And "History" icon in "Edit Tooltip" dialog should be visible
    And no errors should have been logged

  Scenario: The search narrows the list whatever the case of the text
    When user types "race" into "Search" input in "Edit Tooltip" dialog
    Then the column list of "Edit Tooltip" dialog should be exactly "RACE"
    And the column list of "Edit Tooltip" dialog should start with "RACE"
    When user clears "Search" input in "Edit Tooltip" dialog
    And user types "SEVer" into "Search" input in "Edit Tooltip" dialog
    Then the column list of "Edit Tooltip" dialog should be exactly "SEVERITY"
    When user clears "Search" input in "Edit Tooltip" dialog
    And user types "ht" into "Search" input in "Edit Tooltip" dialog
    Then the column list of "Edit Tooltip" dialog should be exactly "HEIGHT, WEIGHT"
    When user clears "Search" input in "Edit Tooltip" dialog
    Then the column list of "Edit Tooltip" dialog should be exactly "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    And no errors should have been logged

  @realizes:public#1208
  Scenario: The columns picked with OK become the tooltip of the plots and the grid, in a cloned view too (public#1208)
    When user unchecks checkbox in "Edit Tooltip" dialog
    Then the "AGE" column should not be checked in the column list of "Edit Tooltip" dialog
    And the "RACE" column should not be checked in the column list of "Edit Tooltip" dialog
    When user toggles the "AGE" column in the column list of "Edit Tooltip" dialog
    And user toggles the "SEX" column in the column list of "Edit Tooltip" dialog
    And user toggles the "WEIGHT" column in the column list of "Edit Tooltip" dialog
    And user clicks on OK button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    And properties of scatter plot viewer should be:
      | X | HEIGHT |
      | Y | WEIGHT |
    When user hovers over the "marker of row 11" area of scatter plot viewer
    Then the tooltip should show columns in the order "HEIGHT, WEIGHT, AGE, SEX"
    When user moves the pointer away from grid
    Then tooltip should be hidden
    When user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns in the order "AGE, SEX, WEIGHT"
    When user moves the pointer away from grid
    Then tooltip should be hidden
    When user hovers over the "cell 11 of AGE" area of grid
    Then the tooltip should show columns in the order "AGE, SEX, WEIGHT"
    When user moves the pointer away from grid
    And user picks "View > Layout > Clone View" from the top menu
    Then the "demog-1000 copy" view should be current
    When user hovers over the "marker of row 11" area of scatter plot viewer
    Then the tooltip should show columns in the order "HEIGHT, WEIGHT, AGE, SEX"
    When user moves the pointer away from grid
    And user closes the current view
    And user switches to the "demog-1000" view
    And user hovers over the "marker of row 11" area of scatter plot viewer
    Then the tooltip should show columns in the order "HEIGHT, WEIGHT, AGE, SEX"
    When user moves the pointer away from grid
    Then no errors should have been logged

  @realizes:public#373
  Scenario: Edit Tooltip on a viewer bound to another table lists that table's columns (public#373)
    Given user opens spgi dataset
    And user switches to the "demog-1000" table view
    When user sets "table" property of scatter plot viewer to "spgi-100"
    Then scatter plot viewer should be bound to table "spgi-100"
    When user picks "Tooltip > Edit..." from the context menu of scatter plot viewer
    Then "Edit Tooltip" dialog should be visible
    And the column list of "Edit Tooltip" dialog should start with "Id"
    And the column list of "Edit Tooltip" dialog should not list "USUBJID"
    And the column list of "Edit Tooltip" dialog should not list "AGE"
    When user clicks on CANCEL button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    And no errors should have been logged

  @realizes:core#474
  Scenario: View > Edit Tooltip keeps setting the table's tooltip for the plots (core#474)
    When user picks "View > Edit Tooltip..." from the top menu
    Then "Edit Tooltip" dialog should be visible
    And "Use tooltip" input in "Edit Tooltip" dialog should have value "Table"
    When user toggles the "WEIGHT" column in the column list of "Edit Tooltip" dialog
    And user clicks on OK button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    When user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns in the order "AGE, SEX"
    When user moves the pointer away from grid
    Then no errors should have been logged

  @edge-cases
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
    Then the tooltip should show columns in the order "AGE, SEX"
    And the tooltip should not show columns "RACE"
    When user moves the pointer away from grid
    Then no errors should have been logged

  @edge-cases
  Scenario: Every column unticked and confirmed is kept, and one ticked column brings the tooltip back
    When user picks "Tooltip > Edit..." from the context menu of box plot viewer
    Then "Edit Tooltip" dialog should be visible
    When user checks checkbox in "Edit Tooltip" dialog
    Then the "RACE" column should be checked in the column list of "Edit Tooltip" dialog
    When user unchecks checkbox in "Edit Tooltip" dialog
    Then the "AGE" column should not be checked in the column list of "Edit Tooltip" dialog
    And the "RACE" column should not be checked in the column list of "Edit Tooltip" dialog
    When user clicks on OK button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    When user picks "Tooltip > Edit..." from the context menu of box plot viewer
    Then "Edit Tooltip" dialog should be visible
    And the "AGE" column should not be checked in the column list of "Edit Tooltip" dialog
    And the "SEX" column should not be checked in the column list of "Edit Tooltip" dialog
    When user toggles the "AGE" column in the column list of "Edit Tooltip" dialog
    And user clicks on OK button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    When user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns in the order "AGE"
    When user moves the pointer away from grid
    And user picks "Tooltip > Edit..." from the context menu of box plot viewer
    Then "Edit Tooltip" dialog should be visible
    When user toggles the "SEX" column in the column list of "Edit Tooltip" dialog
    And user clicks on OK button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    When user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns in the order "AGE, SEX"
    When user moves the pointer away from grid
    Then no errors should have been logged

  @edge-cases
  Scenario: A tooltip column removed from the table is skipped and the others still show
    Given user adds a calculated column "SEX_copy" with formula "${SEX}"
    When user picks "Tooltip > Edit..." from the context menu of box plot viewer
    Then "Edit Tooltip" dialog should be visible
    When user checks checkbox in "Edit Tooltip" dialog
    And user unchecks checkbox in "Edit Tooltip" dialog
    Then the "WEIGHT" column should not be checked in the column list of "Edit Tooltip" dialog
    When user toggles the "AGE" column in the column list of "Edit Tooltip" dialog
    And user toggles the "SEX_copy" column in the column list of "Edit Tooltip" dialog
    And user toggles the "WEIGHT" column in the column list of "Edit Tooltip" dialog
    And user clicks on OK button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    When user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns "AGE, SEX_copy, WEIGHT"
    When user moves the pointer away from grid
    And user scrolls the grid to the "SEX_copy" column
    And user picks "Remove" from the context menu of the "header SEX_copy" area of grid
    And user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns in the order "AGE, WEIGHT"
    When user moves the pointer away from grid
    Then no errors should have been logged
