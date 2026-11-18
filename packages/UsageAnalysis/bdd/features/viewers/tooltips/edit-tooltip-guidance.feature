@journey @viewers @realizes:viewers.tooltips
Feature: The Edit Tooltip dialog says what it edits and what the tooltip will show
  On a scatter plot the dialog's Use tooltip offers two tooltips: "Table (shared)", the one column list
  of the table that every viewer set to inherit from table shows (up to 7 columns), and "Scatter plot
  only", the plot's own list (up to 20). Hovering the captions of Use tooltip and Show column names
  explains them. Under the column list a status line counts the checked columns and warns, before OK,
  about what the tooltip will not show: the columns past its 7 or 20, or nothing at all. Switching Use
  tooltip back and forth keeps what was checked in each. "Scatter plot only" with every column checked
  shows all 11 columns of demog-1000, where the table's tooltip shows 7. Design custom tooltip...
  designs the tooltip Use tooltip names when it is clicked, and the dialog reopened over that design
  says OK replaces it. Remove Group Tooltip is offered disabled while no viewer is the group tooltip;
  once one is, the dialog offers Reset group tooltip, which removes it on OK.

  On demog-1000 (11 columns). The scatter plot adds its X and Y (HEIGHT and WEIGHT) to any tooltip, and
  its Show Column Names is Always, so the tooltip's columns are read by name. What a layout keeps of the
  table's tooltip (its columns, Show column names and Hide) is the xamgle client test "A layout keeps the
  table's Show column names and Hide"; the 7 and 20 columns are "The table's tooltip lists 7 of its
  columns, a viewer's own 20".

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer
    And user sets properties of scatter plot viewer:
      | showLabels | Always |

  Scenario: Use tooltip and Show column names explain themselves on hover
    When user picks "Tooltip > Edit..." from the context menu of scatter plot viewer
    Then "Edit Tooltip" dialog should be visible
    And "Use tooltip" input in "Edit Tooltip" dialog should have value "Table (shared)"
    When user hovers over "Use tooltip" label in "Edit Tooltip" dialog
    Then tooltip should contain text "every viewer whose Show Tooltip is"
    And tooltip should contain text "Scatter plot only"
    When user hovers over "Show column names" label in "Edit Tooltip" dialog
    Then tooltip should contain text "only where values could be mixed up"
    When user clicks on CANCEL button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    And no errors should have been logged

  Scenario: The status line counts the checked columns, warns past the table's 7 and on none, and each tooltip keeps its checks
    When user picks "Tooltip > Edit..." from the context menu of scatter plot viewer
    Then "Edit Tooltip" dialog should contain text "Shows the first 7 of 11 checked"
    When user unchecks checkbox in "Edit Tooltip" dialog
    Then "Edit Tooltip" dialog should contain text "Nothing checked: no row tooltip"
    When user toggles the "AGE" column in the column list of "Edit Tooltip" dialog
    Then "Edit Tooltip" dialog should contain text "1 of 11 checked"
    When user selects "Scatter plot only" in "Use tooltip" input in "Edit Tooltip" dialog
    Then "Edit Tooltip" dialog should contain text "Nothing checked: only X and Y show"
    When user checks checkbox in "Edit Tooltip" dialog
    Then "Edit Tooltip" dialog should contain text "11 of 11 checked"
    When user selects "Table (shared)" in "Use tooltip" input in "Edit Tooltip" dialog
    Then "Edit Tooltip" dialog should contain text "1 of 11 checked"
    And the "AGE" column should be checked in the column list of "Edit Tooltip" dialog
    And the "SEX" column should not be checked in the column list of "Edit Tooltip" dialog
    When user clicks on CANCEL button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    And no errors should have been logged

  Scenario: "Scatter plot only" with every column checked shows all 11, where the table's tooltip shows 7
    When user picks "Tooltip > Edit..." from the context menu of scatter plot viewer
    And user selects "Scatter plot only" in "Use tooltip" input in "Edit Tooltip" dialog
    And user checks checkbox in "Edit Tooltip" dialog
    And user clicks on OK button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    And "Show Tooltip" property of scatter plot viewer should be "show custom tooltip"
    When user hovers over the "marker of row 11" area of scatter plot viewer
    Then the tooltip should show columns "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    When user moves the pointer away from grid
    And user picks "Tooltip > Edit..." from the context menu of scatter plot viewer
    And user selects "Table (shared)" in "Use tooltip" input in "Edit Tooltip" dialog
    Then "Edit Tooltip" dialog should contain text "Shows the first 7 of 11 checked"
    When user clicks on OK button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    And "Show Tooltip" property of scatter plot viewer should be "inherit from table"
    When user hovers over the "marker of row 11" area of scatter plot viewer
    Then the tooltip should show columns "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT"
    When user moves the pointer away from grid
    Then no errors should have been logged

  Scenario: Design custom tooltip... designs the scatter plot's own tooltip once Use tooltip says so, and the dialog then says OK replaces it
    When user picks "Tooltip > Edit..." from the context menu of scatter plot viewer
    And user selects "Scatter plot only" in "Use tooltip" input in "Edit Tooltip" dialog
    And user clicks on "Design custom tooltip..." label in "Edit Tooltip" dialog
    And user clicks on "CLOSE AND APPLY" button
    Then the "demog-1000" view should be current
    And "Show Tooltip" property of scatter plot viewer should be "show custom tooltip"
    When user picks "Tooltip > Edit..." from the context menu of scatter plot viewer
    Then "Use tooltip" input in "Edit Tooltip" dialog should have value "Scatter plot only"
    And "Edit Tooltip" dialog should contain text "A designed tooltip is in use: OK replaces it"
    When user clicks on CANCEL button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    And no errors should have been logged

  Scenario: Remove Group Tooltip is disabled until a viewer is the group tooltip, which Reset group tooltip removes on OK
    When user opens the context menu of scatter plot viewer
    Then the open menu should list "Tooltip > Remove Group Tooltip" as disabled
    When user closes the context menu
    And user picks "Tooltip > Use as Group Tooltip" from the context menu of scatter plot viewer
    And user opens the context menu of scatter plot viewer
    Then the open menu should list "Tooltip > Remove Group Tooltip" as enabled
    When user closes the context menu
    And user picks "Tooltip > Edit..." from the context menu of scatter plot viewer
    And user clicks on "Reset group tooltip" label in "Edit Tooltip" dialog
    Then "Edit Tooltip" dialog should contain text "The group tooltip is removed on OK"
    When user clicks on OK button in "Edit Tooltip" dialog
    Then the "Edit Tooltip" dialog should close
    When user opens the context menu of scatter plot viewer
    Then the open menu should list "Tooltip > Remove Group Tooltip" as disabled
    When user closes the context menu
    Then no errors should have been logged
