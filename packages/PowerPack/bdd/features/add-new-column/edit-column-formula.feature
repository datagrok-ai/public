@journey @realizes:powerpack.dialogs.add-new-column
Feature: Edit Column Formula: a calculated column renamed and given a new formula in the dialog, and recalculated by a grid edit
  A three-row table written in the feature, with dbl = ${x} * 2, a column over it, next = ${dbl} + 1,
  and a running total, run = CumSum(${x}), all made through the platform's API (the dialog that adds
  them is add-new-column.feature's subject). The Formula pane of a column's context panel opens the
  Edit Column Formula dialog with the column's name and formula. A new name and a new formula there
  rename and recalculate the column, and the column over it follows both. A value edited in the grid
  recalculates the running total (GROK-21113). A scatter plot on the column redraws with the values a
  new formula gives, and plots them again when Y is picked anew (public#434). Translated from the QA round of 13 Oct 2026 on
  calculated columns; the Formula pane's own Apply is formula-refreshing.feature.

  The first scenario picks dbl by its header; the later ones make the column current through the
  platform's API, since the context panel drops a click made while it is frozen after a dialog or an
  edit (GROK-21024, kept by design; there is no signal for its end) and the pick is not their
  subject. The dialogs that only look are cancelled, and the one that changes the column comes last
  of the dialogs on dbl. The Edit in dialog button is in the pane before the pane's editor is, and a
  click on it then opens nothing, so the pane's formula is claimed first. The last scenario makes a
  column in the Add New Column dialog and claims that the Edit Column Formula preview, once computed
  (it drops the column the formula does not use), shows it under its own name, not "<name> (2)"
  (public#4293). It comes last because Escape does not close the dialog while the formula editor
  has the focus, so a failure there would leave it open over the next scenario.

  Background:
    Given user is logged in
    And user opens a table "calc" with:
      | x | s |
      | 1 | a |
      | 2 | b |
      | 4 | c |
    And user adds a calculated column "dbl" with formula "${x} * 2"
    And user adds a calculated column "next" with formula "${dbl} + 1"
    And user adds a calculated column "run" with formula "CumSum(${x})"
    And the context panel is open

  Scenario: The Formula pane opens the dialog with the column's name and formula
    When user clicks on the "header dbl" area of grid
    Then the context panel should show "dbl"
    Given Formula pane in context panel is expanded
    Then formula pane editor in Formula pane in context panel should hold the code "${x} * 2"
    When user clicks on "Edit in dialog" button in context panel
    Then "Edit Column Formula" dialog should be visible
    And column name input should have value "dbl"
    And formula editor should hold the code "${x} * 2"
    When user clicks on CANCEL button in "Edit Column Formula" dialog
    Then "Edit Column Formula" dialog should be hidden
    And no errors should have been logged

  Scenario: The dialog of the column over it starts from that column's formula, and CANCEL drops an edit
    Given the "next" column is the current object
    Then the context panel should show "next"
    Given Formula pane in context panel is expanded
    Then formula pane editor in Formula pane in context panel should hold the code "${dbl} + 1"
    When user clicks on "Edit in dialog" button in context panel
    Then "Edit Column Formula" dialog should be visible
    And column name input should have value "next"
    And formula editor should hold the code "${dbl} + 1"
    When user replaces the code of formula editor with "${x}"
    And user clicks on CANCEL button in "Edit Column Formula" dialog
    Then "Edit Column Formula" dialog should be hidden
    And "next" column should have tag "formula" equal to "${dbl} + 1"
    And no errors should have been logged

  Scenario: A new name and a new formula rename and recalculate the column, and the column over it follows
    Given the "dbl" column is the current object
    Then the context panel should show "dbl"
    Given Formula pane in context panel is expanded
    Then formula pane editor in Formula pane in context panel should hold the code "${x} * 2"
    When user clicks on "Edit in dialog" button in context panel
    Then "Edit Column Formula" dialog should be visible
    When user types "triple" into column name input
    And user replaces the code of formula editor with "${x} * 3"
    And user clicks on OK button in "Edit Column Formula" dialog
    Then "Edit Column Formula" dialog should be hidden
    And the table should have a column "triple"
    And the table should not have a column "dbl"
    And "triple" column should have tag "formula" equal to "${x} * 3"
    And the value of "triple" column in row 3 should be "12"
    And "next" column should have tag "formula" equal to "${triple} + 1"
    And the value of "next" column in row 3 should be "13"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A value edited in the grid recalculates the running total
    Then the value of "run" column in row 3 should be "7"
    When user double-clicks on the "cell 1 of x" area of grid
    And user presses Control+A in cell editor
    And user types "10" into cell editor
    And user presses Enter
    Then the value of "x" column in row 1 should be "10"
    And the value of "run" column in row 1 should be "10"
    And the value of "run" column in row 3 should be "16"
    And the value of "triple" column in row 1 should be "30"
    And the value of "next" column in row 1 should be "31"
    And no errors should have been logged

  @realizes:public#434
  Scenario: A scatter plot on a calculated column redraws after its formula is changed in the dialog, and keeps the new values (public#434)
    Given user adds a scatter plot viewer with:
      | X | x      |
      | Y | triple |
    When user remembers the "y axis max" reading of scatter plot viewer
    Given the "triple" column is the current object
    Then the context panel should show "triple"
    Given Formula pane in context panel is expanded
    Then formula pane editor in Formula pane in context panel should hold the code "${x} * 3"
    When user clicks on "Edit in dialog" button in context panel
    Then "Edit Column Formula" dialog should be visible
    When user replaces the code of formula editor with "${x} * 30"
    And user clicks on OK button in "Edit Column Formula" dialog
    Then "Edit Column Formula" dialog should be hidden
    And the "y axis max" reading of scatter plot viewer should not be as remembered
    When user remembers the "y axis max" reading of scatter plot viewer
    And user sets "Y" property of scatter plot viewer to "x"
    Then the "y axis max" reading of scatter plot viewer should not be as remembered
    When user sets "Y" property of scatter plot viewer to "triple"
    Then the "y axis max" reading of scatter plot viewer should be as remembered
    And no errors should have been logged

  Scenario: The preview of a column made in the dialog shows it under its own name
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    When user types "half" into column name input
    And user types "${x} / 2" into formula editor
    And user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And the table should have a column "half"
    Given the "half" column is the current object
    Then the context panel should show "half"
    Given Formula pane in context panel is expanded
    Then formula pane editor in Formula pane in context panel should hold the code "${x} / 2"
    When user clicks on "Edit in dialog" button in context panel
    Then "Edit Column Formula" dialog should be visible
    And the preview grid should not show "s"
    And the preview grid should not show "half (2)"
    And the preview grid should show "half" computed as numbers
    When user clicks on CANCEL button in "Edit Column Formula" dialog
    Then "Edit Column Formula" dialog should be hidden
    And no errors should have been logged
