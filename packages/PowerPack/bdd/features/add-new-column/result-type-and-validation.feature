@journey @realizes:powerpack.dialogs.add-new-column
Feature: Add New Column: the result type, a formula the dialog refuses, and a name already taken
  The Add New Column dialog on a three-row table written in the feature. The Type choice shows the
  type the formula computes while it is on auto, and a type picked by hand is the type of the added
  column; a picked type the formula cannot give is named under the editor and OK is disabled. A
  formula naming a column the table has not got is refused the same way. A vector function inside
  arithmetic (CumSum(...) * 2, public#4253) is a valid formula and adds its column. A name the table
  already has is kept apart: the dialog says which name the column will get, and OK adds it under
  that name. Plain text makes the formula a text template ("${s} + ${s}" gives "a + a", where auto
  would compute the expression). Over a formula that computes a number (public#4295), string adds
  the numbers as text and plain text the formula's text; the type is picked before the formula, so
  OK is pressed once the preview computed under it is shown. With the formula left empty, the
  preview shows the empty column and OK adds it with the picked type (that it carries no formula is
  the package test's). Translated from help/transform/add-new-column.md (Validation, the type
  picker) and the QA round of 13 Oct 2026 on calculated columns.

  Not translated: what each formula computes on every row, and the cases of the validator, are the
  PowerPack package tests "Add new column: validation" and ApiTests "DataFrame: Calculated columns";
  this feature claims what the dialog shows and what its OK adds, once.

  Background:
    Given user is logged in
    And user opens a table "calc" with:
      | x | s |
      | 1 | a |
      | 2 | b |
      | 4 | c |
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible

  Scenario: On auto, the Type choice names the type the formula computes
    When user types "${x} / 2" into formula editor
    Then the preview grid should show "${x} / 2" computed as numbers
    And column type input should have value "auto (double)"
    When user clears formula editor
    And user types "${x} > 1" into formula editor
    Then column type input should have value "auto (bool)"
    And no errors should have been logged

  Scenario: A type picked by hand is the type of the added column
    When user clears formula editor
    And user types "dbl" into column name input
    And user types "${x} * 2" into formula editor
    And user selects "double" in column type input
    And user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And "dbl" column should have type "double"
    And the value of "dbl" column in row 3 should be "8"
    And no errors should have been logged

  Scenario: A picked type the formula cannot give is named, and OK is disabled
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    When user types "${s}" into formula editor
    And user selects "int" in column type input
    Then formula error should contain text "doesn't match with current column type (int)"
    And OK button in "Add New Column" dialog should be disabled
    When user selects "string" in column type input
    Then OK button in "Add New Column" dialog should be enabled
    When user clicks on CANCEL button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And no errors should have been logged

  Scenario: A formula naming a column the table has not got is refused
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    When user types "${nope} + 1" into formula editor
    Then formula error should contain text "nope"
    And OK button in "Add New Column" dialog should be disabled
    When user clicks on CANCEL button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And no errors should have been logged

  Scenario: A vector function inside arithmetic is a valid formula and adds its column
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    When user types "running" into column name input
    And user types "CumSum(${x}) * 2" into formula editor
    Then the preview grid should show "running" computed as numbers
    When user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And the value of "running" column in row 1 should be "2"
    And the value of "running" column in row 3 should be "14"
    And no errors should have been logged

  Scenario: Plain text makes the formula a text template filled in for every row
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    When user types "template" into column name input
    And user types "${s} + ${s}" into formula editor
    And user selects "plain text" in column type input
    And user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And "template" column should have type "string"
    And the value of "template" column in row 1 should be "a + a"
    And the value of "template" column in row 3 should be "c + c"
    And no errors should have been logged

  Scenario: A name the table already has is kept apart, and OK adds it under the new name
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    When user types "${x} + 10" into formula editor
    Then the preview grid should show "${x} + 10" computed as numbers
    When user types "x" into column name input
    Then formula error should contain text "Column named x already exists. Column will be named x (2)"
    When user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And the table should have a column "x (2)"
    And the value of "x (2)" column in row 3 should be "14"
    And the value of "x" column in row 3 should be "4"
    And no errors should have been logged

  Scenario: A number the formula computes is added as text with string, and as the formula's text with plain text
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    When user types "as text" into column name input
    And user selects "string" in column type input
    And user types "${x} + 1" into formula editor
    Then the preview grid should show "as text" computed as numbers
    When user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And "as text" column should have type "string"
    And the value of "as text" column in row 3 should be "5"
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    When user types "as template" into column name input
    And user selects "plain text" in column type input
    And user types "${x} + 1" into formula editor
    Then the "text of cell 3 of as template" reading of preview grid viewer should be "4 + 1"
    When user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And "as template" column should have type "string"
    And the value of "as template" column in row 3 should be "4 + 1"
    And no errors should have been logged

  Scenario: With the formula left empty, OK adds an empty column of the picked type
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    When user selects "int" in column type input
    And user types "notes" into column name input
    Then the "text of cell 1 of notes" reading of preview grid viewer should be ""
    When user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And "notes" column should have type "int"
    And the value of "notes" column in row 1 should be ""
    And the value of "notes" column in row 3 should be ""
    And no errors should have been logged
