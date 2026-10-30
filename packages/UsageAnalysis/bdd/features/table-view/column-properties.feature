@journey
Feature: The Column Properties dialog edits a column's tags, and the grid shows what they mean
  F2 on a column (shortcuts.md, visualize/viewers/grid.md "Column properties") or Column Properties... in the
  column header's menu (datagrok/concepts/table.md "Metadata") opens a dialog titled with the column's name:
  New name, Tags (a key-value row per tag the column has, the platform's own included, with + inserting an empty row
  above and - removing the row), the Edit category colors... and Change type... links.
  The standard tags it edits are in govern/catalog/tags.md: format (a number format such as "0.000") changes
  how the grid shows the values, description and units show in the column header's tooltip. OK applies, and
  Ctrl+Z puts the column back as it was (help/datagrok/navigation/undo.md: "editing its properties").

  Left out on purpose: the rename itself (PowerPack add-new-column/persistence-northwind and -sources,
  viewers/form, viewers/nx) and its undo (undo-redo.feature); a malformed format (core#1198, Olena's
  edge cases); Change type... (change-type.feature); the format matrix (ddt's format tests). Nothing is put
  on the server.

  Background:
    Given user is logged in
    And user opens a table "props" with:
      | name | city | weight |
      | Ann  | Oslo | 1.5    |
      | Bob  | Rome | 2.25   |
      | Cy   | Oslo | 3      |

  Scenario: F2 opens the dialog of the current column with its name and the format the platform gave it
    When user clicks on the "cell 1 of weight" area of grid
    And user presses F2
    Then "weight" dialog should be visible
    And "New name:" input in "weight" dialog should have value "weight"
    And "Key 1" text input in "weight" dialog should have value "format"
    And "Value 1" text input in "weight" dialog should have value "0.00"
    When user presses Escape
    Then "weight" dialog should be hidden
    And no errors should have been logged

  Scenario: A format tag changes how the grid shows the values, and Ctrl+Z and redo follow it
    Then the "weight" cell of row 1 should be displayed as "1.50"
    When user picks "Column Properties..." from the context menu of the "header weight" area of grid
    Then "Key 1" text input in "weight" dialog should have value "format"
    When user enters "0.000" into "Value 1" text input in "weight" dialog
    And user clicks on OK button in "weight" dialog
    Then the "weight" dialog should close
    And the "weight" cell of row 1 should be displayed as "1.500"
    And the "weight" cell of row 2 should be displayed as "2.250"
    And the "name" cell of row 1 should be displayed as "Ann"
    When user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then the "weight" cell of row 1 should be displayed as "1.50"
    When user presses Control+Shift+Z
    Then the "weight" cell of row 1 should be displayed as "1.500"
    When user presses Control+Z
    Then the "weight" cell of row 1 should be displayed as "1.50"
    And no errors should have been logged

  Scenario: Description and units added as two tag rows show in the header tooltip
    When user picks "Column Properties..." from the context menu of the "header weight" area of grid
    And user clicks on first "Add row" button in "weight" dialog
    And user enters "description" into "Key 1" text input in "weight" dialog
    And user enters "Body weight" into "Value 1" text input in "weight" dialog
    And user clicks on first "Add row" button in "weight" dialog
    Then "Key 2" text input in "weight" dialog should have value "description"
    When user enters "units" into "Key 1" text input in "weight" dialog
    And user enters "kg" into "Value 1" text input in "weight" dialog
    And user clicks on OK button in "weight" dialog
    Then the "weight" dialog should close
    And "weight" column should have units "kg"
    When user hovers over the "header weight" area of grid
    Then tooltip should contain text "Body weight"
    And tooltip should contain text "units: kg"
    When user moves the pointer away from grid
    And user hovers over the "header city" area of grid
    Then tooltip should contain text "city"
    And tooltip should not contain text "Body weight"
    When user moves the pointer away from grid
    And user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then "weight" column should have units ""
    And no errors should have been logged

  Scenario: The dialog shows a tag set earlier, and removing its row removes the tag
    When user picks "Column Properties..." from the context menu of the "header city" area of grid
    And user clicks on first "Add row" button in "city" dialog
    And user enters "units" into "Key 1" text input in "city" dialog
    And user enters "town" into "Value 1" text input in "city" dialog
    And user clicks on OK button in "city" dialog
    Then "city" column should have units "town"
    When user picks "Column Properties..." from the context menu of the "header city" area of grid
    Then "Key 1" text input in "city" dialog should have value "units"
    And "Value 1" text input in "city" dialog should have value "town"
    When user clicks on first "Remove row" button in "city" dialog
    And user clicks on OK button in "city" dialog
    Then the "city" dialog should close
    And "city" column should have units ""
    And no errors should have been logged

  Scenario: A name another column already has is refused
    When user picks "Column Properties..." from the context menu of the "header weight" area of grid
    And user enters "City" into "New name:" input in "weight" dialog
    Then "New name:" input in "weight" dialog should be invalid
    When user hovers over "New name:" input in "weight" dialog
    Then tooltip should contain text "Column 'City' already exists."
    When user enters "mass" into "New name:" input in "weight" dialog
    Then "New name:" input in "weight" dialog should be valid
    When user presses Escape
    Then "weight" dialog should be hidden
    And the table should have the columns "name, city, weight"
    And no errors should have been logged

  Scenario: Edit category colors... is offered for a text column and not for a number column
    When user picks "Column Properties..." from the context menu of the "header weight" area of grid
    Then "Edit category colors..." label in "weight" dialog should be disabled
    When user presses Escape
    And user picks "Column Properties..." from the context menu of the "header city" area of grid
    Then "Edit category colors..." label in "city" dialog should be enabled
    When user presses Escape
    Then "city" dialog should be hidden
    And no errors should have been logged
