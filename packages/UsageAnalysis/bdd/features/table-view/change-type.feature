@journey
Feature: Change Type converts a column from its header menu or from Column Properties, and Ctrl+Z puts it back
  Follows help/transform/change-column-type.md: right-click a column header and pick Change Type... to open
  the Change column type dialog (Columns, the current Type read-only, New type, Format, Not convertible, and
  the Convertible / Not convertible / Empty counts with icons to filter or select those rows); OK converts
  the column, a value that cannot be converted takes the Not convertible value, and Ctrl+Z undoes the change.
  The same dialog opens from the Change type... link of the Column Properties dialog. Each scenario puts the
  column back with Ctrl+Z (and checks it), so the next one starts from the table as written here.

  Left out on purpose: a matrix of source and target types and of parse formats is the ddt tests of
  ChangeColumnsType (core/shared/ddt); the list New type offers for several columns at once is that
  matrix too. The counts beside Convertible / Not convertible / Empty have no name to read them by: the
  filter icon beside them is claimed instead, and the counts are in MISSING.md. The select icon is not
  claimed: it renders with no glyph (MISSING.md, Findings 30 Oct 2026). Nothing is put on the server.

  Background:
    Given user is logged in
    And user opens a table "ct" with:
      | name | code |
      | a    | 10   |
      | b    | 20   |
      | c    | n/a  |
      | d    |      |

  Scenario: The header menu opens the dialog with the column's current type, and OK converts it to the new type
    Then "code" column should have type "string"
    When user picks "Change Type..." from the context menu of the "header code" area of grid
    Then "Change column type" dialog should be visible
    And "Type" input in "Change column type" dialog should have value "string"
    And "Type" input in "Change column type" dialog should be disabled
    When user selects "int" in "New type" input in "Change column type" dialog
    And user clicks on OK button in "Change column type" dialog
    Then the "Change column type" dialog should close
    And "code" column should have type "int"
    And the value of "code" column in row 2 should be "20"
    And the value of "name" column in row 2 should be "b"
    When user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then "code" column should have type "string"
    And the value of "code" column in row 3 should be "n/a"
    And no errors should have been logged

  Scenario: A value that cannot be converted takes the Not convertible value
    When user picks "Change Type..." from the context menu of the "header code" area of grid
    And user selects "int" in "New type" input in "Change column type" dialog
    And user enters "-1" into "Not convertible" input in "Change column type" dialog
    And user clicks on OK button in "Change column type" dialog
    Then the "Change column type" dialog should close
    And "code" column should have type "int"
    And the value of "code" column in row 3 should be "-1"
    And the value of "code" column in row 1 should be "10"
    When user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then "code" column should have type "string"
    And the value of "code" column in row 3 should be "n/a"
    And no errors should have been logged

  Scenario: The filter icon beside Not convertible keeps only the rows that cannot be converted
    When user picks "Change Type..." from the context menu of the "header code" area of grid
    And user selects "int" in "New type" input in "Change column type" dialog
    And user clicks on second "Filter" icon in "Change column type" dialog
    Then 1 row should pass the filter
    And all rows where "code" is "n/a" should pass the filter
    When user clicks on CANCEL button in "Change column type" dialog
    Then the "Change column type" dialog should close
    And "code" column should have type "string"
    When user resets the filter
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: Change type... in Column Properties opens the same dialog, and redo converts the column again
    When user picks "Column Properties..." from the context menu of the "header code" area of grid
    Then "code" dialog should be visible
    When user clicks on "Change type..." label in "code" dialog
    Then "Change column type" dialog should be visible
    When user selects "double" in "New type" input in "Change column type" dialog
    And user clicks on OK button in "Change column type" dialog
    Then the "Change column type" dialog should close
    And "code" column should have type "double"
    When user presses Escape
    Then "code" dialog should be hidden
    When user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then "code" column should have type "string"
    When user presses Control+Shift+Z
    Then "code" column should have type "double"
    And the value of "code" column in row 1 should be "10"
    When user presses Control+Z
    Then "code" column should have type "string"
    And no errors should have been logged
