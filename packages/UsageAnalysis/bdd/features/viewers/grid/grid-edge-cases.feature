@viewers @realizes:viewers.grid
Feature: Grid edge cases
  Guards for the grid at the edges of its settings and its data, on demog-1000 (AGE in row 4 is 45,
  SEVERITY in row 3 is None and in row 4 Medium). An HTML summary column keeps rendering after a
  column it shows is removed, and the grid scrolls to its last row with no error (TestTrack
  `Viewers/Grid/grid.md`, "Summary Column Surviving a Removed Source Column", the GROK-19942 contract).
  The rest claims what the grid should do at the edges of its settings, its pins, its sort and its
  editing; a scenario the grid does not meet yet is a known failure with its ticket: core#1137
  (Frozen Columns above the column count or below 0), core#1191 (rows pinned beyond ten), core#1193
  (Shift+Delete in a read-only grid), core#1194 (Shift+Enter), core#1195 (a blank group name),
  core#1196 (a pinned row brought back by Control+Z), core#1197 (Default HTML Form with no column)
  and core#1198 (a malformed number format). A value refused in the last row keeps the cell as it was,
  and a paste into a read-only grid changes nothing with no notice: a paste needs the cell in edit
  mode, and a double-click on a read-only cell already says it is read-only.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then grid should show 1000 rows

  @edge-cases
  Scenario: An HTML summary column keeps rendering after a column it shows is removed
    When user picks "Add > Summary Columns > Default HTML Form" from the context menu of the "cell 2 of USUBJID" area of grid
    And user clicks on All link in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then the "cell type of html" reading of grid should be "html"
    When user picks "Remove" from the context menu of the "header CONTROL" area of grid
    Then the table should not have a column "CONTROL"
    When user clicks on the "cell 2 of AGE" area of grid
    And user presses Control+End
    Then grid should have a "cell 1000 of html" area
    And the "cell type of html" reading of grid should be "html"
    And no errors should have been logged

  @edge-cases @known-failure:42 @realizes:core#1137
  Scenario: Unpin Column with Frozen Columns above the column count logs no error
    core#1137: with Frozen Columns 20 on an 11-column table every column counts as frozen, and Unpin
    Column inserts the column at position 19 of a shorter list, which throws "RangeError: Value not in
    range". Remove the tag when core#1137 is fixed and this scenario passes.
    When user sets "Frozen Columns" property of grid to "20"
    And user picks "Pin > Unpin Column" from the context menu of the "header AGE" area of grid
    Then no errors should have been logged

  @edge-cases
  Scenario Outline: Frozen Columns typed in the context panel as <typed> is not taken
    When user hovers over the "header AGE" area of grid
    And user clicks on "Edit properties (F4)" icon in grid
    Then context panel should be visible
    Given "Columns" category in context panel is expanded
    When user enters "<typed>" into "Frozen Columns" property in context panel
    Then "Frozen Columns" property of grid should be "1"
    And no errors should have been logged

    @known-failure:51 @realizes:core#1137
    Examples: The typed value is taken
      | typed |
      | 20    |
      | -1    |

  @edge-cases @known-failure:72 @realizes:core#1137
  Scenario: Frozen Columns is brought down to the column count when columns are removed
    When user sets "Frozen Columns" property of grid to "6"
    And user removes "DEMOG" column
    And user removes "CONTROL" column
    And user removes "STARTED" column
    And user removes "SEVERITY" column
    And user removes "WEIGHT" column
    And user removes "HEIGHT" column
    And user removes "DIS_POP" column
    And user removes "RACE" column
    Then the "column order" reading of grid should be "USUBJID, AGE, SEX"
    And "Frozen Columns" property of grid should be "4"
    And no errors should have been logged

  @edge-cases
  Scenario: A value refused in the last row keeps that cell as it was and keeps the added row
    When user sets "Add New Row On Last Row Edit" property of grid to "true"
    And user clicks on the "cell 3 of AGE" area of grid
    And user presses Control+End
    And user double-clicks on the "cell 1000 of AGE" area of grid
    And user presses Control+A in cell editor
    And user types "abc" into cell editor
    And user presses Enter
    Then the table should have 1001 rows
    And the value of "AGE" column in row 1000 should be "63"
    And no errors should have been logged

  @edge-cases
  Scenario: A read-only grid does not take a paste
    When user sets "Allow Edit" property of grid to "false"
    And user clicks on the "cell 2 of AGE" area of grid
    And user presses Control+C
    And user clicks on the "cell 11 of AGE" area of grid
    And user presses Control+V
    Then the value of "AGE" column in row 11 should be "46"
    And no errors should have been logged

  @edge-cases @known-failure:104 @realizes:core#1193
  Scenario: A read-only grid keeps its rows on Shift+Delete
    When user sets "Allow Edit" property of grid to "false"
    And user clicks on the "row header 5" area of grid holding Control
    Then 1 row should be selected
    When user presses Shift+Delete
    Then the table should have 1000 rows
    And no errors should have been logged

  @edge-cases
  Scenario: Empty values go last in a sort, descending and ascending
    Given user opens a table "sort-empty" with:
      | id | v |
      | a  | 2 |
      | b  |   |
      | c  | 1 |
      | d  | 3 |
    When user double-clicks on the "header v" area of grid
    Then the "sort direction" reading of grid should be "descending"
    And the "row order" reading of grid should be "4, 1, 3, 2"
    When user double-clicks on the "header v" area of grid
    Then the "sort direction" reading of grid should be "ascending"
    And the "row order" reading of grid should be "3, 1, 4, 2"
    And no errors should have been logged

  @edge-cases @known-failure:133 @realizes:core#1195
  Scenario: A column group named with spaces only is not made
    When user clicks on the "header AGE" area of grid holding Control
    And user clicks on the "header SEX" area of grid holding Control
    Given the context panel is open
    And Actions accordion header in context panel is expanded
    When user clicks on "Group columns..." text in context panel
    Then Group input in dialog should be visible
    When user types "   " into Group input in dialog
    And user clicks on OK button in dialog
    Then "AGE" column should have tag "group" equal to ""
    And "SEX" column should have tag "group" equal to ""
    And no errors should have been logged

  @edge-cases @known-failure:142 @realizes:core#1194
  Scenario: Shift+Enter selects nothing while Allow Row Selection is off
    When user sets "Allow Row Selection" property of grid to "false"
    And user clicks on the "cell 1 of SEX" area of grid
    And user presses Shift+Enter
    Then no rows should be selected
    And no errors should have been logged

  @edge-cases @known-failure:150 @realizes:core#1194
  Scenario: Shift+Enter selects only the same rows that pass the filter
    When user filters rows where "AGE" is between 20 and 30
    And user clicks on the "cell 1 of SEX" area of grid
    And user presses Shift+Enter
    Then every selected row should pass the filter
    And no errors should have been logged

  @edge-cases @known-failure:160 @realizes:core#1191
  Scenario: Twelve rows pinned show ten, and the grid says that pinned rows are hidden
    When user clicks on the "row header 1" area of grid
    And user clicks on the "row header 12" area of grid holding Shift
    Then 12 rows should be selected
    When user picks "Pin > Pin Selected Rows" from the context menu of the "cell 1 of USUBJID" area of grid
    Then the "pinned rows" reading of grid should be 10
    And a warning balloon containing "10" should have been shown

  @edge-cases @known-failure:172 @realizes:core#1191
  Scenario: Unpin All Rows unpins the pinned rows the grid does not show
    When user clicks on the "row header 1" area of grid
    And user clicks on the "row header 12" area of grid holding Shift
    And user picks "Pin > Pin Selected Rows" from the context menu of the "cell 1 of USUBJID" area of grid
    Then the "pinned rows" reading of grid should be 10
    When user picks "Pin > Unpin All Rows" from the context menu of the "cell 1 of USUBJID" area of grid
    Then the "pinned rows" reading of grid should be 0
    When user saves the layout of the current table view
    And user loads the saved layout
    Then the "pinned rows" reading of grid should be 0
    And no errors should have been logged

  @edge-cases @known-failure:184 @realizes:core#1196
  Scenario: A pinned row deleted and brought back by Control+Z is pinned again
    When user picks "Pin > Pin Row" from the context menu of the "cell 3 of USUBJID" area of grid
    Then the "pinned rows" reading of grid should be 1
    When user clicks on the "row header 3" area of grid holding Control
    And user presses Shift+Delete
    Then the table should have 999 rows
    When user presses Control+Z
    Then the table should have 1000 rows
    And the "pinned rows" reading of grid should be 1
    And no errors should have been logged

  @edge-cases @known-failure:191 @realizes:core#1197
  Scenario: Default HTML Form with no column ticked cannot be confirmed
    When user picks "Add > Summary Columns > Default HTML Form" from the context menu of the "cell 2 of USUBJID" area of grid
    Then "Select columns..." dialog should be visible
    And OK button in "Select columns..." dialog should be disabled

  @edge-cases @known-failure:198 @realizes:core#1198
  Scenario: A malformed number format marks the format field invalid
    When user picks "Format > Custom..." from the context menu of the "header AGE" area of grid
    Then Format AGE dialog should be visible
    When user enters "0.0.0" into Custom input in Format AGE dialog
    Then Custom input in Format AGE dialog should be invalid

  @edge-cases @known-failure:206 @realizes:core#1198
  Scenario: A malformed number format is not taken and the old format stays
    When user picks "Format > Custom..." from the context menu of the "header AGE" area of grid
    Then Format AGE dialog should be visible
    When user enters "0.0.0" into Custom input in Format AGE dialog
    And user clicks on OK button in Format AGE dialog
    Then "AGE" column should have tag "format" equal to ""
    And the "text of cell 1 of AGE" reading of grid should be "26"
    And no errors should have been logged

  @edge-cases
  Scenario: The sort follows a renamed column
    When user double-clicks on the "header AGE" area of grid
    Then the "sort column" reading of grid should be "AGE"
    When user renames "AGE" column to "YEARS"
    Then the "sort column" reading of grid should be "YEARS"
    And no errors should have been logged

  @edge-cases
  Scenario: A pinned row follows a renamed column through a layout
    When user picks "Pin > Pin Row" from the context menu of the "cell 3 of USUBJID" area of grid
    And user renames "USUBJID" column to "ID"
    And user saves the layout of the current table view
    And user loads the saved layout
    Then the "pinned rows" reading of grid should be 1
    And no errors should have been logged

  @edge-cases
  Scenario: A column group follows a renamed column
    When user clicks on the "header AGE" area of grid holding Control
    And user clicks on the "header SEX" area of grid holding Control
    Given the context panel is open
    And Actions accordion header in context panel is expanded
    When user clicks on "Group columns..." text in context panel
    And user types "Person" into Group input in dialog
    And user clicks on OK button in dialog
    Then grid should have a "group Person" area
    When user renames "AGE" column to "YEARS"
    Then "YEARS" column should have tag "group" equal to "Person"
    And grid should have a "group Person" area
    And no errors should have been logged

  @edge-cases
  Scenario Outline: Row Height typed in the context panel as <typed> is not taken
    When user hovers over the "header AGE" area of grid
    And user clicks on "Edit properties (F4)" icon in grid
    Then context panel should be visible
    Given "Rows" category in context panel is expanded
    When user enters "<typed>" into "Row Height" property in context panel
    Then "Row Height" property of grid should be "28"
    And no errors should have been logged

    @known-failure:249 @realizes:core#1200
    Examples: The typed value is taken although the property declares 1..500
      | typed |
      | 0     |
      | 501   |
