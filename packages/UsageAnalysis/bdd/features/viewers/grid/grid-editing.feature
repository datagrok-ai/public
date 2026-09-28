@journey @viewers @realizes:viewers.grid
Feature: Grid cell editing and the clipboard
  Editing a cell in place and moving values around: a double-click opens the editor the grid names
  `cell-editor`, Enter commits and Escape cancels, Delete clears the current cell, Allow Edit off
  refuses the edit with a read-only balloon, a typed digit opens the editor on the current cell,
  Add New Row On Last Row Edit appends a row after an edit of the last one and only while it is on,
  Editable by (Context Panel > Advanced > Permissions) lets the users it names edit the column and
  refuses the rest with "only editable by", Physical Control+Shift+C copies the current cell, Control+C
  copies the selected rows with a header line, tab-separated, in the grid's column order - a moved
  column where it was moved, a column hidden from the grid still in (the operator confirmed the
  hidden column belongs in the copy) - Control+V moves a value between cells, Control+A
  then copy then paste leaves the table alone, and Shift+Delete removes the selected rows while
  Control+Z brings them back. The permissions scenario types the login of the account the run
  signs in with, then that of the second account the library's setup keeps (`bddsecond`, or
  DATAGROK_SHARING_LOGIN), on SEVERITY where `grid-ui.md` uses Chemist of spgi-100 (operator
  decision D8). One journey on demog-1000; every scenario puts the value back. The auto-append and
  permissions scenarios run before the clipboard ones: placed after them, the grid drew no cells
  once the appended row was deleted (seen on dev in two runs, not reproduced outside the journey;
  reported to the operator).

  Not translated, and why: from `grid-ui.md` "Context Panel - Permissions", the "Pin if editable" switch is not claimed, and the refusal
  message is matched by its text, not by the login it names.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then grid should show 1000 rows
    And "Allow Edit" property of grid should be "true"
    And cell editor should be hidden

  Scenario: Enter commits the edited value and Escape cancels the next edit
    Given user listens for "d4-grid-cell-value-edited" event on grid
    When user double-clicks on the "cell 4 of AGE" area of grid
    Then cell editor should be visible
    When user presses Control+A in cell editor
    And user types "99" into cell editor
    And user presses Enter
    Then "d4-grid-cell-value-edited" event should have fired on grid
    And the value of "AGE" column in row 4 should be "99"
    And cell editor should be hidden
    When user double-clicks on the "cell 4 of AGE" area of grid
    Then cell editor should be visible
    When user presses Control+A in cell editor
    And user types "77" into cell editor
    And user presses Escape
    Then the value of "AGE" column in row 4 should be "99"
    And cell editor should be hidden
    When user double-clicks on the "cell 4 of AGE" area of grid
    And user presses Control+A in cell editor
    And user types "45" into cell editor
    And user presses Enter
    Then the value of "AGE" column in row 4 should be "45"
    And no errors should have been logged

  Scenario: Delete clears the current cell without opening the editor
    When user clicks on the "cell 5 of AGE" area of grid
    And user presses Delete
    Then "AGE" column should have missing values
    And cell editor should be hidden
    When user double-clicks on the "cell 5 of AGE" area of grid
    And user types "51" into cell editor
    And user presses Enter
    Then the value of "AGE" column in row 5 should be "51"
    And "AGE" column should have no missing values
    And no errors should have been logged

  Scenario: A read-only grid refuses the edit and warns
    When user sets "Allow Edit" property of grid to "false"
    And user double-clicks on the "cell 6 of AGE" area of grid
    Then a warning balloon containing "read-only" should have been shown
    And cell editor should be hidden
    When user sets "Allow Edit" property of grid to "true"
    And user clicks on the "cell 6 of AGE" area of grid
    And user presses 7
    Then cell editor should be visible
    When user presses Enter
    Then the value of "AGE" column in row 6 should be "7"
    When user double-clicks on the "cell 6 of AGE" area of grid
    And user presses Control+A in cell editor
    And user types "49" into cell editor
    And user presses Enter
    Then the value of "AGE" column in row 6 should be "49"
    And no errors should have been logged

  Scenario: Add New Row On Last Row Edit appends a row after an edit in the last row, and only while it is on
    When user sets "Add New Row On Last Row Edit" property of grid to "true"
    And user clicks on the "cell 3 of AGE" area of grid
    And user presses Control+End
    And user double-clicks on the "cell 1000 of AGE" area of grid
    And user presses Control+A in cell editor
    And user types "33" into cell editor
    And user presses Enter
    Then the table should have 1001 rows
    And the value of "AGE" column in row 1000 should be "33"
    When user sets "Add New Row On Last Row Edit" property of grid to "false"
    And user presses Control+End
    And user double-clicks on the "cell 1001 of AGE" area of grid
    And user types "34" into cell editor
    And user presses Enter
    Then the value of "AGE" column in row 1001 should be "34"
    And the table should have 1001 rows
    When user clicks on the "row header 1001" area of grid holding Control
    Then 1 row should be selected
    When user presses Shift+Delete
    Then the table should have 1000 rows
    When user presses Control+Home
    And user presses Control+End
    Then grid should have a "cell 1000 of AGE" area
    When user double-clicks on the "cell 1000 of AGE" area of grid
    And user presses Control+A in cell editor
    And user types "63" into cell editor
    And user presses Enter
    Then the value of "AGE" column in row 1000 should be "63"
    When user clears the row selection
    And user presses Control+Home
    Then no errors should have been logged

  Scenario: Editable by lets the users it names edit the column and refuses everyone else
    When user clicks on the "header SEVERITY" area of grid
    Given the context panel is open
    Then the context panel should show "SEVERITY"
    Given Advanced accordion header in context panel is expanded
    And Permissions accordion header in context panel is expanded
    When user enters the sharing user's login into "Editable by" input in context panel
    And user clicks on the "cell 4 of SEVERITY" area of grid
    And user presses Enter
    Then cell editor should be hidden
    When user presses Delete
    Then a warning balloon containing "only editable by" should have been shown
    And the value of "SEVERITY" column in row 4 should be "Medium"
    When user clicks on the "header SEVERITY" area of grid
    Then the context panel should show "SEVERITY"
    When user enters the current user's login into "Editable by" input in context panel
    And user clicks on the "cell 4 of SEVERITY" area of grid
    And user presses Enter
    Then cell editor should be visible
    When user presses Escape
    Then cell editor should be hidden
    When user clicks on the "header SEVERITY" area of grid
    And user clears "Editable by" input in context panel
    And user presses Enter in "Editable by" input in context panel
    And user clicks on the "cell 4 of SEVERITY" area of grid
    And user presses Enter
    Then cell editor should be visible
    When user presses Escape
    Then no errors should have been logged

  Scenario: Control+Shift+C copies the current cell
    When user clicks on the "cell 4 of AGE" area of grid
    And user presses ControlLeft+Shift+C
    Then the clipboard should have the text "45"
    And no errors should have been logged

  Scenario: Control+C copies the selected rows in the grid's column order, a hidden column included
    When user drags the "header HEIGHT" area of grid to the "header DEMOG" area
    And user picks "Hide" from the context menu of the "header RACE" area of grid
    Then the "column order" reading of grid should be "USUBJID, AGE, SEX, DIS_POP, WEIGHT, DEMOG, HEIGHT, CONTROL, STARTED, SEVERITY"
    When user clicks on the "row header 1" area of grid
    And user clicks on the "row header 5" area of grid holding Shift
    Then 5 rows should be selected
    When user presses Control+C
    Then the clipboard should hold 6 lines
    And line 1 of the clipboard should hold the tab-separated values "USUBJID, AGE, SEX, RACE, DIS_POP, WEIGHT, DEMOG, HEIGHT, CONTROL, STARTED, SEVERITY"
    And line 2 of the clipboard should hold the tab-separated values "X0273T21000300003, 26, F, Caucasian, Indigestion, 74.10, 26 C F, 174.705, false, 8/2/1990, High"
    And line 6 of the clipboard should hold the tab-separated values "X0273T21000500006, 51, F, Caucasian, Indigestion, 60.00, 51 C F, 164.980, false, 5/28/1990, None"
    When user clears the row selection
    And user drags the "header HEIGHT" area of grid to the "header DIS_POP" area
    And user picks "Order or Hide Columns..." from the context menu of the "cell 2 of AGE" area of grid
    And user clicks on the "cell 4 of x" area of Grid viewer in Order or Hide Columns dialog
    Then the "text of cell 4 of x" reading of Grid viewer in Order or Hide Columns dialog should be "true"
    When user clicks on CLOSE button in Order or Hide Columns dialog
    Then the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    And no rows should be selected
    And no errors should have been logged

  Scenario: Select all, copy and paste leaves the table as it was
    When user clicks on the "cell 1 of USUBJID" area of grid
    And user presses Control+A
    And user presses Control+C
    Then the clipboard should hold 1001 lines
    When user presses Control+V
    Then the table should have 1000 rows
    And the value of "USUBJID" column in row 1 should be "X0273T21000300003"
    And the value of "AGE" column in row 1 should be "26"
    When user clears the row selection
    Then no errors should have been logged

  Scenario: Control+V pastes a copied value into another cell
    When user clicks on the "cell 2 of AGE" area of grid
    And user presses Control+C
    And user clicks on the "cell 11 of AGE" area of grid
    And user presses Control+V
    Then the value of "AGE" column in row 11 should be "30"
    When user double-clicks on the "cell 11 of AGE" area of grid
    And user presses Control+A in cell editor
    And user types "46" into cell editor
    And user presses Enter
    Then the value of "AGE" column in row 11 should be "46"
    And no errors should have been logged

  Scenario: Shift+Delete removes the selected rows and Control+Z brings them back
    When user clicks on the "row header 4" area of grid holding Control
    And user clicks on the "row header 6" area of grid holding Control
    And user clicks on the "row header 8" area of grid holding Control
    Then 3 rows should be selected
    When user presses Shift+Delete
    Then the table should have 997 rows
    When user presses Control+Z
    Then the table should have 1000 rows
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged
