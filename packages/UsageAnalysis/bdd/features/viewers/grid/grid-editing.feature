@journey @viewers @realizes:viewers.grid
Feature: Grid cell editing and the clipboard
  Editing a cell in place and moving values around: a double-click opens the editor the grid names
  `cell-editor`, Enter commits and Escape cancels, Delete clears the current cell, Allow Edit off
  refuses the edit with a read-only balloon, a typed digit opens the editor on the current cell,
  Physical Control+Shift+C copies the current cell, Control+C / Control+V move a value between cells,
  Control+A then copy then paste leaves the table alone, and Shift+Delete removes the selected rows
  while Control+Z brings them back. One journey on demog-1000; every scenario puts the value back.

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

  Scenario: Control+Shift+C copies the current cell
    When user clicks on the "cell 4 of AGE" area of grid
    And user presses ControlLeft+Shift+C
    Then the clipboard should have the text "45"
    And no errors should have been logged

  Scenario: Control+C copies the selected rows
    When user clicks on the "row header 1" area of grid
    And user clicks on the "row header 5" area of grid holding Shift
    Then 5 rows should be selected
    When user presses Control+C
    Then the clipboard should contain the text "X0273T21000300003"
    And the clipboard should contain the text "X0273T21000500006"
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged

  Scenario: Select all, copy and paste leaves the table as it was
    When user clicks on the "cell 1 of USUBJID" area of grid
    And user presses Control+A
    And user presses Control+C
    And user presses Control+V
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
