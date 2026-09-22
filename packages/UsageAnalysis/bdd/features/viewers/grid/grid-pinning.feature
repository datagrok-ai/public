@viewers @realizes:viewers.grid
Feature: Grid pinned rows and pinned columns
  Everything the Pin menu does, from both of its places: on a data cell it pins the row by that
  cell's value (Pin Row, Unpin Row, Unpin All Rows, and Pin Selected Rows when the row is part of
  the selection); on a header it pins the column (Pin Column, Pin <n> Columns for a column
  selection, Unpin Column, Unpin All Columns). The grid reports the pinned rows as `pinned rows` and
  lists the cells it drew, pinned ones included, as `cell <r> of <col>` / `row header <r>` in
  table rows, so a claim about which rows are pinned is a claim about which rows the grid still
  draws after it has scrolled to the end; the frozen columns are the `Frozen Columns` property
  (the row header counts as one) and the current cell is `current row` / `current column`. The
  pinned-rows selection case is `grid-rows-select-filter-navigate.md` scenario 6 (the rows are
  pinned from a data cell's menu, which pins by the cell's value, where the spec right-clicks the
  row number); the rest is what
  the product offers beyond the TestTrack specs. Not a journey: every scenario starts on a fresh
  demog-1000, because Unpin All Columns only moves the freeze line back and leaves a pinned column
  where pinning put it, so one scenario's pinning would decide the next one's column order. USUBJID
  is unique, SEX is not; rows 990, 995 and 998 are X0273T51070200002, X0273T51080100021 and
  X0273T51080200011.

  Not claimed: which row is current after Control+End or ArrowDown once rows are pinned. With N
  rows pinned the keyboard stops at table row 1000 - N (998 with two pinned) although rows 999 and
  1000 are drawn - measured on dev 2026-09-11 with 0, 1 and 2 pinned rows, reported to the operator
  as a suspected product defect rather than written down here as a known failure. Nor how many rows
  a layout pins back when one was pinned by a non-unique value: the warning says the value "won't
  be applied from the layout", yet the layout pins the first row holding it (row 1, the first F,
  in place of row 8) - measured on dev with the layout kept in memory and saved to the server
  alike, and put to the operator; the scenario claims only what both readings agree on, that row 8
  does not come back and row 5, pinned by its unique USUBJID, does. After the layout is loaded the
  grid is scrolled to its end with `user makes the last row current`: a click and a Control+End right
  after the load left the grid unscrolled on dev (reported to the operator as a possible product
  symptom), and which gesture scrolls is not the subject.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then grid should show 1000 rows
    And "Frozen Columns" property of grid should be "1"
    And the "pinned rows" reading of grid should be 0

  Scenario: Pin Column freezes one more column and Pin Row pins two unique rows
    Given user listens for "d4-grid-pinned_rows-changed" event on grid
    When user picks "Pin > Pin Column" from the context menu of the "header SEX" area of grid
    Then "Frozen Columns" property of grid should be "2"
    And the "column order" reading of grid should be "SEX, USUBJID, AGE, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    When user picks "Pin > Pin Row" from the context menu of the "cell 1 of USUBJID" area of grid
    Then "d4-grid-pinned_rows-changed" event should have fired on grid
    And the "pinned rows" reading of grid should be 1
    When user picks "Pin > Pin Row" from the context menu of the "cell 3 of USUBJID" area of grid
    Then the "pinned rows" reading of grid should be 2
    And no error or warning balloon should have been shown
    When user picks "Pin > Unpin All Rows" from the context menu of the "cell 1 of USUBJID" area of grid
    Then the "pinned rows" reading of grid should be 0
    When user picks "Pin > Unpin All Columns" from the context menu of the "header SEX" area of grid
    Then "Frozen Columns" property of grid should be "1"
    And no errors should have been logged

  Scenario: Control+clicks under two pinned rows select exactly the rows clicked
    When user picks "Pin > Pin Row" from the context menu of the "cell 1 of USUBJID" area of grid
    And user picks "Pin > Pin Row" from the context menu of the "cell 2 of USUBJID" area of grid
    Then the "pinned rows" reading of grid should be 2
    When user clicks on the "cell 5 of AGE" area of grid
    And user presses Control+End
    Then grid should have a "row header 1000" area
    And grid should have a "row header 1" area
    And grid should have a "row header 2" area
    And grid should not have a "row header 5" area
    When user clicks on the "row header 990" area of grid holding Control
    And user clicks on the "row header 995" area of grid holding Control
    And user clicks on the "row header 998" area of grid holding Control
    Then 3 rows should be selected
    And only rows where "USUBJID" is one of "X0273T51070200002, X0273T51080100021, X0273T51080200011" should be selected
    And no errors should have been logged

  Scenario: Unpin Row releases one pinned row and keeps the other
    When user picks "Pin > Pin Row" from the context menu of the "cell 3 of USUBJID" area of grid
    And user picks "Pin > Pin Row" from the context menu of the "cell 6 of USUBJID" area of grid
    Then the "pinned rows" reading of grid should be 2
    When user picks "Pin > Unpin Row" from the context menu of the "cell 3 of USUBJID" area of grid
    Then the "pinned rows" reading of grid should be 1
    When user clicks on the "cell 8 of AGE" area of grid
    And user presses Control+End
    Then grid should have a "row header 1000" area
    And grid should have a "row header 6" area
    And grid should not have a "row header 3" area
    And no errors should have been logged

  Scenario: Pin Selected Rows pins every selected row at once
    When user clicks on the "row header 5" area of grid
    And user clicks on the "row header 7" area of grid holding Shift
    Then rows 5 to 7 should be selected
    And 3 rows should be selected
    When user picks "Pin > Pin Selected Rows" from the context menu of the "cell 6 of USUBJID" area of grid
    Then the "pinned rows" reading of grid should be 3
    When user clicks on the "cell 9 of AGE" area of grid
    And user presses Control+End
    Then grid should have a "row header 1000" area
    And grid should have a "row header 5" area
    And grid should have a "row header 6" area
    And grid should have a "row header 7" area
    And no errors should have been logged

  Scenario: A row pinned by a non-unique value warns and is not the row the layout brings back
    When user picks "Pin > Pin Row" from the context menu of the "cell 8 of SEX" area of grid
    Then a warning balloon containing "non-unique" should have been shown
    And the "pinned rows" reading of grid should be 1
    When user picks "Pin > Pin Row" from the context menu of the "cell 5 of USUBJID" area of grid
    Then the "pinned rows" reading of grid should be 2
    When user saves the layout of the current table view
    And user picks "Pin > Unpin All Rows" from the context menu of the "cell 5 of USUBJID" area of grid
    Then the "pinned rows" reading of grid should be 0
    When user loads the saved layout
    And user makes the last row current
    Then grid should have a "row header 1000" area
    And grid should have a "row header 5" area
    And grid should not have a "row header 8" area
    And no errors should have been logged

  Scenario: Pin 2 Columns pins a column selection and Unpin Column releases one of them
    When user clicks on the "header RACE" area of grid holding Control
    And user clicks on the "header DIS_POP" area of grid holding Control
    Then columns "RACE, DIS_POP" should be selected
    When user picks "Pin > Pin 2 Columns" from the context menu of the "header RACE" area of grid
    Then "Frozen Columns" property of grid should be "3"
    And the "column order" reading of grid should be "RACE, DIS_POP, USUBJID, AGE, SEX, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    When user presses Escape
    And user picks "Pin > Unpin Column" from the context menu of the "header RACE" area of grid
    Then "Frozen Columns" property of grid should be "2"
    And the "column order" reading of grid should be "DIS_POP, RACE, USUBJID, AGE, SEX, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    When user picks "Pin > Unpin All Columns" from the context menu of the "header DIS_POP" area of grid
    Then "Frozen Columns" property of grid should be "1"
    And no errors should have been logged

  Scenario: A header dropped inside the pinned columns is pinned, one dropped outside is only moved
    When user picks "Pin > Pin Column" from the context menu of the "header SEX" area of grid
    And user picks "Pin > Pin Column" from the context menu of the "header RACE" area of grid
    Then "Frozen Columns" property of grid should be "3"
    When user drags the "header DEMOG" area of grid to the "header SEX" area
    Then "Frozen Columns" property of grid should be "4"
    And the "column order" reading of grid should be "SEX, DEMOG, RACE, USUBJID, AGE, DIS_POP, HEIGHT, WEIGHT, CONTROL, STARTED, SEVERITY"
    When user drags the "header HEIGHT" area of grid to the "header AGE" area
    Then "Frozen Columns" property of grid should be "4"
    And the "column order" reading of grid should be "SEX, DEMOG, RACE, USUBJID, AGE, HEIGHT, DIS_POP, WEIGHT, CONTROL, STARTED, SEVERITY"
    And no errors should have been logged

  Scenario: The arrow keys cross the pinned-column boundary both ways
    When user picks "Pin > Pin Column" from the context menu of the "header SEX" area of grid
    Then "Frozen Columns" property of grid should be "2"
    When user clicks on the "cell 3 of SEX" area of grid
    Then the "current column" reading of grid should be "SEX"
    When user presses ArrowRight
    Then the "current column" reading of grid should be "USUBJID"
    When user presses ArrowRight
    Then the "current column" reading of grid should be "AGE"
    When user presses ArrowLeft
    And user presses ArrowLeft
    Then the "current column" reading of grid should be "SEX"
    And row 3 should be current
    And no errors should have been logged
