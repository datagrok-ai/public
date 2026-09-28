@journey @viewers @realizes:viewers.grid
Feature: Grid row selection, navigation and row source
  How the grid picks rows out: a click sets the current row without selecting, a drag down the
  row-number strip selects a range with Shift, Control+A takes everything and Escape drops rows,
  columns and the current row together, Control+click on headers selects columns, the arrows and
  Control+Home / Control+End / PageDown move the current row, Tab and Shift+Tab walk the cells and
  wrap to the next and the previous row at the edges, Space and Shift+Enter select by value, Allow
  Row Selection gates all of it (Shift+ArrowDown moves the caret and selects nothing), and Row
  Source decides how many rows the grid shows at all. Under a filter and a sort the arrows walk the
  grid's own order and Shift+ArrowDown selects the rows it leaves - the AGE-31 rows come first, in
  row order, since the sort keeps tied rows as they were: 95, 115, 117, 210, 320, 430, 497. Also the mouse gestures of `grid-ui.md`: a Shift+click
  range on the row strip; a Control+click on a row number inverts that row's selection without
  making it current, and a plain click moves the current row and leaves the selection as it was; a
  drag down the strip; Shift+click and Control+click on the headers. The Tab wrap, the plain click
  that keeps the selection and the Control+click that does not move the current row are the
  product's behaviour as the operator confirmed it, where `grid-ui.md` and
  `grid-rows-select-filter-navigate.md` scenario 1 expect otherwise. One journey on demog-1000 (SEX
  M 447, AGE above 30 844).

  Escape is pressed with the focus given to the grid's overlay canvas and the focus claimed first:
  the grid clears the current row only for a key that reaches that canvas, and clears the selection
  alone for one that reaches its root, so without the focus step the claim depended on where the
  click had left it (measured on dev). Rows pinned under a scroll are `grid-pinning.feature`'s;
  Control+F and the block drag across the cells are out of scope by operator decision.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then grid should show 1000 rows
    And no rows should be selected

  Scenario: A click sets the current row and selects nothing
    When user clicks on the "cell 6 of AGE" area of grid
    Then row 6 should be current
    And the "current row" reading of grid should be 6
    And the "current column" reading of grid should be "AGE"
    And no rows should be selected
    And no errors should have been logged

  Scenario: Shift-dragging the row-number strip selects a range
    When user drags a selection box from the "row header 6" area to the "row header 11" area of grid
    Then 6 rows should be selected
    And rows 6 to 11 should be selected
    And no columns should be selected
    And grid should have repainted
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged

  Scenario: Control+A selects every row and every column and Escape clears both
    When user clicks on the "cell 1 of AGE" area of grid
    And user presses Control+A
    Then all rows should be selected
    And columns "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY" should be selected
    When user focuses on grid overlay
    Then grid overlay should be focused
    When user presses Escape
    Then no rows should be selected
    And no columns should be selected
    And the "current row" reading of grid should be 0
    And no errors should have been logged

  Scenario: Six Control+clicks on headers select six columns without an error
    When user clicks on the "header AGE" area of grid holding Control
    And user clicks on the "header HEIGHT" area of grid holding Control
    And user clicks on the "header WEIGHT" area of grid holding Control
    And user clicks on the "header SEX" area of grid holding Control
    And user clicks on the "header RACE" area of grid holding Control
    And user clicks on the "header DIS_POP" area of grid holding Control
    Then columns "AGE, HEIGHT, WEIGHT, SEX, RACE, DIS_POP" should be selected
    And no errors should have been logged
    When user presses Escape
    Then no columns should be selected

  Scenario: Arrows and Control+Home, Control+End and PageDown move the current row
    When user clicks on the "cell 5 of AGE" area of grid
    And user presses ArrowDown
    And user presses ArrowDown
    And user presses ArrowRight
    And user presses ArrowRight
    Then row 7 should be current
    And the "current column" reading of grid should be "RACE"
    When user presses Control+Home
    Then row 1 should be current
    When user presses Control+End
    Then row 1000 should be current
    When user presses Control+Home
    And user takes a snapshot of grid
    And user presses PageDown
    Then the "current row" reading of grid should be higher than before
    When user takes a snapshot of grid
    And user presses PageDown
    Then the "current row" reading of grid should be higher than before
    When user presses Control+Home
    Then row 1 should be current
    And no errors should have been logged

  Scenario: Space selects the current row and Shift+Enter every row of the same value
    When user clicks on the "cell 4 of AGE" area of grid
    And user presses Space
    Then 1 row should be selected
    And only rows where "USUBJID" is "X0273T21000400002" should be selected
    When user presses Space
    Then no rows should be selected
    When user clicks on the "cell 4 of SEX" area of grid
    And user presses Shift+Enter
    Then only rows where "SEX" is "M" should be selected
    And 447 rows should be selected
    When user clears the row selection
    Then no errors should have been logged

  Scenario: Allow Row Selection off keeps the caret moving and the selection empty
    When user sets "Allow Row Selection" property of grid to "false"
    And user clicks on the "cell 5 of AGE" area of grid
    Then row 5 should be current
    When user presses ArrowDown
    Then row 6 should be current
    When user presses Shift+ArrowDown
    And user presses Shift+ArrowDown
    And user presses Shift+ArrowDown
    Then row 9 should be current
    And no rows should be selected
    When user presses Space
    Then no rows should be selected
    When user sets "Allow Row Selection" property of grid to "true"
    And user presses Space
    Then 1 row should be selected
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged

  Scenario: Row Source decides how many rows the grid shows
    When user filters rows where "SEX" is "M"
    Then 447 rows should pass the filter
    When user selects the first 5 rows
    Then 5 rows should be selected
    When user sets "Row Source" property of grid to "Filtered"
    Then grid should show 447 rows
    When user sets "Row Source" property of grid to "Selected"
    Then grid should show 5 rows
    When user sets "Row Source" property of grid to "All"
    Then grid should show 1000 rows
    And 447 rows should pass the filter
    When user resets the filter
    And user clears the row selection
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: A filter and a sort share one order
    When user filters rows where "AGE" is between 31 and 999
    Then 844 rows should pass the filter
    When user sets "Row Source" property of grid to "Filtered"
    Then grid should show 844 rows
    When user remembers the "row order" reading of grid
    And user double-clicks on the "header AGE" area of grid
    And user double-clicks on the "header AGE" area of grid
    Then the "sort direction" reading of grid should be "ascending"
    And grid should show 844 rows
    And the "row order" reading of grid should differ from before
    When user presses Control+Home
    Then the "current row" reading of grid should be 95
    And "AGE" of the current row should be "31"
    When user presses ArrowDown
    Then the "current row" reading of grid should be 115
    When user presses ArrowDown
    Then the "current row" reading of grid should be 117
    When user presses Shift+ArrowDown
    And user presses Shift+ArrowDown
    And user presses Shift+ArrowDown
    And user presses Shift+ArrowDown
    Then the "current row" reading of grid should be 497
    And only rows where "USUBJID" is one of "X0273T22000400018, X0273T26004000004, X0273T29013500101, X0273T37001500013" should be selected
    And every selected row should pass the filter
    And 844 rows should pass the filter
    When user clears the row selection
    And user picks "Sort > Reset" from the context menu of the "header AGE" area of grid
    Then the "row order" reading of grid should be as remembered
    When user sets "Row Source" property of grid to "All"
    And user resets the filter
    Then all rows should pass the filter
    And grid should show 1000 rows
    And no errors should have been logged

  Scenario: Shift+click, Control+click, a plain click and a drag on the row strip
    When user clicks on the "row header 5" area of grid
    And user clicks on the "row header 10" area of grid holding Shift
    Then rows 5 to 10 should be selected
    And 6 rows should be selected
    And the "current row" reading of grid should be 10
    When user clicks on the "row header 15" area of grid holding Control
    Then 7 rows should be selected
    And only rows where "USUBJID" is one of "X0273T21000500006, X0273T21000500008, X0273T21000700006, X0273T21000800002, X0273T21000800004, X0273T21000800006, X0273T21001100003" should be selected
    And the "current row" reading of grid should be 10
    When user clicks on the "row header 15" area of grid holding Control
    Then 6 rows should be selected
    And rows 5 to 10 should be selected
    And the "current row" reading of grid should be 10
    When user clicks on the "row header 15" area of grid holding Control
    Then 7 rows should be selected
    When user clicks on the "row header 8" area of grid
    Then 7 rows should be selected
    And only rows where "USUBJID" is one of "X0273T21000500006, X0273T21000500008, X0273T21000700006, X0273T21000800002, X0273T21000800004, X0273T21000800006, X0273T21001100003" should be selected
    And the "current row" reading of grid should be 8
    When user drags the "row header 12" area of grid to the "row header 16" area
    Then rows 12 to 16 should be selected
    And 5 rows should be selected
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged

  Scenario: Shift+click adds headers to the selection and Control+click takes one back
    When user clicks on the "header AGE" area of grid holding Shift
    And user clicks on the "header HEIGHT" area of grid holding Shift
    Then columns "AGE, HEIGHT" should be selected
    When user clicks on the "header AGE" area of grid holding Control
    Then columns "HEIGHT" should be selected
    When user presses Escape
    Then no columns should be selected
    And no errors should have been logged

  Scenario: Tab and Shift+Tab walk the cells and wrap at the row edges
    When user clicks on the "cell 3 of AGE" area of grid
    And user presses Tab
    Then the "current column" reading of grid should be "SEX"
    And the "current row" reading of grid should be 3
    When user presses Shift+Tab
    Then the "current column" reading of grid should be "AGE"
    And the "current row" reading of grid should be 3
    When user clicks on the "cell 3 of SEVERITY" area of grid
    And user presses Tab
    Then the "current column" reading of grid should be "USUBJID"
    And the "current row" reading of grid should be 4
    When user presses Shift+Tab
    Then the "current column" reading of grid should be "SEVERITY"
    And the "current row" reading of grid should be 3
    When user clicks on the "cell 3 of USUBJID" area of grid
    And user presses Shift+Tab
    Then the "current column" reading of grid should be "SEVERITY"
    And the "current row" reading of grid should be 2
    When user presses Tab
    Then the "current column" reading of grid should be "USUBJID"
    And the "current row" reading of grid should be 3
    And no errors should have been logged
