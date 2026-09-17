@journey @viewers @realizes:viewers.grid
Feature: Grid row selection, navigation and row source
  How the grid picks rows out: a click sets the current row without selecting, a drag down the
  row-number strip selects a range with Shift, Control+A takes everything and Escape drops rows, columns and
  the current row together, Control+click on headers selects columns, the arrows and Control+Home /
  Control+End / PageDown move the current row, Space and Shift+Enter select by value, Allow Row
  Selection gates all of it, and Row Source decides how many rows the grid shows at all. One journey
  on demog-1000 (SEX M 447, AGE above 30 844).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then grid should show 1000 rows
    And no rows should be selected

  Scenario: A click sets the current row and selects nothing
    When user clicks on the "cell 6 of USUBJID" area of grid
    Then row 6 should be current
    And the "current row" reading of grid should be 6
    And the "current column" reading of grid should be "USUBJID"
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
    When user presses Escape
    Then no rows should be selected
    And no columns should be selected
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
    When user picks "Sort > Reset" from the context menu of the "header AGE" area of grid
    Then the "row order" reading of grid should be as remembered
    When user sets "Row Source" property of grid to "All"
    And user resets the filter
    Then all rows should pass the filter
    And grid should show 1000 rows
    And no errors should have been logged
