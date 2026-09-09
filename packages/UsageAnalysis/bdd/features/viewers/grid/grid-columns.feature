@journey @viewers @realizes:viewers.grid
Feature: Grid column geometry
  Everything that decides where a column sits and how wide it is: the Column Sizing presets,
  sorting by a double-click on the header (the first one sorts descending), resizing a column and
  the rows by dragging their handles, reordering a column onto another, the Order or Hide Columns
  dialog, and pinning a column and two rows. The grid reports its own geometry — `column width of
  <col>`, `column order`, `pinned rows`, `sort column`, `sort direction` — and the drag handles as
  the hit areas `column resizer <col>` and `row resizer <row>`, so nothing is measured from private
  bounds. One journey on demog-1000, whose table order the sorting never touches.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then grid should show 1000 rows
    And "Frozen Columns" property of grid should be "1"
    And the value of "AGE" column in row 1 should be "26"

  Scenario: Column Sizing presets order the widths Minimal, Optimal, Maximal
    When user picks "Column Sizing > Optimal" from the context menu of the "cell 4 of AGE" area of grid
    Then grid should have repainted
    When user picks "Column Sizing > Minimal" from the context menu of the "cell 4 of AGE" area of grid
    Then the "column width of AGE" reading of grid should be lower than before
    And the "column width of RACE" reading of grid should be lower than before
    When user picks "Column Sizing > Maximal" from the context menu of the "cell 4 of AGE" area of grid
    Then the "column width of AGE" reading of grid should be higher than before
    And the "column width of SEVERITY" reading of grid should be higher than before
    When user picks "Column Sizing > Optimal" from the context menu of the "cell 4 of AGE" area of grid
    Then the "column width of SEVERITY" reading of grid should be lower than before
    And no errors should have been logged

  Scenario: A double-click on the header sorts descending, then ascending, then off
    Given user listens for "d4-grid-rows-sorted" event on grid
    When user double-clicks on the "header AGE" area of grid
    Then "d4-grid-rows-sorted" event should have fired on grid
    And the "sort column" reading of grid should be "AGE"
    And the "sort direction" reading of grid should be "descending"
    And grid should have a "cell 600 of AGE" area
    And the value of "AGE" column in row 1 should be "26"
    When user double-clicks on the "header AGE" area of grid
    Then the "sort direction" reading of grid should be "ascending"
    And grid should have a "cell 692 of AGE" area
    And grid should not have a "cell 600 of AGE" area
    When user double-clicks on the "header AGE" area of grid
    Then the "sort column" reading of grid should be ""
    And grid should have a "cell 1 of AGE" area
    And grid should not have a "cell 692 of AGE" area
    And the value of "AGE" column in row 1 should be "26"
    And no errors should have been logged

  Scenario: The Sort dialog sorts by its first level
    When user picks "Sort..." from the context menu of the "cell 4 of AGE" area of grid
    Then Sort Table dialog should be visible
    When user clicks on OK button in Sort Table dialog
    Then Sort Table dialog should be hidden
    And the "sort column" reading of grid should be "USUBJID"
    And the "sort direction" reading of grid should be "ascending"
    And the value of "AGE" column in row 1 should be "26"
    When user picks "Sort > Reset" from the context menu of the "header AGE" area of grid
    Then the "sort column" reading of grid should be ""
    And no errors should have been logged

  Scenario: Dragging the header handle resizes one column and the row handle every row
    When user drags the "column resizer AGE" area of grid by 40 pixels to the right
    Then the "column width of AGE" reading of grid should be higher than before
    And the "column width of SEX" reading of grid should be the same as before
    When user drags the "row resizer 1" area of grid by 12 pixels to the down
    Then the "cell 1 of AGE" area of grid should be taller than before
    And the "column width of AGE" reading of grid should be the same as before
    And grid should have repainted
    When user sets "Row Height" property of grid to "28"
    And user picks "Column Sizing > Optimal" from the context menu of the "cell 4 of AGE" area of grid
    Then no errors should have been logged

  Scenario: Dragging a header onto another reorders the columns
    When user drags the "header HEIGHT" area of grid to the "header DEMOG" area
    Then the "column order" reading of grid should differ from before
    And grid should have a "header HEIGHT" area
    And grid should have a "header DEMOG" area
    And the table should have 11 columns
    And no errors should have been logged
    When user drags the "header HEIGHT" area of grid to the "header DIS_POP" area

  Scenario: Order or Hide Columns opens on the grid and keeps its title
    When user picks "Order or Hide Columns..." from the context menu of the "cell 4 of AGE" area of grid
    Then Order or Hide Columns dialog should be visible
    And title of Order or Hide Columns dialog should contain text "Order or Hide Columns"
    When user closes Order or Hide Columns dialog
    Then Order or Hide Columns dialog should be hidden
    And no errors should have been logged

  Scenario: Pin Column freezes one more column and Pin Row pins two unique rows
    Given user listens for "d4-grid-pinned_rows-changed" event on grid
    When user picks "Pin > Pin Column" from the context menu of the "header SEX" area of grid
    Then "Frozen Columns" property of grid should be "2"
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
