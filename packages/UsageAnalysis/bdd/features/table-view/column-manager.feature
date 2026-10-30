@journey
Feature: Column Manager lists the table's columns and acts on the table view's grid
  Alt+C toggles the Columns pane (Column Manager), as `help/datagrok/navigation/shortcuts.md` and
  `help/datagrok/navigation/panels/column-manager.md` say: a grid with a row per column of the current table
  (its type letter and its name) and a search box. Its usage table is claimed row by row: a click jumps to the
  column (the table's current column), Ctrl+click toggles the column's selection in the grid and Esc clears
  it, dragging a row reorders the columns of the table view's grid, and the popup menu acts on the column —
  Hide takes it out of the grid and Show brings it back — and Statistics adds a statistics column to the list.

  Left out on purpose, one line each: the type filter and the pane following another table are
  `viewers/grid/grid-columns.feature`'s (GROK-19332); Shift+drag is the grid's own row selection, pressed by
  the grid features; applying a command to "all selected tables" needs several tables selected in the Tables
  pane, which no step selects (MISSING.md); the other commands of the popup (Remove, Extract, Change Type...,
  Normalize...) are the column's own commands, claimed by their features from the grid's header menu.

  One journey on a four-column table written here; nothing goes to the server. Every scenario closes the pane
  and puts back the selection, order and visibility it changed.

  Background:
    Given user is logged in
    And user opens a table "bdd-colman" with:
      | name | num | price | date       |
      | a    | 1   | 1.5   | 2020-01-01 |
      | b    | 2   | 2.5   | 2020-01-02 |
      | c    | 3   | 3.5   | 2020-01-03 |

  Scenario: Alt+C docks the Columns pane with a row per column, and Alt+C closes it
    Then column manager should be absent
    When user clicks on the "cell 1 of name" area of grid
    And user presses Alt+C
    Then column manager should be visible
    And the "rows shown" reading of Grid viewer in column manager should be 4
    And the "text of cell 1 of __name" reading of Grid viewer in column manager should be "name"
    And the "text of cell 4 of __name" reading of Grid viewer in column manager should be "date"
    And the "text of cell 1 of __t" reading of Grid viewer in column manager should be "s"
    And the "text of cell 2 of __t" reading of Grid viewer in column manager should be "i"
    And the "text of cell 3 of __t" reading of Grid viewer in column manager should be "d"
    When user presses Alt+C
    Then column manager should be absent
    And no errors should have been logged

  Scenario: The search box narrows the list to the matching columns
    When user clicks on the "cell 1 of name" area of grid
    And user presses Alt+C
    Then the "rows shown" reading of Grid viewer in column manager should be 4
    When user enters "pr" into "Search" input in column manager
    Then the "rows shown" reading of Grid viewer in column manager should be 1
    And the "text of cell 3 of __name" reading of Grid viewer in column manager should be "price"
    When user enters "a" into "Search" input in column manager
    Then the "rows shown" reading of Grid viewer in column manager should be 2
    When user presses Alt+C
    Then column manager should be absent
    And no errors should have been logged

  Scenario: A click on a column's row makes it the current column
    When user clicks on the "cell 1 of name" area of grid
    Then the current column should be "name"
    When user presses Alt+C
    And user clicks on the "cell 3 of __name" area of Grid viewer in column manager
    Then the current column should be "price"
    When user presses Alt+C
    Then column manager should be absent
    And no errors should have been logged

  Scenario: Ctrl+click toggles a column's selection in the grid, and Esc clears it
    When user clicks on the "cell 1 of name" area of grid
    And user presses Alt+C
    Then no columns should be selected
    When user clicks on the "cell 2 of __name" area of Grid viewer in column manager holding Control
    And user clicks on the "cell 4 of __name" area of Grid viewer in column manager holding Control
    Then columns "num, date" should be selected
    When user clicks on the "cell 4 of __name" area of Grid viewer in column manager holding Control
    Then columns "num" should be selected
    When user presses Escape
    Then no columns should be selected
    When user presses Alt+C
    Then column manager should be absent
    And no errors should have been logged

  Scenario: Dragging a column's row reorders the columns of the grid
    Then the "column order" reading of grid should be "name, num, price, date"
    When user clicks on the "cell 1 of name" area of grid
    And user presses Alt+C
    And user drags the "cell 4 of __name" area of Grid viewer in column manager to the "cell 1 of __name" area
    And user presses Alt+C
    Then column manager should be absent
    And the "column order" reading of grid should be "date, name, num, price"
    When user presses Alt+C
    Then the "text of cell 1 of __name" reading of Grid viewer in column manager should be "date"
    When user drags the "cell 1 of __name" area of Grid viewer in column manager to the "cell 4 of __name" area
    And user presses Alt+C
    Then column manager should be absent
    And the "column order" reading of grid should be "name, num, price, date"
    And no errors should have been logged

  Scenario: Hide in the popup menu takes the column out of the grid, and Show brings it back
    Then the "column order" reading of grid should be "name, num, price, date"
    When user clicks on the "cell 1 of name" area of grid
    And user presses Alt+C
    And user picks "'num' column > Hide" from the context menu of the "cell 2 of __name" area of Grid viewer in column manager
    Then the "rows shown" reading of Grid viewer in column manager should be 4
    When user presses Alt+C
    Then column manager should be absent
    And the "column order" reading of grid should be "name, price, date"
    When user presses Alt+C
    And user picks "'num' column > Show" from the context menu of the "cell 2 of __name" area of Grid viewer in column manager
    And user presses Alt+C
    Then column manager should be absent
    And the "column order" reading of grid should be "name, num, price, date"
    And no errors should have been logged

  Scenario: Statistics in the popup menu adds a statistics column to the list
    When user clicks on the "cell 1 of name" area of grid
    And user presses Alt+C
    Then the "column order" reading of Grid viewer in column manager should be "__t, __name"
    When user picks "Statistics > count" from the context menu of the "cell 1 of __name" area of Grid viewer in column manager
    Then the "column order" reading of Grid viewer in column manager should be "__t, __name, count"
    And the "text of cell 2 of count" reading of Grid viewer in column manager should be "3.00"
    When user picks "Statistics > count" from the context menu of the "cell 1 of __name" area of Grid viewer in column manager
    Then the "column order" reading of Grid viewer in column manager should be "__t, __name"
    When user presses Alt+C
    Then column manager should be absent
    And no errors should have been logged
