@journey
Feature: The table view's keyboard shortcuts toggle the shell's panels, and Alt+C's Column Manager acts on the grid
  The shortcuts of `help/datagrok/navigation/shortcuts.md` (and the Console and Variables sections of
  `navigation.md`) that no other feature presses, each claimed by what the user sees open or change: Alt+X
  hides and shows the Toolbox, Alt+V docks and closes the Variables panel, the tilde opens the Console, Alt+S
  opens Settings, and Ctrl+Shift+V opens a second view of the same table. Alt+C toggles the Columns pane
  (Column Manager, `help/datagrok/navigation/panels/column-manager.md`): a grid with a row per column of the
  current table (its type letter and its name) and a search box; a click on a row makes that column current,
  Ctrl+click toggles the column's selection in the grid and Esc clears it, dragging a row reorders the
  columns of the table view's grid, Hide and Show in its popup menu take the column out of the grid and
  bring it back, and Statistics adds a statistics column to the list.

  Left out on purpose, one line each: F4, F7, Alt+T, Ctrl+O, Ctrl+H and Ctrl+Z are pressed by their own
  features; Ctrl+Shift+R is `undo-redo.feature`'s; Ctrl+S saves the layout the way the Layouts pane's Save does
  (`layouts/layout-names-tags.feature`); Ctrl+L docks two panels no step can close yet (MISSING.md); F1, F11
  and Alt+F need steps that do not exist (MISSING.md); the Column Manager's type filter and the pane following
  another table are `viewers/grid/grid-columns.feature`'s (GROK-19332); its other popup commands (Remove,
  Extract, Change Type..., Normalize...) are the column's own, claimed from the grid's header menu.

  One journey on a four-column table written here; nothing goes to the server. Every scenario closes what it
  opened and puts back the selection, order and visibility it changed.

  Background:
    Given user is logged in
    And user opens a table "bdd-colman" with:
      | name | num | price | date       |
      | a    | 1   | 1.5   | 2020-01-01 |
      | b    | 2   | 2.5   | 2020-01-02 |
      | c    | 3   | 3.5   | 2020-01-03 |

  Scenario: Alt+X, Alt+V and the tilde toggle the Toolbox, the Variables panel and the Console
    Given the toolbox pane is shown
    Then toolbox should be visible
    When user presses Alt+X
    Then toolbox should be hidden
    And "Toolbox" status bar toggle should not be selected
    When user presses Alt+X
    Then toolbox should be visible
    And "Toolbox" status bar toggle should be selected
    And "Variables" dock panel should be absent
    When user presses Alt+V
    Then "Variables" dock panel should be visible
    And "Variables" status bar toggle should be selected
    When user presses Alt+V
    Then "Variables" dock panel should be absent
    And console should be absent
    When user presses Backquote
    Then console should be visible
    And "Console" status bar toggle should be selected
    When user closes the console
    Then console should be hidden
    And no errors should have been logged

  Scenario: Alt+S opens the Settings
    When user presses Alt+S
    Then the current view should be a settings view
    And the "Settings - General" view should be current
    When user closes the current view
    Then the "bdd-colman" view should be current
    And no errors should have been logged

  Scenario: Ctrl+Shift+V opens a second view on the same table
    When user presses Control+Shift+V
    Then the "bdd-colman copy" view should be current
    And the open table views should be exactly "bdd-colman, bdd-colman"
    When user closes the current view
    Then the open table views should be exactly "bdd-colman"
    And no errors should have been logged

  Scenario: Alt+C docks the Columns pane with a row per column, a click on a row makes that column current, and Alt+C closes it
    Then column manager should be absent
    When user clicks on the "cell 1 of name" area of grid
    Then the current column should be "name"
    When user presses Alt+C
    Then column manager should be visible
    And the "rows shown" reading of Grid viewer in column manager should be 4
    And the "text of cell 1 of __name" reading of Grid viewer in column manager should be "name"
    And the "text of cell 4 of __name" reading of Grid viewer in column manager should be "date"
    And the "text of cell 1 of __t" reading of Grid viewer in column manager should be "s"
    And the "text of cell 2 of __t" reading of Grid viewer in column manager should be "i"
    And the "text of cell 3 of __t" reading of Grid viewer in column manager should be "d"
    When user clicks on the "cell 3 of __name" area of Grid viewer in column manager
    Then the current column should be "price"
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
    When user clears "Search" input in column manager
    Then the "rows shown" reading of Grid viewer in column manager should be 4
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
    And user presses Alt+C
    Then column manager should be absent
    And the "column order" reading of grid should be "name, price, date"
    When user presses Alt+C
    Then the "rows shown" reading of Grid viewer in column manager should be 4
    When user picks "'num' column > Show" from the context menu of the "cell 2 of __name" area of Grid viewer in column manager
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
