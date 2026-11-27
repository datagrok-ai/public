@guide @tutorial
Feature: How-to · Basics
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). The grid's everyday gestures on demog: the keys that move the
  current cell (arrows, Ctrl+End, Ctrl+Home), sorting by a double-click on a header (descending, ascending, reset),
  a column widened by its header's edge, a column moved by dragging its header, a column hidden from the header's
  menu and shown again from Order or Hide Columns, and a taller row height dragged between row numbers. Then (it
  absorbed the former "Find a column in a wide table" clip) a click on the spgi-100 tab switches to a table wider than
  the screen (88 columns): Pin > Pin Column from the Structure header's menu keeps the structures in sight,
  "Columns: 88" in the status bar docks the Columns pane, its search narrows the list to "hERG Classification", and a
  click on it takes the grid there, the pinned Structure column still beside it. spgi and demog are opened off camera
  (demog in front, spgi-100's tab beside it from the start); nothing is made on the server.

  Scenario: Work with the grid
    Move around, sort, resize, reorder, hide and show columns, change the row height, and find a column in a wide table.

    Given user is logged in
    And simple mode is off
    And user opens spgi dataset
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    # caption: Click a cell; the arrow keys move from there
    When user clicks on the "cell 1 of AGE" area of grid
    And user presses ArrowDown
    And user presses ArrowDown
    And user presses ArrowRight
    Then row 3 should be current
    # caption: Ctrl+End jumps to the last row, Ctrl+Home back to the first
    When user presses Control+End
    Then row 5850 should be current
    When user presses Control+Home
    Then row 1 should be current
    # caption: Double-click a header to sort by it
    When user double-clicks on the "header AGE" area of grid
    Then the "sort column" reading of grid should be "AGE"
    And the "sort direction" reading of grid should be "descending"
    # caption: Again for the other direction, a third time to unsort
    When user double-clicks on the "header AGE" area of grid
    Then the "sort direction" reading of grid should be "ascending"
    When user double-clicks on the "header AGE" area of grid
    Then the "sort column" reading of grid should be ""
    # caption: Drag a header's right edge to resize the column
    When user drags the "column resizer RACE" area of grid by 60 pixels to the right
    Then the "column width of RACE" reading of grid should be higher than before
    # caption: Drag a header to move the column
    When user drags the "header HEIGHT" area of grid to the "header AGE" area
    Then the "column order" reading of grid should include the text "AGE, HEIGHT, SEX"
    # caption: Right-click a header to hide its column
    When user picks "Hide" from the context menu of the "header WEIGHT" area of grid
    Then grid should not have a "header WEIGHT" area
    # caption: Order or Hide Columns lists them all: tick it to show it again
    When user picks "Order or Hide Columns..." from the context menu of the "header AGE" area of grid
    Then Order or Hide Columns dialog should be visible
    When user toggles the "WEIGHT" column in the column list of Order or Hide Columns dialog
    And user clicks on CLOSE button in Order or Hide Columns dialog
    Then grid should have a "header WEIGHT" area
    # caption: Drag between row numbers to make rows taller
    When user drags the "row resizer 3" area of grid by 20 pixels to the down
    Then the "row header 3" area of grid should be taller than before
    # caption: Now a wider table: click its tab
    When user clicks on the tab of the "spgi-100" view
    Then the "spgi-100" view should be current
    # caption: 88 columns: more than fit on the screen
    When user hovers over "Columns: 88" text in status bar
    # caption: Pin a column to keep it in sight
    And user picks "Pin > Pin Column" from the context menu of the "header Structure" area of grid
    Then "Frozen Columns" property of grid should be "2"
    # caption: Columns in the status bar lists them all
    When user clicks on "Columns: 88" text in status bar
    Then column manager should be visible
    # caption: Type to find one by name
    When user types "herg" into "Search" input in column manager
    Then the "rows shown" reading of Grid viewer in column manager should be 1
    # caption: Click it: the grid takes you there
    When user clicks on the "cell 36 of __name" area of Grid viewer in column manager
    Then the current column should be "hERG Classification"
    And no errors should have been logged
