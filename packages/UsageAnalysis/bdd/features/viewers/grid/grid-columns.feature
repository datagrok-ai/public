@journey @viewers @realizes:viewers.grid
Feature: Grid column geometry
  Everything that decides where a column sits and how wide it is: the Column Sizing presets,
  sorting by a double-click on the header (the first one sorts descending) and by two levels in the
  Sort dialog, resizing a column and the rows by dragging their handles - two selected columns
  resize together until Escape breaks the link, and a column dragged past its left edge collapses
  to a hairline and stays in the column order - reordering a column onto another, widening columns
  until the grid scrolls sideways (GROK-19753), and the Order or Hide Columns dialog: its type
  filter, Reset filter (GROK-19333) and its title after both (GROK-20167), and the column manager
  the status bar docks - a click on "Columns: 11" - keeping its type filter while it rebinds from
  one table to another and back (GROK-19332). Pinning lives in
  `grid-pinning.feature`. The grid reports its own geometry - `column width of <col>`, `column
  order`, `sort column`, `sort direction`, `x scroll span` - and the drag handles as the hit areas
  `column resizer <col>`, `row resizer <row>` and `x scroll handle`, so nothing is measured from
  private bounds; the dialog lists the columns in a grid of its own, whose `rows shown` is how many
  columns the filter leaves. One journey on demog-1000, whose table order the sorting never touches
  (F is oldest at 80 in row 720, M youngest at 18 in row 816; the Sort dialog's second level is
  descending by default, and the feature relies on that default rather than setting it). The
  sideways-scroll scenario also resizes a column while the grid is scrolled, the case GROK-19753
  was about.

  The second table is the first 100 rows of demog-1000 with RACE, DIS_POP, DEMOG and SEVERITY
  removed through the API (data preparation), so it holds two string columns where the first holds
  six and the manager's count tells which table it lists. Not translated, and why: from `grid.md` "Header Histogram Strip" - the
  grid reports no reading for the kind of its header (Top > Histogram is a checked menu item, but a
  mark on the menu says nothing about the strip being drawn), and from "Multi-Column and Row-Height Resizing" the
  values of a squeezed column turning into small circles and the hairline column showing no
  content - `text of cell` stays the value and nothing reports how a cell is drawn; both need a
  reading in the core. The row-height drag stays as before: one row header dragged by 12 px, the
  cell taller than before. From
  `grid-dialogs-groups.md` scenario 3 and `grid-ui.md` "Context Panel - Column Hamburger Menu
  Inline Filter": the hamburger icon a header shows on hover is painted on the grid's canvas and is
  neither a hit area nor an element, so its popup, its inline filter, Add filter and the colour
  sync with the Context Panel (GROK-19288) need a `header menu icon <col>` area in the core.

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

  Scenario: The Sort dialog sorts by two levels
    When user picks "Sort..." from the context menu of the "cell 4 of AGE" area of grid
    Then Sort Table dialog should be visible
    When user selects "SEX" in first column input in Sort Table dialog
    And user selects "AGE" in second column input in Sort Table dialog
    And user clicks on OK button in Sort Table dialog
    Then Sort Table dialog should be hidden
    And the "sort column" reading of grid should be "SEX"
    And the "sort direction" reading of grid should be "ascending"
    And the value of "AGE" column in row 1 should be "26"
    When user clicks on the "cell 597 of AGE" area of grid
    Then row 597 should be current
    When user presses Control+Home
    Then row 720 should be current
    And "SEX" of the current row should be "F"
    And "AGE" of the current row should be "80"
    When user presses Control+End
    Then row 816 should be current
    And "SEX" of the current row should be "M"
    And "AGE" of the current row should be "18"
    When user presses Control+Home
    And user picks "Sort > Reset" from the context menu of the "header AGE" area of grid
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
    Then the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, WEIGHT, DEMOG, HEIGHT, CONTROL, STARTED, SEVERITY"
    And no errors should have been logged
    When user drags the "header HEIGHT" area of grid to the "header DIS_POP" area
    Then the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"

  Scenario: Order or Hide Columns filters its list by type and Reset filter clears it (GROK-19333, GROK-20167)
    When user picks "Order or Hide Columns..." from the context menu of the "cell 4 of AGE" area of grid
    Then Order or Hide Columns dialog should be visible
    And the "rows shown" reading of Grid viewer in Order or Hide Columns dialog should be 11
    When user hovers over "Search" input in Order or Hide Columns dialog
    And user clicks on font icon menu icon in Order or Hide Columns dialog
    And user hovers over "Types" menu item in context menu
    And user clicks on "int" menu item in context menu
    Then the "rows shown" reading of Grid viewer in Order or Hide Columns dialog should be 1
    And the "text of cell 2 of __name" reading of Grid viewer in Order or Hide Columns dialog should be "AGE"
    And "int" menu item in context menu should be selected
    When user clicks on "Reset filter" menu item in context menu
    Then the "rows shown" reading of Grid viewer in Order or Hide Columns dialog should be 11
    And "int" menu item in context menu should not be selected
    And title of Order or Hide Columns dialog should contain text "Order or Hide Columns"
    When user closes the context menu
    Then Order or Hide Columns dialog should be visible
    When user clicks on CLOSE button in Order or Hide Columns dialog
    Then Order or Hide Columns dialog should be hidden
    And no errors should have been logged

  Scenario: Two selected columns resize together until Escape, and a column collapses to a hairline
    When user clicks on the "header AGE" area of grid holding Control
    And user clicks on the "header HEIGHT" area of grid holding Control
    Then columns "AGE, HEIGHT" should be selected
    When user remembers the "column width of HEIGHT" reading of grid
    And user drags the "column resizer AGE" area of grid by 40 pixels to the right
    Then the "column width of AGE" reading of grid should be higher than before
    And the "column width of HEIGHT" reading of grid should not be as remembered
    And the "column width of AGE" and "column width of HEIGHT" readings of grid should be the same
    And the "column width of SEX" reading of grid should be the same as before
    When user presses Escape
    Then no columns should be selected
    When user remembers the "column width of HEIGHT" reading of grid
    And user drags the "column resizer AGE" area of grid by 30 pixels to the right
    Then the "column width of AGE" reading of grid should be higher than before
    And the "column width of HEIGHT" reading of grid should be as remembered
    When user drags the "column resizer WEIGHT" area of grid by 200 pixels to the left
    Then the "column width of WEIGHT" reading of grid should be lower than before
    And the "column width of WEIGHT" reading of grid should be between 0 and 3
    And the "column order" reading of grid should contain "WEIGHT"
    When user picks "Column Sizing > Optimal" from the context menu of the "cell 4 of AGE" area of grid
    Then the "column width of WEIGHT" reading of grid should be higher than before
    And no errors should have been logged

  Scenario: Columns widened past the view scroll sideways and draw cleanly (GROK-19753)
    Then grid should not have an "x scroll handle" area
    When user drags the "column resizer AGE" area of grid by 300 pixels to the right
    And user drags the "column resizer SEX" area of grid by 300 pixels to the right
    And user drags the "column resizer RACE" area of grid by 300 pixels to the right
    Then the "column width of RACE" reading of grid should be higher than before
    And grid should have an "x scroll handle" area
    And the "x scroll span" reading of grid should be between 0.1 and 0.95
    And grid should not have a "header SEVERITY" area
    When user drags the "x scroll handle" area of grid by 900 pixels to the right
    Then grid should have a "header SEVERITY" area
    And grid should have a "cell 4 of SEVERITY" area
    And no errors should have been logged
    When user drags the "column resizer STARTED" area of grid by 60 pixels to the right
    Then the "column width of STARTED" reading of grid should be higher than before
    And the "x scroll span" reading of grid should be lower than before
    And grid should have a "header SEVERITY" area
    And no errors should have been logged
    When user picks "Column Sizing > Optimal" from the context menu of the "cell 4 of SEVERITY" area of grid
    Then grid should not have an "x scroll handle" area
    And no errors should have been logged

  Scenario: The column manager keeps its type filter while it rebinds to another table and back (GROK-19332)
    When user opens demog-1000 dataset keeping the first 100 rows as "demog-100"
    And user removes "RACE" column
    And user removes "DIS_POP" column
    And user removes "DEMOG" column
    And user removes "SEVERITY" column
    Then the table should have 7 columns
    When user switches to the "demog-1000" table view
    And user clicks on "Columns: 11" text in status bar
    Then column manager should be visible
    And the "rows shown" reading of Grid viewer in column manager should be 11
    When user hovers over "Search" input in column manager
    And user clicks on font icon menu icon in column manager
    And user hovers over "Types" menu item in context menu
    And user clicks on "string" menu item in context menu
    Then the "rows shown" reading of Grid viewer in column manager should be 6
    When user closes the context menu
    And user switches to the "demog-100" table view
    Then the "rows shown" reading of Grid viewer in column manager should be 2
    When user switches to the "demog-1000" table view
    Then the "rows shown" reading of Grid viewer in column manager should be 6
    When user hovers over "Search" input in column manager
    And user clicks on font icon menu icon in column manager
    And user hovers over "Types" menu item in context menu
    And user clicks on "Reset filter" menu item in context menu
    Then the "rows shown" reading of Grid viewer in column manager should be 11
    When user closes the context menu
    And user switches to the "demog-100" table view
    And user closes the current view
    And user switches to the "demog-1000" table view
    Then the table should have 1000 rows
    And no errors should have been logged
