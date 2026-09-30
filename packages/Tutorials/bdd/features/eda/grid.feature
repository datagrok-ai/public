@tutorials @serial @realizes:tutorials.grid
Feature: The Grid Customization tutorial
  Walks Exploratory Data Analysis > Grid Customization from its card to the end: keyboard
  navigation, a row added and edited, row selection by clicks and keys, the rows with one value
  selected and deleted, column selection, sorting, reordering, resizing, formats, colour coding and
  a summary column. Each step is claimed as ticked and as done on the grid and the table — the
  current cell, the values typed, the exact rows selected, 2550 rows left without a "None", the
  column order, the sizes against before, the format tags, the new summary column.
  Translated from playwright-tests/e2e/tutorials/grid.test.ts, which located the grid's cells by
  hit-testing its canvas pixel by pixel.

  Reordering needs no signal of its own: a header dropped anywhere over the row-number column
  inserts the dragged columns first, and the selected columns travel together.
  Two hints the tutorial captured when their step began (the new-row icon, the remove-rows icon)
  follow their controls now.

  Serial, because a finished tutorial writes its completion record into the account's settings,
  which every page syncs whole.

  Background:
    Given user is logged in
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "Grid Customization" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner completes the Grid Customization tutorial
    When user starts the "Grid Customization" tutorial
    Then the tutorial progress should be 1 of 31

    When user clicks on the "cell 1 of AGE" area of grid
    And user presses ArrowDown
    And user presses ArrowDown
    And user presses ArrowDown
    Then the tutorial step "Go a few rows down in the grid" should be done
    And row 4 should be current
    When user presses Control+Home
    Then the tutorial step "Jump back to the first row" should be done
    And row 1 should be current
    When user presses Control+End
    Then the tutorial step "Jump to the last row" should be done
    And row 5850 should be current
    When user presses End
    Then the tutorial step "Jump to the last column" should be done
    And the current column should be "SEVERITY"
    When user presses Home
    Then the tutorial step "Jump to the first column" should be done
    And the current column should be "USUBJID"

    Then "Add new row" icon should be hinted
    When user clicks on "Add new row" icon
    Then the tutorial step "Add a new row by clicking the \"+\" icon at the last row" should be done
    And the table should have 5851 rows
    When user double-clicks on the "cell 5851 of USUBJID" area of grid
    And user types "X0273T51080200024" into cell editor
    And user presses Enter
    Then the tutorial step "Enter \"X0273T51080200024\" as USUBJID in the new row (#5851)" should be done
    And the value of "USUBJID" column in row 5851 should be "X0273T51080200024"
    # editing the last row adds one more (Add New Row On Last Row Edit)
    And the table should have 5852 rows
    When user double-clicks on the "cell 5851 of AGE" area of grid
    And user types "37" into cell editor
    And user presses Enter
    Then the tutorial step "Set AGE to \"37\" in the row you edited (#5851)" should be done
    When user double-clicks on the "cell 5851 of SEX" area of grid
    And user types "F" into cell editor
    And user presses Enter
    Then the tutorial step "Set SEX to \"F\" in the row you edited (#5851)" should be done
    And the value of "SEX" column in row 5851 should be "F"
    When user clicks on the "cell 5850 of RACE" area of grid
    And user presses Control+C
    And user clicks on the "cell 5851 of RACE" area of grid
    And user presses Control+V
    Then the tutorial step "Copy the value from row #5850 to row #5851 for the RACE column" should be done
    And the value of "RACE" column in row 5851 should be "Caucasian"

    # a plain click on a row number makes the row current; Shift selects from the current row (5851)
    When user clicks on the "row header 5851" area of grid holding Shift
    Then the tutorial step "Select the edited row: hold Shift and click its number (#5851)" should be done
    And rows 5851 to 5851 should be selected
    When user presses Escape
    Then the tutorial step "Remove selection by pressing \"Esc\"" should be done
    And no rows should be selected
    When user presses Control+A
    Then the tutorial step "Select all rows with \"Ctrl+A\"" should be done
    And all rows should be selected
    When user presses Escape
    And user presses Control+Home
    And user clicks on the "row header 1" area of grid
    And user clicks on the "row header 5" area of grid holding Shift
    And user presses Control+End
    And user clicks on the "row header 5848" area of grid holding Control
    And user clicks on the "row header 5849" area of grid holding Control
    And user clicks on the "row header 5850" area of grid holding Control
    And user clicks on the "row header 5851" area of grid holding Control
    And user clicks on the "row header 5852" area of grid holding Control
    Then the tutorial step "Select the first five and the last five table rows (#1 - #5 and #5848 - #5852)" should be done
    And 10 rows should be selected
    And rows 1 to 5 should all be selected
    And rows 5848 to 5852 should all be selected
    When user presses Escape
    Then the tutorial step "Clear the selection" should be done

    # the tutorial opens the context panel here; beside it, the Tutorials panel and the toolbox, the
    # last columns lie past the grid's edge: the learner makes room
    Given the toolbox pane is hidden
    When user presses F4
    Then context panel should be hidden
    # row 1 of demog is a "None" severity; 3302 rows share it
    Given the tutorial step "Find the SEVERITY column and select all rows with the \"None\" value" should not be done yet
    When user presses Control+Home in grid
    And user presses End in grid
    And user clicks on the "cell 1 of SEVERITY" area of grid
    And user presses Shift+Enter in grid
    Then the tutorial step "Find the SEVERITY column and select all rows with the \"None\" value" should be done
    And 3302 rows should be selected
    And all rows where "SEVERITY" is "None" should be selected
    When user presses Shift+Delete
    Then the tutorial step "Delete these rows (3302)" should be done
    And the table should have 2550 rows
    And the table should have no rows where "SEVERITY" is "None"

    When user presses Home
    And user clicks on the "header HEIGHT" area of grid holding Shift
    And user clicks on the "header WEIGHT" area of grid holding Shift
    Then the tutorial step "Select HEIGHT and WEIGHT columns" should be done
    And columns "HEIGHT, WEIGHT" should be selected
    When user presses Escape
    Then the tutorial step "Clear the selection" should be done 2 times
    And no columns should be selected

    When user clicks on the "cell 1 of SEX" area of grid
    And user presses Home
    And user double-clicks on the "header AGE" area of grid
    Then the "sort column" reading of grid should be "AGE"
    And the "sort direction" reading of grid should be "descending"
    And the tutorial step "Sort subjects by AGE in the descending order" should be done
    When user double-clicks on the "header AGE" area of grid
    Then the "sort direction" reading of grid should be "ascending"
    When user double-clicks on the "header AGE" area of grid
    Then the tutorial step "Reset sorting in the grid" should be done
    And the "sort column" reading of grid should be ""

    When user clicks on the "header HEIGHT" area of grid holding Shift
    And user clicks on the "header WEIGHT" area of grid holding Shift
    And user clicks on the "header STARTED" area of grid holding Shift
    And user drags the "header HEIGHT" area of grid to the "row header 1" area
    Then the tutorial step "Move columns HEIGHT, WEIGHT and STARTED to the beginning of the grid" should be done
    And the "column order" reading of grid should include the text "HEIGHT, WEIGHT, STARTED, USUBJID"
    When user presses Escape
    And user presses F4
    Then context panel should be visible

    When user drags the "row resizer 3" area of grid by 20 pixels to the down
    Then the tutorial step "Increase row height" should be done
    And the "row header 3" area of grid should be taller than before
    When user drags the "column resizer STARTED" area of grid by 40 pixels to the right
    Then the tutorial step "Extend the STARTED column width" should be done
    And the "column width of STARTED" reading of grid should be higher than before

    # the format items show the value of the first row in that format
    When user picks "Format > scientific (1.58E2)" from the context menu of the "header HEIGHT" area of grid
    Then the tutorial step "Change number formatting of HEIGHT to scientific notation" should be done
    And "HEIGHT" column should have tag "format" equal to "scientific"
    When user picks "Format > max two digits after comma (62.7)" from the context menu of the "header WEIGHT" area of grid
    Then the tutorial step "Set number formatting of WEIGHT to \"max two digits after comma\"" should be done
    And "WEIGHT" column should have tag "format" equal to "max two digits after comma"
    When user picks "Format > dd.MM.yyyy (15.05.1990)" from the context menu of the "header STARTED" area of grid
    Then the tutorial step "Set date formatting of STARTED to \"dd.MM.yyyy\"" should be done
    And "STARTED" column should have tag "format" equal to "dd.MM.yyyy"
    When user picks "Color Coding > Linear" from the context menu of the "header HEIGHT" area of grid
    Then the tutorial step "Set linear color coding for HEIGHT" should be done
    And "HEIGHT" column should be color-coded linearly

    When user picks "Add > Summary Columns > Bar Chart" from the context menu of the "cell 3 of AGE" area of grid
    Then the tutorial step "Add a summary column with inline bar chart" should be done
    And the "column order" reading of grid should include the text "Bar Chart"
    When user clicks on the "header Bar Chart" area of grid
    Then the context panel should show "Bar Chart"
    And the tutorial step "Click the \"Bar Chart\" column header" should be done

    And the "Grid Customization" tutorial should be completed
    And the tutorial should have listed 30 steps
    And the tutorial progress should be 31 of 31
    And no hint should be shown
    And no errors should have been logged
