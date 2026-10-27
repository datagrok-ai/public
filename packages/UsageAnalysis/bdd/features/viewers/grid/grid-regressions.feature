@viewers @realizes:viewers.grid
Feature: Grid regressions from closed issues
  Closed grid issues, each claimed where the bug showed, on demog-1000. With Sync New
  Columns off a new calculated column is in the table but not in the grid, and with it on the grid
  shows it (public#2348). A column narrowed below its name still shows its header tooltip, and the
  header menu takes the tooltip away (public#1130). A pinned column sorts by a double-click on its
  header, descending, ascending and off, and takes a resize (public#1159). A project with a second
  grid on a second table, sorted, opens again with the sort and no error (public#2718).
  Control+Shift+Home and Control+Shift+End select from the current row to the first and to the last
  row (public#2663). Enter and a typed digit open the editor on the current cell, and F2 opens the
  column's Column Properties dialog, as the help says (public#629). A cut-off cell's value in its
  tooltip (public#989, public#996) shows an empty tooltip on the current master: see MISSING.md.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then grid should show 1000 rows

  # public#2348
  @realizes:public#2348
  Scenario: Sync New Columns decides whether a new calculated column is drawn (public#2348)
    When user sets "Sync New Columns" property of grid to "false"
    And user adds a calculated column "H2" with formula "${HEIGHT} * 2"
    Then the table should have a column "H2"
    And the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    When user sets "Sync New Columns" property of grid to "true"
    And user adds a calculated column "H3" with formula "${HEIGHT} * 3"
    Then the table should have a column "H3"
    And the "column order" reading of grid should include the text "H3"
    And the "column order" reading of grid should not include the text "H2"
    And no errors should have been logged

  # public#1130
  @realizes:public#1130
  Scenario: A column narrowed below its name still shows its header tooltip, and the header menu hides it (public#1130)
    When user drags the "column resizer RACE" area of grid by 60 pixels to the left
    Then the "column width of RACE" reading of grid should be lower than before
    When user hovers over the "header RACE" area of grid
    Then tooltip should be visible
    And tooltip should contain text "RACE"
    When user right-clicks on the "header RACE" area of grid
    Then tooltip should be hidden
    When user closes the context menu
    Then no errors should have been logged

  # public#1159
  @realizes:public#1159
  Scenario: A pinned column sorts by its header, descending, ascending and off, and takes a resize (public#1159)
    When user picks "Pin > Pin Column" from the context menu of the "header AGE" area of grid
    Then "Frozen Columns" property of grid should be "2"
    When user double-clicks on the "header AGE" area of grid
    Then the "sort column" reading of grid should be "AGE"
    And the "sort direction" reading of grid should be "descending"
    When user double-clicks on the "header AGE" area of grid
    Then the "sort direction" reading of grid should be "ascending"
    When user double-clicks on the "header AGE" area of grid
    Then the "sort column" reading of grid should be ""
    When user drags the "column resizer AGE" area of grid by 40 pixels to the right
    Then the "column width of AGE" reading of grid should be higher than before
    And no errors should have been logged

  # public#2718
  @realizes:public#2718
  Scenario: A view with a sorted grid of a second table opens again from its project (public#2718)
    Given the user's own project "bdd grid two tables {time}" is removed now and at feature end
    And user opens demog-1000 dataset keeping the first 100 rows as "demog-100"
    And user switches to the "demog-1000" table view
    And user adds a grid viewer
    When user sets "table" property of second grid viewer to "demog-100"
    Then second grid viewer should be bound to table "demog-100"
    When user double-clicks on the "header AGE" area of second grid viewer
    Then the "sort column" reading of second grid viewer should be "AGE"
    When user saves all open table views as project "bdd grid two tables {time}"
    And user closes all views
    And user opens the "bdd grid two tables {time}" project
    Then second grid viewer should be bound to table "demog-100"
    And the "sort column" reading of second grid viewer should be "AGE"
    And no viewer of the current view should report an error
    And no errors should have been logged

  # public#2663
  @realizes:public#2663
  Scenario: Control+Shift+Home and Control+Shift+End select from the current row to the edge (public#2663)
    Given user selects no rows
    When user clicks on the "cell 5 of AGE" area of grid
    Then row 5 should be current
    When user presses Control+Shift+Home
    Then 5 rows should be selected
    And rows 1 to 5 should all be selected
    When user clears the row selection
    And user clicks on the "cell 5 of AGE" area of grid
    And user presses Control+Shift+End
    Then 996 rows should be selected
    And rows 5 to 1000 should all be selected
    And no errors should have been logged

  # public#629
  @realizes:public#629
  Scenario: Enter and a typed digit open the editor on the current cell (public#629)
    When user clicks on the "cell 4 of AGE" area of grid
    Then cell editor should be hidden
    When user presses Enter
    Then cell editor should be visible
    When user presses Escape
    Then cell editor should be hidden
    And the value of "AGE" column in row 4 should be "45"
    When user presses 7
    Then cell editor should be visible
    When user presses Escape
    Then cell editor should be hidden
    And the value of "AGE" column in row 4 should be "45"
    When user presses F2
    Then "AGE" dialog should be visible
    And cell editor should be hidden
    When user presses Escape
    Then "AGE" dialog should be absent
    And no errors should have been logged
