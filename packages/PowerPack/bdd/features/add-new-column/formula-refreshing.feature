@journey @realizes:powerpack.cp.add-new-column-persists @realizes:GROK-17109
Feature: A chain of calculated columns recalculates on a formula edit and survives a project round trip
  demog opened from Browse > Files > Demo, so the table has a file behind it. Three columns are added
  through the Add New Column dialog, each on the one before: Weight2 = ${WEIGHT} + 100, Weight3 =
  ${Weight2} + 100, Weight4 = Log10(${Weight3}) - 0.2, and each is checked row by row against its
  source. Each formula is then edited in the Formula pane of the column's context panel, and Apply
  recalculates that column and everything downstream of it, leaving what is upstream as it was. The
  view is saved as a project from the toolbar with Data sync on, closed, and reopened: the three
  columns come back with the last formulas, their Formula panes, and values recomputed from the file
  (GROK-17109). Translated from TestTrack PowerPack/formula-refreshing.md.

  Before the save one WEIGHT cell is edited in the grid; the reopened project reads the file again,
  so that cell holds the file's value and the whole chain is recomputed from it. The views are closed
  through the platform's API (the workspace then holds no table), and the project is reopened from
  Browse > Dashboards.

  Background:
    Given user is logged in
    And no project named "bdd-anc-chain-{run}" is on the server
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks on Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    And the table should have 5850 rows

  Scenario: Three columns, each computed from the one before
    When user clicks on "Add New Column..." icon
    And user types "Weight2" into column name input
    And user types "${WEIGHT} + 100" into formula editor
    And user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And every value of "Weight2" column should equal "WEIGHT" column plus 100
    When user clicks on "Add New Column..." icon
    And user types "Weight3" into column name input
    And user types "${Weight2} + 100" into formula editor
    And user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And every value of "Weight3" column should equal "WEIGHT" column plus 200
    When user clicks on "Add New Column..." icon
    And user types "Weight4" into column name input
    And user types "Log10(${Weight3}) - 0.2" into formula editor
    And user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And every value of "Weight4" column should be the decimal log of "Weight3" column minus 0.2
    And no errors should have been logged

  Scenario: An edit of Weight2 in its Formula pane recalculates Weight3 and Weight4
    Given the context panel is open
    When user clicks on the "header Weight2" area of grid
    Then the context panel should show "Weight2"
    Given Formula pane in context panel is expanded
    Then formula pane editor in Formula pane in context panel should hold the formula "${WEIGHT} + 100"
    When user types "${WEIGHT} + 200" into formula pane editor in Formula pane in context panel
    And user clicks on Apply button in Formula pane in context panel
    Then "Weight2" column should have tag "formula" equal to "${WEIGHT} + 200"
    And every value of "Weight2" column should equal "WEIGHT" column plus 200
    And every value of "Weight3" column should equal "WEIGHT" column plus 300
    And every value of "Weight4" column should be the decimal log of "Weight3" column minus 0.2
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: An edit of Weight3 recalculates Weight4 and leaves Weight2
    When user clicks on the "header Weight3" area of grid
    Then the context panel should show "Weight3"
    Given Formula pane in context panel is expanded
    When user types "${Weight2} + 50" into formula pane editor in Formula pane in context panel
    And user clicks on Apply button in Formula pane in context panel
    Then "Weight3" column should have tag "formula" equal to "${Weight2} + 50"
    And every value of "Weight3" column should equal "Weight2" column plus 50
    And every value of "Weight2" column should equal "WEIGHT" column plus 200
    And every value of "Weight4" column should be the decimal log of "Weight3" column minus 0.2
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: An edit of Weight4 changes Weight4 only
    When user clicks on the "header Weight4" area of grid
    Then the context panel should show "Weight4"
    Given Formula pane in context panel is expanded
    When user types "Log10(${Weight3}) - 0.1" into formula pane editor in Formula pane in context panel
    And user clicks on Apply button in Formula pane in context panel
    Then "Weight4" column should have tag "formula" equal to "Log10(${Weight3}) - 0.1"
    And every value of "Weight4" column should be the decimal log of "Weight3" column minus 0.1
    And every value of "Weight3" column should equal "WEIGHT" column plus 250
    And every value of "Weight2" column should equal "WEIGHT" column plus 200
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Saved with Data sync, closed and reopened, the chain comes back as last edited
    When user double-clicks on the "cell 1 of WEIGHT" area of grid
    And user presses Control+A in cell editor
    And user types "500" into cell editor
    And user presses Enter
    Then the value of "WEIGHT" column in row 1 should be "500"
    And the value of "Weight2" column in row 1 should be "700"
    And the value of "Weight3" column in row 1 should be "750"
    When user clicks on Save button in toolbar
    Then "Save project" dialog should be visible
    When user enters "bdd-anc-chain-{run}" into Name text input in "Save project" dialog
    And user switches on Data sync input in "Save project" dialog
    Then Data sync input in "Save project" dialog should be switched on
    When user clicks on OK button in "Save project" dialog
    Then "Save project" dialog should be hidden
    And no error or warning balloon should have been shown
    And 1 project named "bdd-anc-chain-{run}" should be on the server
    When user presses Escape
    And user closes all views
    Then no table should be open
    When user clicks on Dashboards tree node inside browse tree
    Then the "Projects" view should be current
    When user types "bdd-anc-chain-{run}" into gallery search
    Then "bdd-anc-chain-{run}" project card should become visible within 60 seconds
    When user double-clicks on "bdd-anc-chain-{run}" project card
    Then the "demog" view should be current
    And the table should have 5850 rows
    # Data sync reads the file again: the edited cell holds what the file holds, and the chain follows it
    And the value of "WEIGHT" column in row 1 should be "73.19999694824219"
    And the value of "Weight2" column in row 1 should be "273.20001220703125"
    And "Weight2" column should have tag "formula" equal to "${WEIGHT} + 200"
    And "Weight3" column should have tag "formula" equal to "${Weight2} + 50"
    And "Weight4" column should have tag "formula" equal to "Log10(${Weight3}) - 0.1"
    And every value of "Weight2" column should equal "WEIGHT" column plus 200
    And every value of "Weight3" column should equal "Weight2" column plus 50
    And every value of "Weight4" column should be the decimal log of "Weight3" column minus 0.1
    Given the context panel is open
    When user clicks on the "header Weight3" area of grid
    Then the context panel should show "Weight3"
    Given Formula pane in context panel is expanded
    Then formula pane editor in Formula pane in context panel should hold the formula "${Weight2} + 50"
    And no error or warning balloon should have been shown
    And no errors should have been logged
