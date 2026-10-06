@journey
Feature: The Table Manager lists the open tables and moves between them
  Alt+T docks the Table Manager, a grid of the open tables in a panel titled Tables: it lists every
  open table in the order they opened, a click on a row brings that table's view to front and makes
  the table the current object, Open as table makes a table of the list, a closed table leaves it, and
  Alt+T closes the manager again. Translated from TestTrack General/table-manager.md,
  table-manager-spec.ts and the manual-only table-manager-ui.md; the manager's grid reports the rows
  and cells of its inner grid.

  Not translated, and why: Show > All, whose second pick removes the name column along with the
  attributes (GROK-17558, reopened); the Shift multi-select, which the manager's grid (no row header)
  takes as a drag, not a click. Nothing is put on the server; the docked panel is the browser's state,
  and the journey closes it.

  Background:
    Given user is logged in
    And the context panel is open
    And user opens cars dataset
    And user opens iris dataset
    And user opens beer dataset
    Then the open table views should be exactly "cars, iris, beer"

  Scenario: Alt+T docks the manager with every open table, in the order they opened
    Then "Tables" dock panel should be absent
    When user presses Alt+T
    Then the "rows" reading of Grid viewer in "Tables" dock panel should be 3
    And the "text of cell 1 of name" reading of Grid viewer in "Tables" dock panel should be "cars"
    And the "text of cell 2 of name" reading of Grid viewer in "Tables" dock panel should be "iris"
    And the "text of cell 3 of name" reading of Grid viewer in "Tables" dock panel should be "beer"
    And no errors should have been logged

  Scenario: A click on a row brings its table's view to front and makes the table current
    When user clicks on the "cell 1 of name" area of Grid viewer in "Tables" dock panel
    Then the "cars" view should be current
    And the context panel should show "cars"
    When user clicks on the "cell 2 of name" area of Grid viewer in "Tables" dock panel
    Then the "iris" view should be current
    And the context panel should show "iris"
    And no errors should have been logged

  Scenario: Open as table makes a table of the manager's list, and a closed table leaves the list
    When user picks "Open as table" from the context menu of the "cell 1 of name" area of Grid viewer in "Tables" dock panel
    Then the "rows" reading of Grid viewer in "Tables" dock panel should be 4
    And the table should have 3 rows
    And the value of "name" column in row 1 should be "cars"
    And the value of "name" column in row 3 should be "beer"
    When user closes the current view
    Then the open table views should be exactly "cars, iris, beer"
    And the "rows" reading of Grid viewer in "Tables" dock panel should be 3
    And no errors should have been logged

  Scenario: Alt+T closes the manager again
    When user presses Alt+T
    Then "Tables" dock panel should be absent
    And no errors should have been logged
