@guide @tutorial
Feature: How-to · Viewers
  For the maintainer: a pivot table set up from its three rows of chips. It opens on DIS_POP by SEVERITY,
  avg(AGE); RACE added to Group by with its + and DIS_POP's chip closed leave four rows; SEVERITY's chip
  closed and SEX added to Pivot make one column per sex (the cells checked against an independent
  group-by of demog); HEIGHT's header dragged from the grid onto Aggregate adds avg(HEIGHT), and the
  context menu of the avg(AGE) chip changes it to the maximum. ADD publishes the result as a new table,
  "demog aggregation" (4 rows), in its own view; it lives only in the browser. The + opens the column
  picker that `user picks {string} in the open column picker` (howto-viewers.ts) types into. The saved
  pivot parameters are cleared and the viewer added off camera (it docks above the grid); demog is
  opened off camera; nothing is made on the server.

  Scenario: Summarize a table with a pivot table
    Group, pivot and aggregate the rows of a table with a pivot table, then add the result as a new table.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user clears the saved pivot table parameters
    And user adds a pivot table viewer
    Then the "group by" reading of pivot table viewer should be "DIS_POP"
    And the "pivot" reading of pivot table viewer should be "SEVERITY"
    When the video says "A pivot table sums up the rows of a table by group"
    # caption: Group by: one row per value
    And user clicks on the "add group by" area of pivot table viewer
    And user picks "RACE" in the open column picker
    Then the "group by" reading of pivot table viewer should be "DIS_POP, RACE"
    # caption: The x on a chip takes it out
    When user clicks on the "remove group by chip DIS_POP" area of pivot table viewer
    Then the "group by" reading of pivot table viewer should be "RACE"
    And the "aggregated rows" reading of pivot table viewer should be 4
    # caption: Pivot: one column per value
    When user clicks on the "remove pivot chip SEVERITY" area of pivot table viewer
    And user clicks on the "add pivot" area of pivot table viewer
    And user picks "SEX" in the open column picker
    Then the "pivot" reading of pivot table viewer should be "SEX"
    And the aggregated values of pivot table viewer should match "avg(AGE)" grouped by "RACE" pivoted on "SEX"
    # caption: Aggregate: what to compute. Or drag a column from the grid
    When user drags the "header HEIGHT" area of grid onto the "aggregate row" area of pivot table viewer
    Then the "aggregate" reading of pivot table viewer should be "avg(AGE), avg(HEIGHT)"
    # caption: Right-click a chip to change its function
    When user picks "Aggregation > max" from the context menu of the "aggregate chip avg(AGE)" area of pivot table viewer
    And user closes the context menu
    Then the "aggregate" reading of pivot table viewer should be "max(AGE), avg(HEIGHT)"
    # caption: ADD puts the result into a new table
    When user clicks on the "add to workspace" area of pivot table viewer
    Then table "demog aggregation" should be open
    And the "demog aggregation" view should be current
    And the table should have 4 rows
    And no errors should have been logged
