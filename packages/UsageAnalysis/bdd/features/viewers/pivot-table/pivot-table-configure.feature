@journey @viewers @realizes:viewers.pivot-table
Feature: Pivot table — configuring the tag rows
  Every way the cross tab is configured: the + picker of a row, a column header dragged from the main
  grid onto a row, the chip's own context menu (which stays open after a pick — GROK-16899 — and
  rebuilds its Aggregation group when the chip is repointed at another column), Remove others, the
  Pivot row that disappears with the last aggregate and takes the pivot columns with it, and the
  Refresh icon that reseeds the type-driven defaults. One journey on demog-1000, which starts as
  DIS_POP / avg(AGE) / SEVERITY: adding SEX to Group by makes 12 groups, RACE makes 24, and Refresh
  leaves no grouping at all with the first two numerical columns, avg(AGE) and avg(HEIGHT).
  Not translated: the property panel's three column-list editors (`prop-group-by`, `prop-aggregate`,
  `prop-pivot`) as DOM widgets — the last scenario writes the same lists as properties and claims the
  chips the viewer rebuilt from them, which is what the editors are there to do; and the old spec's
  "the panel exposes the three editors", a presence check the same scenario subsumes.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pivot table viewer
    Then the "group by" reading of pivot table viewer should be "DIS_POP"
    And the "aggregate" reading of pivot table viewer should be "avg(AGE)"
    And the "pivot" reading of pivot table viewer should be "SEVERITY"
    And the "aggregated rows" reading of pivot table viewer should be 6

  Scenario: The + picker of the Group by row adds a second key column
    When user adds "SEX" to the "group by" row of pivot table viewer
    Then the "group by" reading of pivot table viewer should be "DIS_POP, SEX"
    And the "key columns" reading of pivot table viewer should be "DIS_POP, SEX"
    And pivot table viewer should have a "group by chip SEX" area
    And the "aggregated rows" reading of pivot table viewer should be 12
    And the "text of grid cell 1 of SEX" reading of pivot table viewer should be "F"
    When user clicks on the "remove group by chip SEX" area of pivot table viewer
    Then the "group by" reading of pivot table viewer should be "DIS_POP"
    And the "aggregated rows" reading of pivot table viewer should be 6
    And no errors should have been logged

  Scenario: A column header dragged from the main grid onto the Group by row becomes a key
    When user drags the "header RACE" area of grid onto the "group by row" area of pivot table viewer
    Then the "group by" reading of pivot table viewer should be "DIS_POP, RACE"
    And pivot table viewer should have a "group by chip RACE" area
    And the "aggregated rows" reading of pivot table viewer should be 24
    When user clicks on the "remove group by chip RACE" area of pivot table viewer
    Then the "group by" reading of pivot table viewer should be "DIS_POP"
    And the "aggregated rows" reading of pivot table viewer should be 6
    And no errors should have been logged

  Scenario: A column header dragged onto the Aggregate row becomes an aggregation
    When user drags the "header HEIGHT" area of grid onto the "aggregate row" area of pivot table viewer
    Then the "aggregate" reading of pivot table viewer should be "avg(AGE), avg(HEIGHT)"
    And the "aggregated columns" reading of pivot table viewer should be 11
    When user clicks on the "remove aggregate chip avg(HEIGHT)" area of pivot table viewer
    Then the "aggregate" reading of pivot table viewer should be "avg(AGE)"
    And the "aggregated columns" reading of pivot table viewer should be 6
    And no errors should have been logged

  Scenario: The chip menu takes two aggregation picks without closing
    When user picks "Aggregation > sum" from the context menu of the "aggregate chip avg(AGE)" area of pivot table viewer
    Then the "aggregate" reading of pivot table viewer should be "sum(AGE)"
    And the "aggregations" reading of pivot table viewer should be "sum(AGE)"
    And the aggregated values of pivot table viewer should match "sum(AGE)" grouped by "DIS_POP" pivoted on "SEVERITY"
    And the open menu should list "Aggregation > med"
    When user clicks on "med" menu item in context menu
    Then the "aggregate" reading of pivot table viewer should be "med(AGE)"
    And the aggregated values of pivot table viewer should match "med(AGE)" grouped by "DIS_POP" pivoted on "SEVERITY"
    When user closes the context menu
    And user picks "Aggregation > avg" from the context menu of the "aggregate chip med(AGE)" area of pivot table viewer
    Then the "aggregate" reading of pivot table viewer should be "avg(AGE)"
    When user closes the context menu
    Then no errors should have been logged

  Scenario: Repointing the chip at another column rebuilds the offered aggregations
    When user picks "Column > HEIGHT" from the context menu of the "aggregate chip avg(AGE)" area of pivot table viewer
    Then the "aggregate" reading of pivot table viewer should be "avg(HEIGHT)"
    And the open menu should list "Aggregation > geomean"
    And the open menu should list "Aggregation > stdev"
    And "Aggregation > avg" menu item in context menu should be selected
    And "Aggregation > sum" menu item in context menu should not be selected
    When user closes the context menu
    Then the aggregated values of pivot table viewer should match "avg(HEIGHT)" grouped by "DIS_POP" pivoted on "SEVERITY"
    When user picks "Column > AGE" from the context menu of the "aggregate chip avg(HEIGHT)" area of pivot table viewer
    And user closes the context menu
    Then the "aggregate" reading of pivot table viewer should be "avg(AGE)"
    And no errors should have been logged

  Scenario: Remove others leaves the chip it was opened on
    When user adds "HEIGHT" to the "aggregate" row of pivot table viewer
    Then the "aggregate" reading of pivot table viewer should be "avg(AGE), avg(HEIGHT)"
    When user picks "Remove others" from the context menu of the "aggregate chip avg(AGE)" area of pivot table viewer
    Then the "aggregate" reading of pivot table viewer should be "avg(AGE)"
    And the "aggregated columns" reading of pivot table viewer should be 6
    And pivot table viewer should not have a "aggregate chip avg(HEIGHT)" area
    And no errors should have been logged

  Scenario: Removing the last aggregate hides the Pivot row and clears the pivot columns
    When user clicks on the "remove aggregate chip avg(AGE)" area of pivot table viewer
    Then the "aggregate" reading of pivot table viewer should be ""
    And the "pivot" reading of pivot table viewer should be ""
    And the "pivot row shown" reading of pivot table viewer should be "false"
    And pivot table viewer should not have a "pivot row" area
    And the "aggregated columns" reading of pivot table viewer should be 1
    And the "text of grid cell 5 of DIS_POP" reading of pivot table viewer should be "RA"
    And no errors should have been logged

  Scenario: Re-adding an aggregate brings the Pivot row back, empty
    When user adds "AGE" to the "aggregate" row of pivot table viewer
    Then the "aggregate" reading of pivot table viewer should be "avg(AGE)"
    And the "pivot row shown" reading of pivot table viewer should be "true"
    And pivot table viewer should have a "pivot row" area
    And the "pivot" reading of pivot table viewer should be ""
    And the "aggregated columns" reading of pivot table viewer should be 2
    When user adds "SEVERITY" to the "pivot" row of pivot table viewer
    Then the "pivot" reading of pivot table viewer should be "SEVERITY"
    And the "aggregated columns" reading of pivot table viewer should be 6
    And no errors should have been logged

  Scenario: Refresh reseeds the first two numerical columns and clears the rest
    When user clicks on the "refresh" area of pivot table viewer
    Then the "group by" reading of pivot table viewer should be ""
    And the "pivot" reading of pivot table viewer should be ""
    And the "aggregate" reading of pivot table viewer should be "avg(AGE), avg(HEIGHT)"
    And the "aggregated rows" reading of pivot table viewer should be 1
    And the "aggregated columns" reading of pivot table viewer should be 2
    And no errors should have been logged

  Scenario: The column lists written as properties rebuild the chips
    When user sets "Group By Column Names" property of pivot table viewer to "DIS_POP"
    And user sets "Aggregate Column Names" property of pivot table viewer to "WEIGHT"
    And user sets "Aggregate Agg Types" property of pivot table viewer to "avg"
    And user sets "Pivot Column Names" property of pivot table viewer to "SEVERITY"
    Then the "group by" reading of pivot table viewer should be "DIS_POP"
    And the "aggregate" reading of pivot table viewer should be "avg(WEIGHT)"
    And the "pivot" reading of pivot table viewer should be "SEVERITY"
    And pivot table viewer should have a "aggregate chip avg(WEIGHT)" area
    And the "aggregated rows" reading of pivot table viewer should be 6
    And the aggregated values of pivot table viewer should match "avg(WEIGHT)" grouped by "DIS_POP" pivoted on "SEVERITY"
    When user sets "Aggregate Column Names" property of pivot table viewer to "AGE"
    And user sets "Aggregate Agg Types" property of pivot table viewer to "avg"
    Then the "aggregate" reading of pivot table viewer should be "avg(AGE)"
    And the "text of grid cell 5 of None avg(AGE)" reading of pivot table viewer should be "52.30"
    And no errors should have been logged
