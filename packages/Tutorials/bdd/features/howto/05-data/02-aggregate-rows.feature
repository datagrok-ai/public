@guide @tutorial
Feature: How-to · Data
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). Data > Aggregate Rows adds a pivot table to the view,
  set up on its own from the column types (DIS_POP by SEVERITY, avg(AGE)). The clip regroups it by RACE, pivots
  it on SEX, adds WEIGHT and switches that aggregate to the median from the chip's menu, then ADD publishes the
  result as a table view of its own ("demog aggregation": 4 races). demog is opened off camera and the saved
  pivot parameters are cleared first; nothing is made on the server.

  Scenario: Aggregate rows
    Group rows by one column, spread them over another, and choose what to compute for each group.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And user clears the saved pivot table parameters
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: Aggregate Rows: summarize the table group by group
    When user picks "Data > Aggregate Rows..." from the top menu
    Then the open tableview should have 1 pivot table viewer
    And the "group by" reading of pivot table viewer should be "DIS_POP"
    # caption: Group by: one result row per value
    When user clicks on the "remove group by chip DIS_POP" area of pivot table viewer
    And user picks "RACE" with the plus of the "group by" row of pivot table viewer
    Then the "group by" reading of pivot table viewer should be "RACE"
    # caption: Pivot: one result column per value
    When user clicks on the "remove pivot chip SEVERITY" area of pivot table viewer
    And user picks "SEX" with the plus of the "pivot" row of pivot table viewer
    Then the "pivot" reading of pivot table viewer should be "SEX"
    # caption: Aggregate: what to compute for each group
    When user picks "WEIGHT" with the plus of the "aggregate" row of pivot table viewer
    Then the "aggregate" reading of pivot table viewer should be "avg(AGE), avg(WEIGHT)"
    # caption: Right-click an aggregate to change its function
    When user picks "Aggregation > med" from the context menu of the "aggregate chip avg(WEIGHT)" area of pivot table viewer
    And user closes the context menu
    Then the "aggregate" reading of pivot table viewer should be "avg(AGE), med(WEIGHT)"
    And the "aggregated rows" reading of pivot table viewer should be 4
    When the video says "ADD turns the result into a table of its own" at the bottom
    And user hovers over the "add to workspace" area of pivot table viewer
    And user clicks on the "add to workspace" area of pivot table viewer
    Then table "demog aggregation" should be open
    And the "demog aggregation" view should be current
    And the table should have 4 rows
    And no errors should have been logged
