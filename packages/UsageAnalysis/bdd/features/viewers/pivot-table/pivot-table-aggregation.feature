@journey @viewers @realizes:viewers.pivot-table
Feature: Pivot table — the aggregation it publishes and the parameters it remembers
  What the pivot computes, checked against an independent groupBy of the same table rather than
  against its own look: the cross tab on the grid, the table ADD publishes into the workspace (whose
  key column keeps the type of the column it came from — GROK-16074), an identifier column that makes
  one group per row (GROK-16201), and the configurations the command bar's history saves, offers
  again and drops once a column they name is gone. One journey on demog-1000: DIS_POP has 6
  categories and USUBJID 1000, ADD publishes "demog-1000 aggregation", and avg(WEIGHT) by RACE reads
  70.52 for Asian, the first of the four alphabetical groups.
  Not translated: the old spec's raw `localStorage['grok-aggregation-history']` read — the history is
  claimed here through the entries the menu actually offers, which is the same store filtered the way
  the product filters it.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user clears the saved pivot table parameters
    And user adds a pivot table viewer
    Then the "group by" reading of pivot table viewer should be "DIS_POP"
    And the "aggregate" reading of pivot table viewer should be "avg(AGE)"
    And the "aggregated rows" reading of pivot table viewer should be 6

  Scenario: ADD publishes the cross tab it drew
    Then the aggregated values of pivot table viewer should match "avg(AGE)" grouped by "DIS_POP" pivoted on "SEVERITY"
    When user clicks on the "add to workspace" area of pivot table viewer
    Then table "demog-1000 aggregation" should be open
    And table "demog-1000 aggregation" should have 6 rows
    And table "demog-1000 aggregation" should have no missing values in "DIS_POP" column
    And the "demog-1000 aggregation" view should be current
    And the value of "DIS_POP" column in row 1 should be "AS"
    And the value of "DIS_POP" column in row 5 should be "RA"
    And table "demog-1000 aggregation" should have no missing values in "None avg(AGE)" column
    And every value of "None avg(AGE)" column should lie between 36.91 and 52.31
    When user closes the current view
    And user switches to the "demog-1000" table view
    Then the "aggregated rows" reading of pivot table viewer should be 6
    And no errors should have been logged

  Scenario: The published key column keeps the type of the column it groups
    When user sets "Pivot Column Names" property of pivot table viewer to ""
    Then the "aggregated columns" reading of pivot table viewer should be 2
    When user clicks on the "add to workspace" area of pivot table viewer
    Then the "demog-1000 aggregation" view should be current
    And "DIS_POP" column should have type "string"
    And "avg(AGE)" column should have type "double"
    And the table should have 6 rows
    When user closes the current view
    And user switches to the "demog-1000" table view
    Then the aggregated values of pivot table viewer should match "avg(AGE)" grouped by "DIS_POP"
    And no errors should have been logged

  Scenario: Grouping by an identifier column makes one row per identifier
    When user sets "Group By Column Names" property of pivot table viewer to "USUBJID"
    Then the "aggregated rows" reading of pivot table viewer should be 1000
    And the "limited" reading of pivot table viewer should be "false"
    And the "rows shown" reading of pivot table viewer should be 1000
    And the "aggregated columns" reading of pivot table viewer should be 2
    When user sets "Group By Column Names" property of pivot table viewer to "DIS_POP"
    And user sets "Pivot Column Names" property of pivot table viewer to "SEVERITY"
    Then the "aggregated rows" reading of pivot table viewer should be 6
    And no errors should have been logged

  Scenario: Save parameters, and the history menu puts the configuration back
    When user sets "Group By Column Names" property of pivot table viewer to "RACE"
    And user sets "Aggregate Column Names" property of pivot table viewer to "WEIGHT"
    And user sets "Aggregate Agg Types" property of pivot table viewer to "avg"
    And user sets "Pivot Column Names" property of pivot table viewer to ""
    Then the "aggregate" reading of pivot table viewer should be "avg(WEIGHT)"
    And the "aggregated rows" reading of pivot table viewer should be 4
    When user picks "Save parameters" from the history menu of pivot table viewer
    Then the "history entries" reading of pivot table viewer should be "key(RACE),avg(WEIGHT)"
    When user sets "Group By Column Names" property of pivot table viewer to "DIS_POP"
    And user sets "Aggregate Column Names" property of pivot table viewer to "AGE"
    And user sets "Aggregate Agg Types" property of pivot table viewer to "avg"
    And user sets "Pivot Column Names" property of pivot table viewer to "SEVERITY"
    Then the "aggregate" reading of pivot table viewer should be "avg(AGE)"
    When user picks "key(RACE),avg(WEIGHT)" from the history menu of pivot table viewer
    Then the "group by" reading of pivot table viewer should be "RACE"
    And the "aggregate" reading of pivot table viewer should be "avg(WEIGHT)"
    And the "aggregated rows" reading of pivot table viewer should be 4
    And the "text of grid cell 1 of RACE" reading of pivot table viewer should be "Asian"
    And the "text of grid cell 1 of avg(WEIGHT)" reading of pivot table viewer should be "70.52"
    And the aggregated values of pivot table viewer should match "avg(WEIGHT)" grouped by "RACE"
    When user clears the saved pivot table parameters
    And user sets "Group By Column Names" property of pivot table viewer to "DIS_POP"
    And user sets "Aggregate Column Names" property of pivot table viewer to "AGE"
    And user sets "Aggregate Agg Types" property of pivot table viewer to "avg"
    And user sets "Pivot Column Names" property of pivot table viewer to "SEVERITY"
    Then the "aggregated rows" reading of pivot table viewer should be 6
    And no errors should have been logged

  Scenario: A saved configuration is offered only while the table still has its columns
    When user sets "Group By Column Names" property of pivot table viewer to "RACE"
    And user sets "Aggregate Column Names" property of pivot table viewer to "WEIGHT"
    And user sets "Aggregate Agg Types" property of pivot table viewer to "avg"
    And user sets "Pivot Column Names" property of pivot table viewer to ""
    And user picks "Save parameters" from the history menu of pivot table viewer
    Then the "history entries" reading of pivot table viewer should contain "avg(WEIGHT)"
    When user sets "Group By Column Names" property of pivot table viewer to "DIS_POP"
    And user sets "Aggregate Column Names" property of pivot table viewer to "AGE"
    And user sets "Aggregate Agg Types" property of pivot table viewer to "avg"
    And user sets "Pivot Column Names" property of pivot table viewer to "SEVERITY"
    And user renames "WEIGHT" column to "MASS"
    And user clicks on close icon of pivot table viewer
    And user adds a pivot table viewer
    Then the "aggregated rows" reading of pivot table viewer should be 6
    And the "history entries" reading of pivot table viewer should not contain "WEIGHT"
    When user renames "MASS" column to "WEIGHT"
    And user clicks on close icon of pivot table viewer
    And user adds a pivot table viewer
    Then the "history entries" reading of pivot table viewer should contain "avg(WEIGHT)"
    And the "aggregated rows" reading of pivot table viewer should be 6
    When user clears the saved pivot table parameters
    Then no errors should have been logged
