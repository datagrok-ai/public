@journey @viewers @realizes:viewers.pivot-table
Feature: Pivot table — the frame around the aggregation
  What a fresh pivot table configures itself to, what its counts area says about the result, and the
  chrome that can be taken away: the tag rows and the counts under Show Header, the history and
  refresh icons under Show Command Bar, the title and the description. Closing the viewer and adding
  it again brings the same automatic cross tab back with a clean console (GROK-17122), and the
  Select columns dialog of the Data row leaves the row with the one table it had (github-3414,
  GROK-14995). One journey on demog-1000: `auto()` takes the two categorical columns with the most
  categories, descending — DIS_POP (6) for Group by, SEVERITY (5) for Pivot — and the first
  numerical column, AGE, with avg, so the cross tab is 6 rows by 6 columns and RA's None cell reads
  52.30 while AS has no Critical row at all.
  Not translated: the layout `find`-by-id echoes of the old spec, and its `expect(consoleErrors)`
  deltas, which are the `no errors should have been logged` rider every scenario here ends on.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pivot table viewer
    Then the "group by" reading of pivot table viewer should be "DIS_POP"
    And the "aggregate" reading of pivot table viewer should be "avg(AGE)"
    And the "pivot" reading of pivot table viewer should be "SEVERITY"

  Scenario: A fresh viewer configures itself from the column types
    Then the "key columns" reading of pivot table viewer should be "DIS_POP"
    And the "aggregations" reading of pivot table viewer should be "avg(AGE)"
    And the "data" reading of pivot table viewer should be "demog-1000"
    And the "default aggregation" reading of pivot table viewer should be "avg"
    And the "error" reading of pivot table viewer should be ""
    And pivot table viewer should have a "group by chip DIS_POP" area
    And pivot table viewer should have a "aggregate chip avg(AGE)" area
    And pivot table viewer should have a "pivot chip SEVERITY" area
    And "Group By Column Names" property of pivot table viewer should be "DIS_POP"
    And "Aggregate Column Names" property of pivot table viewer should be "AGE"
    And "Aggregate Agg Types" property of pivot table viewer should be "avg"
    And "Pivot Column Names" property of pivot table viewer should be "SEVERITY"
    And no errors should have been logged

  Scenario: The counts area reports the shape of the aggregation, and the cells hold it
    Then the "aggregated rows" reading of pivot table viewer should be 6
    And the "aggregated columns" reading of pivot table viewer should be 6
    And the "limited" reading of pivot table viewer should be "false"
    And the "rows shown" reading of pivot table viewer should be 1000
    And pivot table viewer should have a "counts" area
    And pivot table viewer should have a "add to workspace" area
    And the "text of grid cell 5 of DIS_POP" reading of pivot table viewer should be "RA"
    And the "text of grid cell 5 of None avg(AGE)" reading of pivot table viewer should be "52.30"
    And the "text of grid cell 2 of Critical avg(AGE)" reading of pivot table viewer should be "29.00"
    And the "text of grid cell 1 of Critical avg(AGE)" reading of pivot table viewer should be ""
    And the aggregated values of pivot table viewer should match "avg(AGE)" grouped by "DIS_POP" pivoted on "SEVERITY"
    And no errors should have been logged

  Scenario: Show Header takes the tag rows and the counts away and gives them back
    When user sets "Show Header" property of pivot table viewer to "false"
    Then the "header shown" reading of pivot table viewer should be "false"
    And pivot table viewer should not have a "group by row" area
    And pivot table viewer should not have a "aggregate row" area
    And pivot table viewer should not have a "pivot row" area
    And pivot table viewer should not have a "counts" area
    And the "aggregated rows" reading of pivot table viewer should be 6
    And the "text of grid cell 5 of DIS_POP" reading of pivot table viewer should be "RA"
    When user sets "Show Header" property of pivot table viewer to "true"
    Then the "header shown" reading of pivot table viewer should be "true"
    And pivot table viewer should have a "group by row" area
    And pivot table viewer should have a "counts" area
    And no errors should have been logged

  Scenario: Show Command Bar takes the history and refresh icons away and gives them back
    Then pivot table viewer should have a "history" area
    And pivot table viewer should have a "refresh" area
    When user sets "Show Command Bar" property of pivot table viewer to "false"
    Then the "command bar shown" reading of pivot table viewer should be "false"
    And pivot table viewer should not have a "command bar" area
    And pivot table viewer should not have a "history" area
    And pivot table viewer should not have a "refresh" area
    And pivot table viewer should have a "group by row" area
    When user sets "Show Command Bar" property of pivot table viewer to "true"
    Then the "command bar shown" reading of pivot table viewer should be "true"
    And pivot table viewer should have a "history" area
    And no errors should have been logged

  Scenario: The title bar shows the title and the description obeys its visibility mode
    When user sets properties of pivot table viewer:
      | Show Title | true      |
      | Title      | Cross tab |
    Then title of pivot table viewer should have text "Cross tab"
    When user sets properties of pivot table viewer:
      | Description                 | Rows by disease |
      | Description Visibility Mode | Always          |
    Then description of pivot table viewer should be visible
    And description of pivot table viewer should have text "Rows by disease"
    When user sets "Description Visibility Mode" property of pivot table viewer to "Never"
    Then description of pivot table viewer should be hidden
    When user sets properties of pivot table viewer:
      | Title                       | |
      | Description                 | |
      | Description Visibility Mode | Auto |
    Then the "aggregated rows" reading of pivot table viewer should be 6
    And no errors should have been logged

  Scenario: Closing the viewer and adding it again brings the same cross tab back
    When user clicks on close icon of pivot table viewer
    Then pivot table viewer should be absent
    When user adds a pivot table viewer
    Then pivot table viewer should be visible
    And the "group by" reading of pivot table viewer should be "DIS_POP"
    And the "aggregate" reading of pivot table viewer should be "avg(AGE)"
    And the "pivot" reading of pivot table viewer should be "SEVERITY"
    And the "aggregated rows" reading of pivot table viewer should be 6
    And the "text of grid cell 5 of None avg(AGE)" reading of pivot table viewer should be "52.30"
    And no errors should have been logged

  Scenario: The Select columns dialog of the Data row leaves the row as it was
    When user clicks on the "data chip demog-1000" area of pivot table viewer
    Then "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then "Select columns..." dialog should be absent
    And the "data" reading of pivot table viewer should be "demog-1000"
    And the "group by" reading of pivot table viewer should be "DIS_POP"
    And the "aggregated rows" reading of pivot table viewer should be 6
    And the "aggregated columns" reading of pivot table viewer should be 6
    And no errors should have been logged
