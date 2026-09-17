@journey @viewers @realizes:viewers.pivot-table
Feature: Pivot table — the filter and the selection a click sends back
  What the aggregated grid does to the table it came from. With Row Source = All and Filtering
  Enabled, moving the current cell filters the source table down to the clicked group and says so in
  the filter's own message, and a second click REPLACES that filter rather than intersecting with it
  (GROK-17726); under any other Row Source, or with Filtering Enabled off, a click leaves the filter
  where it was; and a row selected in the aggregated grid selects that group's source rows through
  the link the aggregation attaches. What Row Source itself aggregates is the sibling feature
  `pivot-table-row-source.feature`, which needs a registered filter card and a view of its own.
  One journey on demog-1000 grouped by DIS_POP with avg(AGE) and no pivot: RA holds 434 rows,
  Psoriasis 204, and the two of them together 638.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pivot table viewer
    And user sets "Pivot Column Names" property of pivot table viewer to ""
    And user sets "Row Source" property of pivot table viewer to "All"
    Then the "aggregated rows" reading of pivot table viewer should be 6
    And the "rows shown" reading of pivot table viewer should be 1000
    And the "text of grid cell 5 of DIS_POP" reading of pivot table viewer should be "RA"

  Scenario: Selecting a row of the aggregated grid selects the group's source rows
    Given no rows should be selected
    When user clicks on the "grid row header 5" area of pivot table viewer holding Control
    Then 434 rows should be selected
    And only rows where "DIS_POP" is "RA" should be selected
    And all rows should pass the filter
    When user clicks on the "grid row header 5" area of pivot table viewer holding Control
    Then no rows should be selected
    And the "text of grid cell 6 of DIS_POP" reading of pivot table viewer should be "UC"
    And no errors should have been logged

  Scenario: A second selected row adds its group to the selection
    When user clicks on the "grid row header 4" area of pivot table viewer holding Control
    Then 204 rows should be selected
    And only rows where "DIS_POP" is "Psoriasis" should be selected
    When user clicks on the "grid row header 5" area of pivot table viewer holding Control
    Then 638 rows should be selected
    And all rows where "DIS_POP" is "RA" should be selected
    And all rows where "DIS_POP" is "Psoriasis" should be selected
    And no rows where "DIS_POP" is "UC" should be selected
    When user clicks on the "grid row header 4" area of pivot table viewer holding Control
    And user clicks on the "grid row header 5" area of pivot table viewer holding Control
    Then no rows should be selected
    And the "text of grid cell 6 of DIS_POP" reading of pivot table viewer should be "UC"
    And no errors should have been logged

  Scenario: A cell click filters the source table down to the clicked group
    Then "Filtering Enabled" property of pivot table viewer should be "true"
    And the "filter label" reading of pivot table viewer should be ""
    When user clicks on the "grid cell 5 of DIS_POP" area of pivot table viewer
    Then 434 rows should pass the filter
    And all rows where "DIS_POP" is "RA" should pass the filter
    And no rows where "DIS_POP" is "UC" should pass the filter
    And the "filter label" reading of pivot table viewer should be "DIS_POP in [RA]"
    And the "rows shown" reading of pivot table viewer should be 1000
    When user sets "Filtering Enabled" property of pivot table viewer to "false"
    And user resets the filter
    Then all rows should pass the filter
    And the "text of grid cell 6 of DIS_POP" reading of pivot table viewer should be "UC"
    And no errors should have been logged

  Scenario: The next click replaces the filter instead of narrowing it further
    When user sets "Filtering Enabled" property of pivot table viewer to "true"
    And user clicks on the "grid cell 5 of DIS_POP" area of pivot table viewer
    Then 434 rows should pass the filter
    When user clicks on the "grid cell 4 of DIS_POP" area of pivot table viewer
    Then 204 rows should pass the filter
    And all rows where "DIS_POP" is "Psoriasis" should pass the filter
    And no rows where "DIS_POP" is "RA" should pass the filter
    And the "filter label" reading of pivot table viewer should be "DIS_POP in [Psoriasis]"
    When user sets "Filtering Enabled" property of pivot table viewer to "false"
    And user resets the filter
    Then all rows should pass the filter
    And the "text of grid cell 6 of DIS_POP" reading of pivot table viewer should be "UC"
    And no errors should have been logged

  Scenario: With Filtering Enabled off a click leaves the filter alone
    Then "Filtering Enabled" property of pivot table viewer should be "false"
    When user clicks on the "grid cell 5 of DIS_POP" area of pivot table viewer
    Then all rows should pass the filter
    And the "filter label" reading of pivot table viewer should be ""
    When user clicks on the "grid cell 4 of DIS_POP" area of pivot table viewer
    Then all rows should pass the filter
    And the "text of grid cell 6 of DIS_POP" reading of pivot table viewer should be "UC"
    And no errors should have been logged

  Scenario: At Row Source Filtered a click leaves the filter alone as well
    When user sets "Row Source" property of pivot table viewer to "Filtered"
    And user sets "Filtering Enabled" property of pivot table viewer to "true"
    And user clicks on the "grid cell 5 of DIS_POP" area of pivot table viewer
    Then all rows should pass the filter
    And the "filter label" reading of pivot table viewer should be ""
    And the "rows shown" reading of pivot table viewer should be 1000
    When user clicks on the "grid cell 2 of DIS_POP" area of pivot table viewer
    Then all rows should pass the filter
    When user sets "Row Source" property of pivot table viewer to "All"
    Then no errors should have been logged
