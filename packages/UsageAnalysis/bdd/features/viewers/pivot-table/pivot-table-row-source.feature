@journey @viewers @realizes:viewers.pivot-table
Feature: Pivot table — what Row Source aggregates
  Row Source decides which rows the aggregation runs over: at All the pivot ignores a filter on the
  source table, at Filtered it re-aggregates over the rows the filter passes and says so through its
  own "rows shown" reading. The filter is a categorical card of the filter panel, a registered filter
  — a bitset written straight into the table does not survive the pivot's next refresh, which asks
  the platform to recompute the filter. One journey on demog-1000 grouped by DIS_POP with avg(AGE):
  447 of the 1000 rows are male, RA averages 51.60 over all of them and 52.10 over the male ones.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pivot table viewer
    And user sets "Pivot Column Names" property of pivot table viewer to ""
    And user sets "Row Source" property of pivot table viewer to "All"
    Then the "aggregated rows" reading of pivot table viewer should be 6
    And the "rows shown" reading of pivot table viewer should be 1000

  Scenario: At Row Source All the source filter changes nothing
    When user adds a categorical filter on "SEX" keeping "M"
    Then 447 rows should pass the filter
    And the "rows shown" reading of pivot table viewer should be 1000
    And the "aggregated rows" reading of pivot table viewer should be 6
    And the "text of grid cell 5 of avg(AGE)" reading of pivot table viewer should be "51.60"
    And the aggregated values of pivot table viewer should match "avg(AGE)" grouped by "DIS_POP"
    And no errors should have been logged

  Scenario: At Row Source Filtered the pivot re-aggregates over the filtered rows
    When user sets "Row Source" property of pivot table viewer to "Filtered"
    Then the "rows shown" reading of pivot table viewer should be 447
    And pivot table viewer should show 447 rows
    And the "aggregated rows" reading of pivot table viewer should be 6
    And the "text of grid cell 5 of avg(AGE)" reading of pivot table viewer should be "52.10"
    And the aggregated values of pivot table viewer should match "avg(AGE)" grouped by "DIS_POP" over the filtered rows
    And no errors should have been logged

  Scenario: Lifting the filter takes the pivot back to the whole table
    When user adds a categorical filter on "SEX" keeping "F, M"
    Then all rows should pass the filter
    And the "rows shown" reading of pivot table viewer should be 1000
    And the "text of grid cell 5 of avg(AGE)" reading of pivot table viewer should be "51.60"
    And the aggregated values of pivot table viewer should match "avg(AGE)" grouped by "DIS_POP"
    When user sets "Row Source" property of pivot table viewer to "All"
    Then the "rows shown" reading of pivot table viewer should be 1000
    And no errors should have been logged
