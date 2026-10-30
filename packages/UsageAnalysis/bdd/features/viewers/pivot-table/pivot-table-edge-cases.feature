@viewers @edge-cases @realizes:viewers.pivot-table
Feature: Pivot table edge cases — its Group by column removed, ADD with nothing passing the filter
  The pivot table on data changes its other features do not reach; the pivot is outside the
  crash-only Data Manipulation test the other viewers get (`viewers_test.dart`, the first 30
  viewers). A Group by column removed from the table while the pivot shows it is replaced the way
  the pivot picks its columns when it opens — the categorical column with the most categories,
  between 1 and 100 of them (`PivotViewerLook.auto`, run on every column removal) — and the pivot
  aggregates again with no error. Filtered down to nothing or with every row
  deleted: `viewer-rules-edge-cases.feature`. demog-1000 opens the pivot grouped
  by DIS_POP, pivoted on SEVERITY, aggregating avg(AGE).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pivot table viewer
    Then the "group by" reading of pivot table viewer should be "DIS_POP"
    And the "aggregate" reading of pivot table viewer should be "avg(AGE)"

  Scenario: The Group by column removed from the table is replaced and the pivot aggregates again
    When user removes "DIS_POP" column
    Then the table should not have a column "DIS_POP"
    And the "group by" reading of pivot table viewer should be "SEVERITY"
    And the "error" reading of pivot table viewer should be ""
    And the "aggregated rows" reading of pivot table viewer should be 5
    And no errors should have been logged

  Scenario: ADD with every row filtered out adds an empty table
    When user filters rows where "AGE" is between 200 and 300
    Then 0 rows should pass the filter
    When user clicks on the "add to workspace" area of pivot table viewer
    Then table "demog-1000 aggregation" should be open
    And table "demog-1000 aggregation" should have 0 rows
    And no errors should have been logged
