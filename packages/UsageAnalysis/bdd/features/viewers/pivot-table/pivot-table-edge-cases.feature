@viewers @edge-cases @realizes:viewers.pivot-table
Feature: Pivot table edge cases — its Group by column removed, an empty row set, every row deleted
  The pivot table on data changes its other features do not reach; the pivot is outside the
  crash-only Data Manipulation test the other viewers get (`viewers_test.dart`, the first 30
  viewers). A Group by column removed from the table while the pivot shows it is replaced the way
  the pivot picks its columns when it opens — the categorical column with the most categories,
  between 1 and 100 of them (`PivotViewerLook.auto`, run on every column removal) — and the pivot
  aggregates again with no error. A table filtered down to nothing and a table whose rows are all
  deleted leave a pivot with no aggregated row and no error. demog-1000 opens the pivot grouped
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

  Scenario: With every row filtered out the pivot holds no aggregated row and logs no error
    When user filters rows where "AGE" is between 200 and 300
    Then 0 rows should pass the filter
    And the "rows shown" reading of pivot table viewer should be 0
    And the "aggregated rows" reading of pivot table viewer should be 0
    And the "error" reading of pivot table viewer should be ""
    And no errors should have been logged

  Scenario: Once every row of the table is deleted the pivot holds no aggregated row and logs no error
    When user selects all rows
    And user deletes the selected rows
    Then the table should have 0 rows
    And the "aggregated rows" reading of pivot table viewer should be 0
    And the "error" reading of pivot table viewer should be ""
    And no errors should have been logged
