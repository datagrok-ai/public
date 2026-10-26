@viewers @edge-cases @realizes:viewers.network-diagram
Feature: Network diagram edge cases — an empty row set, every row deleted, one column for both nodes
  The network diagram on settings and data changes its other features do not reach. A table
  filtered down to nothing hides every node and a table whose rows are all deleted logs no error
  (core#1041 lists the 0-row table among the clean cases). The first crash of the open audit
  core#1041 is written as the behaviour it breaks: with the same column as Node 1 and Node 2 a
  filter change should raise no error. demog-1000 opens the diagram on SEX and CONTROL, 4 nodes
  and 4 edges.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a network diagram viewer
    Then the "nodes" reading of network diagram viewer should be 4
    And the "edges" reading of network diagram viewer should be 4

  Scenario: With every row filtered out the diagram hides every node and logs no error
    When user filters rows where "AGE" is between 200 and 300
    Then 0 rows should pass the filter
    And the "rows shown" reading of network diagram viewer should be 0
    And the "nodes" reading of network diagram viewer should be 0
    And no errors should have been logged
    When user resets the filter
    Then the "nodes" reading of network diagram viewer should be 4

  Scenario: Once every row of the table is deleted the diagram logs no error
    When user selects all rows
    And user deletes the selected rows
    Then the table should have 0 rows
    And the "rows shown" reading of network diagram viewer should be 0
    And network diagram viewer should report no error
    And no errors should have been logged

  @known-failure:41 @realizes:core#1041
  Scenario: With the same column as Node 1 and Node 2 a filter change raises no error (core#1041)
    A filter change logs "Column named 'SEX' already exists" and the graph stops following the
    filter (core#1041, crash 1). Remove the tag once the diagram follows the filter.
    When user sets "node2ColumnName" property of network diagram viewer to "SEX"
    Then the "node 2 column" reading of network diagram viewer should be "SEX"
    When user filters rows where "SEX" is "F"
    Then no errors should have been logged
