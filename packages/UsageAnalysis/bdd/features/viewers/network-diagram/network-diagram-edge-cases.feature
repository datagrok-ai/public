@viewers @edge-cases @realizes:viewers.network-diagram
Feature: Network diagram edge cases — every row deleted, one column for both nodes
  The network diagram on settings and data changes its other features do not reach (a table
  filtered down to nothing: `viewer-rules-edge-cases.feature`). A table whose rows are all deleted
  logs no error (core#1041 lists the 0-row table among the clean cases). The first crash of the open audit
  core#1041 is written as the behaviour it breaks: with the same column as Node 1 and Node 2 a
  filter change should raise no error. demog-1000 opens the diagram on SEX and CONTROL, 4 nodes
  and 4 edges.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a network diagram viewer
    Then the "nodes" reading of network diagram viewer should be 4
    And the "edges" reading of network diagram viewer should be 4

  Scenario: Once every row of the table is deleted the diagram logs no error
    When user selects all rows
    And user deletes the selected rows
    Then the table should have 0 rows
    And the "rows shown" reading of network diagram viewer should be 0
    And network diagram viewer should report no error
    And no errors should have been logged

  @realizes:core#1041
  Scenario: With the same column as Node 1 and Node 2 a filter change raises no error (core#1041)
    Guards core#1041 (crash 1, fixed): a filter change logged "Column named 'SEX' already exists"
    and the graph stopped following the filter.
    When user sets "node2ColumnName" property of network diagram viewer to "SEX"
    Then the "node 2 column" reading of network diagram viewer should be "SEX"
    When user filters rows where "SEX" is "F"
    Then no errors should have been logged

  Scenario: More than 10000 nodes are not drawn and the diagram says so
    demog (5850 rows) with USUBJID as Node 1 and a distinct PARTNER per row as Node 2 makes 11700
    nodes: no table here has more than 10000 categories in one column, so the two columns together
    pass the limit.
    Given user opens demog dataset
    Then the table should have 5850 rows
    When user adds a calculated column "PARTNER" with formula "${USUBJID} + '-b'"
    Then every value of "PARTNER" column should contain "-b"
    When user adds a network diagram viewer with:
      | node1ColumnName | USUBJID |
      | node2ColumnName | PARTNER |
    Then network diagram viewer should contain the text "Too many nodes to display, max 10000"
    And the "nodes" reading of network diagram viewer should be 0
    And no errors should have been logged
