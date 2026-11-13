@viewers @realizes:viewers.network-diagram
Feature: Network diagram — Node 1 picked on the viewer's own column selector
  A person changes the first node column with the combo box drawn on the viewer, not through the
  property panel (`network-diagram.feature` sets the property). The pick has to reach the look and
  rebuild the graph: on demog-1000, Node 1 = RACE and Node 2 = CONTROL gives 6 nodes (4 races, 2
  CONTROL values) and 5 edges, and the Caucasian node stands for its 896 rows.

  Scenario: Picking RACE in the Node 1 selector rebuilds the graph around the races, and SEX puts it back
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a network diagram viewer
    Then the "node 1 column" reading of network diagram viewer should be "SEX"
    And the "nodes" reading of network diagram viewer should be 4
    And network diagram viewer should have a "node 1 selector" area
    When user picks "RACE" in the column selector at the "node 1 selector" area of network diagram viewer
    Then the "node 1 column" reading of network diagram viewer should be "RACE"
    And "node1ColumnName" property of network diagram viewer should be "RACE"
    And the "node 2 column" reading of network diagram viewer should be "CONTROL"
    And the "nodes" reading of network diagram viewer should be 6
    And the "edges" reading of network diagram viewer should be 5
    And the "node names" reading of network diagram viewer should contain "Caucasian"
    And the "node names" reading of network diagram viewer should not contain "F"
    And the "rows of node \"Caucasian\"" reading of network diagram viewer should be 896
    And network diagram viewer should have a "node \"Caucasian\"" area
    When user picks "SEX" in the column selector at the "node 1 selector" area of network diagram viewer
    Then the "node 1 column" reading of network diagram viewer should be "SEX"
    And the "nodes" reading of network diagram viewer should be 4
    And the "rows of node \"F\"" reading of network diagram viewer should be 553
    And no errors should have been logged
