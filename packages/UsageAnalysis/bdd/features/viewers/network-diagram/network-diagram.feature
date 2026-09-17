@journey @viewers @realizes:viewers.network-diagram @realizes:entities.viewer.action.close-viewer
Feature: Network diagram graph shape, node columns, filtering and chrome
  Which nodes and edges vis.js is drawing, what the two node columns are, and what the filter and
  the chrome flags do to them.
  Everything here is a count or a name the viewer reports. The old spec proved the viewer had drawn
  by counting more than a thousand non-white canvas pixels, read the node columns by scraping the
  on-canvas combo box's `innerText`, and proved **Show Filtered Out Nodes** had worked by asserting
  the canvas held more ink than before — which is true of any repaint that adds anything. `nodes`
  is how many vis is drawing with the hidden ones excluded, so bringing the filtered-away nodes
  back is 3 going to 4, and a node the filter hid has no `node "<label>"` area at all.
  There is no wait anywhere in this feature. vis lays the graph out on its own animation loop, and
  `isRenderPending` now covers it, so a step settles on this viewer the way it settles on every
  other one instead of chaining canvas-quiet polls — the old spec had nine.
  **Show Arrows is the one scenario with no honest observer**, and it is written as what it is.
  Nothing on this viewer reports an arrow head: the flag goes into vis's edge options and vis draws
  it, so the only things that can be claimed are that the viewer reports the setting, that the
  canvas repainted when it changed, and that a rebuild of the graph keeps it — which is what
  GROK-20617 was about. The old spec asserted the repaint alone and called it "the arrow heads are
  drawn"; the repaint is kept here next to the two claims that do mean something, and no more is
  read into it. An `arrow heads` count, or arrow geometry among the hit areas, is what would turn
  this into a real claim.
  Fixture: demog-1000 with Node 1 = SEX and Node 2 = CONTROL picked automatically — 4 nodes
  (F, M, false, true) and 4 edges, one per distinct pair.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a network diagram viewer
    Then 1000 rows should pass the filter
    And the "node 1 column" reading of network diagram viewer should be "SEX"
    And the "node 2 column" reading of network diagram viewer should be "CONTROL"
    And the "nodes" reading of network diagram viewer should be 4
    And the "edges" reading of network diagram viewer should be 4
    And the "rows shown" reading of network diagram viewer should be 1000
    And network diagram viewer should be painted
    And network diagram viewer should report no error

  Scenario: Each node stands for the rows of the edges hanging off it
    Then the "node names" reading of network diagram viewer should contain "F"
    And the "node names" reading of network diagram viewer should contain "M"
    And the "rows of node \"F\"" reading of network diagram viewer should be 553
    And the "rows of node \"M\"" reading of network diagram viewer should be 447
    And the "rows of node \"false\"" reading of network diagram viewer should be 994
    And the "rows of node \"true\"" reading of network diagram viewer should be 6
    And network diagram viewer should have a "node \"F\"" area
    And network diagram viewer should have a "node \"M\"" area
    And no errors should have been logged

  Scenario: Changing Node 1 rebuilds the graph around the other column's values
    When user sets "node1ColumnName" property of network diagram viewer to "RACE"
    Then the "node 1 column" reading of network diagram viewer should be "RACE"
    And the "nodes" reading of network diagram viewer should be 6
    And the "edges" reading of network diagram viewer should be 5
    And the "node names" reading of network diagram viewer should contain "Caucasian"
    And the "node names" reading of network diagram viewer should not contain "F"
    And the "rows of node \"Caucasian\"" reading of network diagram viewer should be 896
    And the "rows of node \"Other\"" reading of network diagram viewer should be 62
    And the "rows of node \"Black\"" reading of network diagram viewer should be 27
    And the "rows of node \"Asian\"" reading of network diagram viewer should be 15
    And network diagram viewer should have a "node \"Caucasian\"" area
    When user sets "node1ColumnName" property of network diagram viewer to "SEX"
    Then the "nodes" reading of network diagram viewer should be 4
    And no errors should have been logged

  Scenario: A filtered-away node is not drawn, and Show Filtered Out Nodes brings it back
    When user adds a categorical filter on "SEX" keeping "F"
    Then 553 rows should pass the filter
    And the "rows shown" reading of network diagram viewer should be 553
    And the "nodes" reading of network diagram viewer should be 3
    And the "edges" reading of network diagram viewer should be 2
    And the "node names" reading of network diagram viewer should not contain "M"
    And network diagram viewer should not have a "node \"M\"" area
    And network diagram viewer should have a "node \"F\"" area
    When user sets "showFilteredOutNodes" property of network diagram viewer to "true"
    Then the "filtered out nodes shown" reading of network diagram viewer should be "true"
    And the "nodes" reading of network diagram viewer should be 4
    And the "edges" reading of network diagram viewer should be 4
    And network diagram viewer should have a "node \"M\"" area
    And the "rows shown" reading of network diagram viewer should be 553
    When user sets "showFilteredOutNodes" property of network diagram viewer to "false"
    Then the "nodes" reading of network diagram viewer should be 3
    When user hovers over "SEX" filter card
    And user clicks on close of "SEX" filter card
    Then 1000 rows should pass the filter
    And the "nodes" reading of network diagram viewer should be 4
    And the "rows shown" reading of network diagram viewer should be 1000
    And no errors should have been logged

  Scenario: Show Column Selectors takes the two on-viewer selectors away
    Then network diagram viewer should have a "node 1 selector" area
    And network diagram viewer should have a "node 2 selector" area
    When user sets "showColumnSelectors" property of network diagram viewer to "false"
    Then network diagram viewer should not have a "node 1 selector" area
    And network diagram viewer should not have a "node 2 selector" area
    And the "node 1 column" reading of network diagram viewer should be "SEX"
    And network diagram viewer should have a "node \"F\"" area
    When user sets "showColumnSelectors" property of network diagram viewer to "true"
    Then network diagram viewer should have a "node 1 selector" area
    And no errors should have been logged

  Scenario: Show Arrows is a setting the viewer reports, and a rebuild keeps it
    Then the "arrows" reading of network diagram viewer should be "none"
    When user sets "showArrows" property of network diagram viewer to "to"
    Then the "arrows" reading of network diagram viewer should be "to"
    And network diagram viewer should have repainted
    When user sets "node1ColumnName" property of network diagram viewer to "RACE"
    Then the "nodes" reading of network diagram viewer should be 6
    And the "arrows" reading of network diagram viewer should be "to"
    When user sets properties of network diagram viewer:
      | node1ColumnName | SEX  |
      | showArrows      | none |
    Then the "arrows" reading of network diagram viewer should be "none"
    And the "nodes" reading of network diagram viewer should be 4
    And no errors should have been logged

  Scenario: The title bar closes the network diagram
    When user clicks on close icon of network diagram viewer
    Then network diagram viewer should be absent
    And the open tableview should have 0 network diagram viewers
    And no errors should have been logged
