@journey @viewers @realizes:viewers.network-diagram
Feature: Clicking a network diagram node selects the rows behind it
  A click on a node selects the rows of every edge hanging off it, and the viewer says in advance
  how many that is: `rows of node "<label>"` is `getIdsByGroupIds(getConnectedEdges(id))`, the same
  call the click makes.
  This is where the old spec spent most of its lines. It read the canvas back as pixels, threw away
  anything near-white or low-saturation, bucketed the rest into 40x40 cells, took the eight densest
  and clicked up to six of them until the selection stopped being zero — then asserted only that
  something had been selected. Whether it had hit F, M or an edge, and whether the number was the
  right one, it could not say. Every node vis is drawing is now a hit area named after its label,
  so the click is aimed by name and the count is checked against the node's own.
  The same pixel hunt made the negative meaningless: **Select Rows On Click** off was proved by
  four clicks at those guessed points leaving the selection at zero, which passes whether or not
  any of them landed on a node. Here the click is on `node "F"` — the very node that selected 553
  rows in the scenario above — so the zero means the click was ignored, not that it missed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a network diagram viewer
    And user clears the row selection
    Then the "nodes" reading of network diagram viewer should be 4
    And the "rows of node \"F\"" reading of network diagram viewer should be 553
    And no rows should be selected
    And network diagram viewer should report no error

  Scenario: Clicking a node selects exactly the rows it stands for
    When user clicks on the "node \"F\"" area of network diagram viewer
    Then 553 rows should be selected
    And only rows where "SEX" is "F" should be selected
    And network diagram viewer should have repainted
    When user clears the row selection
    And user clicks on the "node \"M\"" area of network diagram viewer
    Then 447 rows should be selected
    And only rows where "SEX" is "M" should be selected
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged

  Scenario: A node of the second column selects the rows of that column's value
    When user clicks on the "node \"true\"" area of network diagram viewer
    Then 6 rows should be selected
    And only rows where "CONTROL" is "true" should be selected
    When user clears the row selection
    And user clicks on the "node \"false\"" area of network diagram viewer
    Then 994 rows should be selected
    And only rows where "CONTROL" is "false" should be selected
    When user clears the row selection
    Then no errors should have been logged

  Scenario: A node click is announced as an event
    Given user listens for "d4-network-diagram-node-click" event on network diagram viewer
    When user clicks on the "node \"M\"" area of network diagram viewer
    Then "d4-network-diagram-node-click" event should have fired on network diagram viewer
    And 447 rows should be selected
    When user clears the row selection
    Then no errors should have been logged

  Scenario: With Select Rows On Click off the same click selects nothing
    When user clicks on the "node \"F\"" area of network diagram viewer
    Then 553 rows should be selected
    When user clears the row selection
    And user sets properties of network diagram viewer:
      | selectRowsOnClick  | false |
      | selectEdgesOnClick | false |
    And user clicks on the "node \"F\"" area of network diagram viewer
    Then no rows should be selected
    When user sets properties of network diagram viewer:
      | selectRowsOnClick  | true |
      | selectEdgesOnClick | true |
    And user clicks on the "node \"F\"" area of network diagram viewer
    Then 553 rows should be selected
    When user clears the row selection
    Then no errors should have been logged

  Scenario: A filtered-away node cannot be clicked because it is not there
    When user adds a categorical filter on "SEX" keeping "F"
    Then 553 rows should pass the filter
    And network diagram viewer should not have a "node \"M\"" area
    When user clicks on the "node \"F\"" area of network diagram viewer
    Then 553 rows should be selected
    And every selected row should pass the filter
    When user clears the row selection
    And user hovers over "SEX" filter card
    And user clicks on close of "SEX" filter card
    Then 1000 rows should pass the filter
    And network diagram viewer should have a "node \"M\"" area
    And no errors should have been logged
