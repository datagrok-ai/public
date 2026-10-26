@journey @viewers @realizes:viewers.network-diagram
Feature: Network diagram regressions of closed issues
  Closed network diagram issues, each claimed by what the viewer draws and reports. core#869: Node
  Shape = icon drew nothing and logged an error for every node, since vis's icon shape needs a glyph
  the viewer did not give; it now draws all four nodes of demog-1000 (F, M, false, true), each with
  ink in its box, and logs nothing. public#2129: under Auto Layout a viewer smaller than 200 by 200 px
  drops its two on-canvas column selectors, and Auto Layout off keeps them at the same size.
  Fixture: demog-1000, Node 1 = SEX and Node 2 = CONTROL picked automatically.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a network diagram viewer
    Then the "nodes" reading of network diagram viewer should be 4
    And network diagram viewer should be painted

  Scenario: The icon node shape draws every node and logs nothing per node (core#869)
    When user sets "nodeShape" property of network diagram viewer to "icon"
    Then network diagram viewer should have repainted
    And the "nodes" reading of network diagram viewer should be 4
    And the "node names" reading of network diagram viewer should contain "F"
    And the "node \"F\"" area of network diagram viewer should be painted
    And the "node \"true\"" area of network diagram viewer should be painted
    And network diagram viewer should report no error
    And no errors should have been logged
    When user sets "nodeShape" property of network diagram viewer to "dot"
    Then the "nodes" reading of network diagram viewer should be 4

  Scenario: Auto Layout drops the column selectors on a small diagram, and Auto Layout off keeps them (public#2129)
    Then network diagram viewer should have a "node 1 selector" area
    And network diagram viewer should have a "node 2 selector" area
    When user resizes network diagram viewer to 180 by 180
    Then network diagram viewer should not have a "node 1 selector" area
    And network diagram viewer should not have a "node 2 selector" area
    When user sets "autoLayout" property of network diagram viewer to "false"
    Then network diagram viewer should have a "node 1 selector" area
    And network diagram viewer should have a "node 2 selector" area
    When user sets "autoLayout" property of network diagram viewer to "true"
    Then network diagram viewer should not have a "node 1 selector" area
    When user restores the size of network diagram viewer
    Then network diagram viewer should have a "node 1 selector" area
    And no errors should have been logged
