@journey @viewers @realizes:viewers.network-diagram
Feature: Suspend Simulation freezes the layout, and the Layout menu says so too
  vis lays the graph out on a physics loop, and `layout signature` is a hash of every drawn node's
  box in **canvas** coordinates — vis's world space. Pan and zoom are viewport moves and leave it
  alone; only the simulation moves it.
  That reading is the whole point of this feature. The old spec's "Suspend simulation freezes the
  layout" scenario checked the property's checkbox and then waited for the canvas to go quiet — it
  asserted nothing whatsoever, so it passed with the setting wired to nothing at all. "Remember it,
  act, assert it is as remembered" is the claim; the acting is a click that repaints the viewer,
  so a signature that moved would have moved for a reason.
  The converse is asserted next to it: with the simulation running, a rebuild of the graph — a new
  Node 1 column — lands the nodes somewhere else and the signature does move. Without that half,
  "unchanged" would be satisfied by a signature that never changes at all.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a network diagram viewer
    And user clears the row selection
    Then the "nodes" reading of network diagram viewer should be 4
    And "suspendSimulation" property of network diagram viewer should be "false"
    And network diagram viewer should report no error

  Scenario: With the simulation suspended the layout does not move
    When user sets "suspendSimulation" property of network diagram viewer to "true"
    And user remembers the "layout signature" reading of network diagram viewer
    And user clicks on the "node \"F\"" area of network diagram viewer
    Then 553 rows should be selected
    And network diagram viewer should have repainted
    And the "layout signature" reading of network diagram viewer should be as remembered
    When user clears the row selection
    Then the "layout signature" reading of network diagram viewer should be as remembered
    When user sets "suspendSimulation" property of network diagram viewer to "false"
    Then "suspendSimulation" property of network diagram viewer should be "false"
    And no errors should have been logged

  Scenario: A rebuild of the graph moves the layout, so an unchanged signature means something
    When user remembers the "layout signature" reading of network diagram viewer
    And user sets "node1ColumnName" property of network diagram viewer to "RACE"
    Then the "nodes" reading of network diagram viewer should be 6
    And the "layout signature" reading of network diagram viewer should not be as remembered
    When user sets "node1ColumnName" property of network diagram viewer to "SEX"
    Then the "nodes" reading of network diagram viewer should be 4
    And no errors should have been logged

  Scenario: Layout > Suspend simulation is the same setting, reached from the viewer's own menu
    Then "suspendSimulation" property of network diagram viewer should be "false"
    When user picks "Layout > Suspend simulation" from the context menu of network diagram viewer
    Then "suspendSimulation" property of network diagram viewer should be "true"
    When user remembers the "layout signature" reading of network diagram viewer
    And user clicks on the "node \"M\"" area of network diagram viewer
    Then 447 rows should be selected
    And the "layout signature" reading of network diagram viewer should be as remembered
    When user clears the row selection
    And user picks "Layout > Suspend simulation" from the context menu of network diagram viewer
    Then "suspendSimulation" property of network diagram viewer should be "false"
    And no errors should have been logged

  Scenario: Reset View leaves the graph and its counts where they were
    When user remembers the "layout signature" reading of network diagram viewer
    And user picks "Reset View" from the context menu of network diagram viewer
    Then the "nodes" reading of network diagram viewer should be 4
    And the "edges" reading of network diagram viewer should be 4
    And the "rows of node \"F\"" reading of network diagram viewer should be 553
    And the "layout signature" reading of network diagram viewer should be as remembered
    And network diagram viewer should have a "node \"F\"" area
    And no errors should have been logged
