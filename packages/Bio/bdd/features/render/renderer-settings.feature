@realizes:bio.panel.sequence-renderer @realizes:bio.rendering.biln
Feature: Sequence renderer settings on the context panel
  A Macromolecule column's header shows, in the context panel's Bioinformatics group, the Sequence
  Renderer pane, whose Renderer Settings set how the column's cells are drawn. A BILN column gets
  Show Connections there, which starts from the column's .mm.cellRenderer.showConnections tag when
  the column has one, else from the package's ShowConnections setting (on by default). The biln
  column of the BILN sample has no such tag, so the box starts checked; unchecking it writes "false"
  into the tag and checking it again writes "true".

  What the cells then draw, a ring bond as a line between its two monomers instead of the C(1,3)
  labels once the rows are tall enough, is the renderer's and is tested by Bio's package tests
  (src/tests/renderers-connections-tests.ts, category `renderers: connections`); the feature claims the checkbox and the tag it writes (the bdd library's CLAUDE.md,
  "What never becomes a feature"). Not claimed: a tag the column already carries as the starting
  value, and the package setting turned off (a package setting the feature would have to put back).

  Background:
    Given user is logged in
    And the Bio package is initialized
    And user opens BILN dataset
    And the context panel is open
    Then "biln" column should have units "biln"

  Scenario: Show Connections starts checked on a BILN column and toggling it writes the column's tag
    When user clicks on the "header biln" area of grid
    Then the context panel should show "biln"
    When user expands Bioinformatics accordion header in context panel
    And user expands "Sequence Renderer" accordion header in context panel
    And user expands "Renderer Settings" accordion header in context panel
    Then "Show Connections" checkbox in context panel should be checked
    When user unchecks "Show Connections" checkbox in context panel
    Then "biln" column should have tag ".mm.cellRenderer.showConnections" equal to "false"
    When user checks "Show Connections" checkbox in context panel
    Then "biln" column should have tag ".mm.cellRenderer.showConnections" equal to "true"
    And no error or warning balloon should have been shown
    And no errors should have been logged
