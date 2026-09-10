@journey @viewers @realizes:viewers.heat-map
Feature: Heat map colouring
  The two switches that decide how a heat map's cells are filled, read off the pixels of one
  column's band rather than off the whole canvas — so "recoloured" is a claim about the column
  whose scale changed, not about a repaint somewhere on screen.
  Both scenarios make the same claim about the same band, which is what makes the pair worth
  having: Global Color Scaling moves it (measured 48828 pixels differ), Heatmap Colors does not
  move it at all (measured 0).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a heat map viewer
    Then the "is heatmap" reading of heat map viewer should be "true"
    And the "row height" reading of heat map viewer should be between 0 and 8
    And heat map viewer should have a "column AGE" area
    And the "column AGE" area of heat map viewer should be painted

  Scenario: Global Color Scaling recolours every numerical column
    Then the "global color scaling" reading of heat map viewer should be "false"
    When user sets "globalColorScaling" property of heat map viewer to "true"
    Then the "global color scaling" reading of heat map viewer should be "true"
    And the "column AGE" area of heat map viewer should have repainted
    And the "column WEIGHT" area of heat map viewer should have repainted
    And the "rows shown" reading of heat map viewer should be 1000
    When user sets "globalColorScaling" property of heat map viewer to "false"
    Then the "global color scaling" reading of heat map viewer should be "false"
    And the "column AGE" area of heat map viewer should have repainted
    And no errors should have been logged

  @known-failure
  Scenario: Heatmap Colors off stops filling the cells with colour (GROK-20619)
    # `heatmapColors` is declared on the look and the property write lands, but nothing in the
    # render path reads it back for an attached heat map: the AGE column's band comes out
    # pixel-for-pixel identical afterwards — 0 pixels differ, against 48828 for the Global Color
    # Scaling write on the same band in the scenario above. The spec this replaces already carried
    # this as knownOpenBug('GROK-20619'), with a whole-canvas repaint as the claim.
    # Left last: a known failure aborts before its restore step.
    Then the "heatmap colors" reading of heat map viewer should be "true"
    When user sets "heatmapColors" property of heat map viewer to "false"
    Then the "heatmap colors" reading of heat map viewer should be "false"
    And the "column AGE" area of heat map viewer should have repainted
