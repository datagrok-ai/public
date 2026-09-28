@journey @viewers @realizes:viewers.heat-map
Feature: Heat map colouring
  The two switches that decide how a heat map's cells are filled, read off the pixels of one
  column's band rather than off the whole canvas — so "recoloured" is a claim about the column
  whose scale changed, not about a repaint somewhere on screen.
  Global Color Scaling recolours the band (its pixels change); Heatmap Colors off takes the fill
  away, so the band is left with less ink than it had.

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

  Scenario: Heatmap Colors off stops filling the cells with colour (GROK-20619)
    # Fixed 2026-09-21 in grid_core.dart: the dense heat-map path coloured every uncoded column
    # directly, bypassing the `heatmapColors` gate of the cell renderer; it now shares that gate.
    # Global Color Scaling above provides a positive repaint check on the same band.
    Then the "heatmap colors" reading of heat map viewer should be "true"
    When user sets "heatmapColors" property of heat map viewer to "false"
    Then the "heatmap colors" reading of heat map viewer should be "false"
    And the "column AGE" area of heat map viewer should have less ink than before
    And no errors should have been logged
