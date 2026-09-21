@journey @viewers @realizes:viewers.density-plot
Feature: Density plot binning, bin shape and the colour scale
  How many bins the plot holds, how many it actually painted, what the densest one holds, and how
  the colour scale maps those counts.
  The old spec proved every one of these with an ink threshold — `|Δ| > 70000` for a bin count,
  `> 100000` for another, `> 200000` for Invert Color Scheme — which says a lot of pixels moved and
  nothing about what they now mean. `bins` is what the viewer holds, `bins drawn` is what the last
  frame painted, and `max bin count` is the top of the colour scale; a bin count change moves all
  three in a direction the feature can state.
  The X and Y selectors offer numerical columns only, which is what GROK-17118 is about: typing a
  categorical column's name into the selector must leave the column where it was.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a density plot viewer with:
      | xColumnName | AGE    |
      | yColumnName | WEIGHT |
    Then the "rows shown" reading of density plot viewer should be 1000
    And the "bins" reading of density plot viewer should be 50
    And the "bin shape" reading of density plot viewer should be "hexagon"
    And the "color transform" reading of density plot viewer should be "linear"
    And the "color scheme inverted" reading of density plot viewer should be "false"

  Scenario: Fewer bins means fewer, fuller bins; more bins means more, emptier ones
    Then the "bins drawn" reading of density plot viewer should be at least 1
    When user remembers the "bins drawn" reading of density plot viewer
    And user sets "bins" property of density plot viewer to "5"
    Then the "bins" reading of density plot viewer should be 5
    And the "bins drawn" reading of density plot viewer should be lower than before
    And the "max bin count" reading of density plot viewer should be higher than before
    And density plot viewer should have repainted
    When user sets "bins" property of density plot viewer to "200"
    Then the "bins drawn" reading of density plot viewer should be higher than before
    And the "max bin count" reading of density plot viewer should be lower than before
    When user sets "bins" property of density plot viewer to "50"
    Then the "bins drawn" reading of density plot viewer should be as remembered
    And no errors should have been logged

  Scenario: A single rectangular bin holds every row that was binned
    When user sets properties of density plot viewer:
      | binShape | rectangle |
      | bins     | 1         |
    Then the "bins drawn" reading of density plot viewer should be 1
    And the "rows in densest bin" and "rows shown" readings of density plot viewer should be the same
    When user sets properties of density plot viewer:
      | bins     | 50      |
      | binShape | hexagon |
    Then no errors should have been logged

  Scenario: The colour scale spans zero to the fullest bin, and the transform is reported
    Then the "color scale min" reading of density plot viewer should be 0
    And the "color scale max" and "max bin count" readings of density plot viewer should be the same
    When user sets "colorTransformType" property of density plot viewer to "logarithmic"
    Then the "color transform" reading of density plot viewer should be "logarithmic"
    And the "color scale min" reading of density plot viewer should be 1
    And density plot viewer should have repainted
    When user sets "colorTransformType" property of density plot viewer to "linear"
    Then the "color scale min" reading of density plot viewer should be 0
    And no errors should have been logged

  Scenario: Inverting the colour scheme repaints the bins and the scale
    Then the "color scheme inverted" reading of density plot viewer should be "false"
    When user sets "invertColorScheme" property of density plot viewer to "true"
    Then the "color scheme inverted" reading of density plot viewer should be "true"
    And the "color scale" area of density plot viewer should have repainted
    And density plot viewer should have repainted
    And the "rows in densest bin" reading of density plot viewer should be the same as before
    When user sets "invertColorScheme" property of density plot viewer to "false"
    Then no errors should have been logged

  Scenario: Rectangles and hexagons bin the same rows into different shapes
    Then the "bin shape" reading of density plot viewer should be "hexagon"
    When user remembers the "rows shown" reading of density plot viewer
    And user sets "binShape" property of density plot viewer to "rectangle"
    Then the "bin shape" reading of density plot viewer should be "rectangle"
    And the "rows shown" reading of density plot viewer should be as remembered
    And density plot viewer should have repainted
    And density plot viewer should be painted
    When user sets "binShape" property of density plot viewer to "hexagon"
    Then the "bin shape" reading of density plot viewer should be "hexagon"
    And no errors should have been logged

  Scenario: The columns can be re-picked on the plot itself
    When user picks "HEIGHT" in the "y" column selector of density plot viewer
    Then the "y column" reading of density plot viewer should be "HEIGHT"
    And the "rows shown" reading of density plot viewer should be 872
    When user picks "WEIGHT" in the "y" column selector of density plot viewer
    Then the "y column" reading of density plot viewer should be "WEIGHT"
    And the "rows shown" reading of density plot viewer should be 1000
    And no errors should have been logged

  Scenario: The column selector offers numerical columns only (GROK-17118)
    Then the "x column" reading of density plot viewer should be "AGE"
    When user types "SEX" into the "x" column selector of density plot viewer
    Then the "x column" reading of density plot viewer should be "AGE"
    And the "rows shown" reading of density plot viewer should be 1000
    And no errors should have been logged

  Scenario: A calculated column can be binned like any other
    Given user adds a calculated column "AGE_PLUS" with formula "${AGE} + 1"
    When user sets "xColumnName" property of density plot viewer to "AGE_PLUS"
    Then the "x column" reading of density plot viewer should be "AGE_PLUS"
    And the "rows shown" reading of density plot viewer should be 1000
    And density plot viewer should be painted
    When user sets "xColumnName" property of density plot viewer to "AGE"
    And user removes "AGE_PLUS" column
    Then no errors should have been logged
