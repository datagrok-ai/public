@journey @viewers @realizes:viewers.density-plot
Feature: Density plot selection, zoom and axis bounds
  Clicking a bin selects the rows binned into it, the wheel and an Alt-drag move the viewport,
  Reset View puts it back, and the four bound properties pin it by hand.
  The old spec clicked the geometric centre of the canvas and hoped a bin was there — on
  demog-1000 over 50 × 50 hexagons the average bin holds 0.35 rows, so that gesture is a coin
  flip. The plot now reports `densest bin`, and the selection claim is exact: the rows selected
  after the click equal the rows the plot counted in that bin. Hexagons make it exact rather than
  approximate, because `_getHexagonDensities` and the click hit test share one mapping.
  Zoom and Reset View are read as the viewport, not as a direction of ink: the old spec asserted
  only that a reset moved the pixels *less far* from the baseline than the zoom had, which passes
  for any partial restoration.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a density plot viewer with:
      | xColumnName | AGE    |
      | yColumnName | WEIGHT |
    Then the "rows shown" reading of density plot viewer should be 1000
    And the "bin shape" reading of density plot viewer should be "hexagon"
    And density plot viewer should have a "densest bin" area

  Scenario: Clicking the densest bin selects exactly the rows counted in it
    Given user clears the row selection
    Then the "rows in densest bin" reading of density plot viewer should be at least 1
    When user clicks on the "densest bin" area of density plot viewer
    Then the "rows selected" and "rows in densest bin" readings of density plot viewer should be the same
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged

  Scenario: The densest bin's tooltip reports how many rows fell into it
    When user hovers over the "densest bin" area of density plot viewer
    Then tooltip should contain text "rows"
    And no errors should have been logged

  Scenario: An Alt-drag zooms into the box and Reset View restores the viewport
    When user picks "Reset View" from the context menu of density plot viewer
    And user remembers the value range of density plot viewer
    And user drags a zoom box over the "view" area of density plot viewer
    Then density plot viewer should show a narrower value range than before
    And the "x axis span" reading of density plot viewer should be lower than before
    And density plot viewer should have repainted
    When user picks "Reset View" from the context menu of density plot viewer
    Then density plot viewer should show the remembered value range
    And no errors should have been logged

  Scenario: The wheel zooms and Reset View undoes it
    When user picks "Reset View" from the context menu of density plot viewer
    And user remembers the value range of density plot viewer
    And user scrolls the mouse wheel up over the "view" area of density plot viewer
    Then density plot viewer should show a narrower value range than before
    When user picks "Reset View" from the context menu of density plot viewer
    Then density plot viewer should show the remembered value range
    And no errors should have been logged

  Scenario: The four bound properties pin the viewport, and clearing them releases it
    When user sets properties of density plot viewer:
      | xMin | 30 |
      | xMax | 60 |
    Then the "x axis min" reading of density plot viewer should be 30
    And the "x axis max" reading of density plot viewer should be 60
    And the "x axis span" reading of density plot viewer should be 30
    When user sets properties of density plot viewer:
      | yMin | 60  |
      | yMax | 100 |
    Then the "y axis min" reading of density plot viewer should be 60
    And the "y axis max" reading of density plot viewer should be 100
    And density plot viewer should have repainted
    When user sets properties of density plot viewer:
      | xMin |  |
      | xMax |  |
      | yMin |  |
      | yMax |  |
    Then the "x axis span" reading of density plot viewer should be higher than before
    And no errors should have been logged

  Scenario: Bin To Range re-bins over the viewport rather than the column range
    When user drags a zoom box over the "view" area of density plot viewer
    And user remembers the "bins drawn" reading of density plot viewer
    And user sets "binToRange" property of density plot viewer to "true"
    Then the "bins drawn" reading of density plot viewer should not be as remembered
    And density plot viewer should have repainted
    And the "bin to range" reading of density plot viewer should be "true"
    When user sets "binToRange" property of density plot viewer to "false"
    Then the "bin to range" reading of density plot viewer should be "false"
    And density plot viewer should have repainted
    And no errors should have been logged

  Scenario: A logarithmic axis over a column with non-positive values shows the axis warning
    Given user adds a calculated column "BELOW_ZERO" with formula "${WEIGHT} - 200"
    Then the table should have a column "BELOW_ZERO"
    When user sets "yColumnName" property of density plot viewer to "BELOW_ZERO"
    Then density plot viewer should not have a "y warning" area
    When user sets "yAxisType" property of density plot viewer to "logarithmic"
    Then density plot viewer should have a "y warning" area
    And no errors should have been logged
    When user sets "yAxisType" property of density plot viewer to "linear"
    Then density plot viewer should not have a "y warning" area
    When user sets "yColumnName" property of density plot viewer to "WEIGHT"
    And user removes "BELOW_ZERO" column
    Then no errors should have been logged

  Scenario: Inverting the Y axis redraws without moving the value range
    When user remembers the value range of density plot viewer
    And user sets "invertYAxis" property of density plot viewer to "true"
    Then density plot viewer should have repainted
    And density plot viewer should show the remembered value range
    When user sets "invertYAxis" property of density plot viewer to "false"
    Then no errors should have been logged
