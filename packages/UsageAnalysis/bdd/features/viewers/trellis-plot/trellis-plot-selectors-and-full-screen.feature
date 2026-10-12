@journey @viewers @realizes:viewers.trellis-plot
Feature: Trellis plot selectors, control panel and the full-screen cell
  Show Control Panel takes the viewer selector strip off the screen and brings it back; Show X Selectors
  and Show Y Selectors switch their strips off and on, read as the areas the trellis reports, which
  follow the layout's own flag rather than what is visible (MISSING.md). The full-screen icon a hovered
  cell carries opens a full-screen dialog named after the cell's categories, with a painted canvas,
  and Escape closes it. Translated from the "Selectors" (steps 1-4) and "Allow viewer full screen"
  (steps 4-5) sections of TestTrack Viewers/TrellisPlot/trellis-plot.md; the icon's hover and the Off
  setting are in trellis-plot-tiles-and-layout. Step 5 of "Selectors" (an X strip switched off stays
  off through Auto Layout's shrink and restore) is not translated: with the flag behind the area, the
  claim could not fail (MISSING.md). One journey on demog-1000 with SEX by RACE and a scatter plot
  inside; each scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names | SEX          |
      | Y Column Names | RACE         |
      | Viewer Type    | Scatter plot |
    Then the "cells" reading of trellis plot viewer should be 8
    And trellis plot viewer should have a "x selectors" area
    And trellis plot viewer should have a "y selectors" area
    And trellis plot viewer should have a "control panel" area

  Scenario: Each selector strip and the control panel leave the screen and come back
    When user sets "Show X Selectors" property of trellis plot viewer to "false"
    Then trellis plot viewer should not have a "x selectors" area
    And trellis plot viewer should have a "y selectors" area
    When user sets "Show Y Selectors" property of trellis plot viewer to "false"
    Then trellis plot viewer should not have a "y selectors" area
    When user sets "Show Control Panel" property of trellis plot viewer to "false"
    Then trellis plot viewer should not have a "control panel" area
    When user sets properties of trellis plot viewer:
      | Show X Selectors   | true |
      | Show Y Selectors   | true |
      | Show Control Panel | true |
    Then trellis plot viewer should have a "x selectors" area
    And trellis plot viewer should have a "y selectors" area
    And trellis plot viewer should have a "control panel" area
    And the "cells" reading of trellis plot viewer should be 8
    And no errors should have been logged

  Scenario: The full-screen icon of a cell opens a dialog named after the cell
    When user hovers over the "cell body F | Caucasian" area of trellis plot viewer
    Then trellis plot viewer should have a "full screen icon" area
    When user clicks on the "full screen icon" area of trellis plot viewer
    Then "SEX: F, RACE: Caucasian" dialog should be visible
    And the canvases of "SEX: F, RACE: Caucasian" dialog should be painted in at least 2 colors
    When user presses Escape
    Then "SEX: F, RACE: Caucasian" dialog should be absent
    And the "cells" reading of trellis plot viewer should be 8
    And no errors should have been logged
