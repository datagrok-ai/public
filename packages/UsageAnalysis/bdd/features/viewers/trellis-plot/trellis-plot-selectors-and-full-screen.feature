@journey @viewers @realizes:viewers.trellis-plot
Feature: Trellis plot selectors, control panel and the full-screen cell
  Show Control Panel takes the viewer selector strip off the screen and brings it back; Show X Selectors
  and Show Y Selectors switch their strips off and on, read as the areas the trellis reports for the
  strips and their pickers while they are visible — a strip switched off keeps its box and is only
  hidden. With Auto Layout on, an X strip switched off stays off while the viewer is shrunk until Auto
  Layout hides its controls and then restored. The full-screen icon a hovered cell carries opens a
  full-screen dialog named after the cell's categories, with a painted canvas and the cell's rows, and
  Escape closes it. Translated from the "Selectors" (steps 1-5) and "Allow viewer full screen" (steps
  4-5) sections of TestTrack Viewers/TrellisPlot/trellis-plot.md; the icon's hover and the Off setting
  are in trellis-plot-tiles-and-layout. One journey on demog-1000 with SEX by RACE and a scatter plot
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
    And trellis plot viewer should not have a "x selector 1" area
    And trellis plot viewer should have a "y selectors" area
    And trellis plot viewer should have a "y selector 1" area
    When user sets "Show Y Selectors" property of trellis plot viewer to "false"
    Then trellis plot viewer should not have a "y selectors" area
    And trellis plot viewer should not have a "y selector 1" area
    When user sets "Show Control Panel" property of trellis plot viewer to "false"
    Then trellis plot viewer should not have a "control panel" area
    When user sets properties of trellis plot viewer:
      | Show X Selectors   | true |
      | Show Y Selectors   | true |
      | Show Control Panel | true |
    Then trellis plot viewer should have a "x selectors" area
    And trellis plot viewer should have a "x selector 1" area
    And trellis plot viewer should have a "y selectors" area
    And trellis plot viewer should have a "y selector 1" area
    And trellis plot viewer should have a "control panel" area
    And the "cells" reading of trellis plot viewer should be 8
    And no errors should have been logged

  Scenario: An X strip switched off stays off while Auto Layout shrinks the trellis and restores it
    Then "Auto Layout" property of trellis plot viewer should be "true"
    When user sets "Show X Selectors" property of trellis plot viewer to "false"
    Then trellis plot viewer should not have a "x selector 1" area
    And trellis plot viewer should have a "y selector 1" area
    When user resizes trellis plot viewer to 240 by 240
    Then trellis plot viewer should not have a "control panel" area
    And trellis plot viewer should not have a "y selector 1" area
    When user restores the size of trellis plot viewer
    Then trellis plot viewer should have a "control panel" area
    And trellis plot viewer should have a "y selector 1" area
    And trellis plot viewer should not have a "x selector 1" area
    When user sets "Show X Selectors" property of trellis plot viewer to "true"
    Then trellis plot viewer should have a "x selector 1" area
    And no errors should have been logged

  Scenario: The full-screen icon of a cell opens a dialog named after the cell
    When user hovers over the "cell body F | Caucasian" area of trellis plot viewer
    Then trellis plot viewer should have a "full screen icon" area
    When user clicks on the "full screen icon" area of trellis plot viewer
    Then "SEX: F, RACE: Caucasian" dialog should be visible
    And the canvases of "SEX: F, RACE: Caucasian" dialog should be painted in at least 2 colors
    And the "full screen rows" reading of trellis plot viewer should be 480
    When user presses Escape
    Then "SEX: F, RACE: Caucasian" dialog should be absent
    And trellis plot viewer should not report a "full screen rows" reading
    And the "cells" reading of trellis plot viewer should be 8
    And no errors should have been logged
