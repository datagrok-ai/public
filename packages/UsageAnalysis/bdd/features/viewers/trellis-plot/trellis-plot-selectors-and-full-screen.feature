@journey @viewers @realizes:viewers.trellis-plot
Feature: Trellis plot selectors, control panel and the full-screen cell
  Show X Selectors, Show Y Selectors and Show Control Panel each take their strip off the screen and
  bring it back. With Auto Layout on, an X selector strip switched off stays off through a shrink
  that hides every control and the restore that brings the control panel back. The full-screen icon a
  hovered cell carries opens that cell's viewer on its own in a full-screen dialog named after the
  cell's categories, and closing it returns to the trellis.
  The X strip is read as absent after the restore because the restore forces a relayout (the control
  panel and the Y strip coming back are its witnesses).
  Translated from the "Selectors" and "Allow viewer full screen" sections of TestTrack
  Viewers/TrellisPlot/trellis-plot.md (steps 4-5 of the latter; the icon's hover and the Off setting
  are in trellis-plot-tiles-and-layout). One journey on demog-1000 with SEX by RACE and a scatter
  plot inside; every scenario puts back what it changed.

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

  Scenario: An X selector strip switched off stays off through Auto Layout's shrink and restore
    Then "Auto Layout" property of trellis plot viewer should be "true"
    When user sets "Show X Selectors" property of trellis plot viewer to "false"
    And user resizes trellis plot viewer to 240 by 240
    Then trellis plot viewer should not have a "control panel" area
    And trellis plot viewer should not have a "y selectors" area
    When user restores the size of trellis plot viewer
    Then trellis plot viewer should have a "control panel" area
    And trellis plot viewer should have a "y selectors" area
    And trellis plot viewer should not have a "x selectors" area
    When user sets "Show X Selectors" property of trellis plot viewer to "true"
    Then trellis plot viewer should have a "x selectors" area
    And no errors should have been logged

  Scenario: The full-screen icon of a cell opens its viewer on its own
    When user hovers over the "cell body F | Caucasian" area of trellis plot viewer
    Then trellis plot viewer should have a "full screen icon" area
    When user clicks on the "full screen icon" area of trellis plot viewer
    Then "SEX: F, RACE: Caucasian" dialog should be visible
    And the canvases of "SEX: F, RACE: Caucasian" dialog should be painted in at least 2 colors
    When user presses Escape
    Then "SEX: F, RACE: Caucasian" dialog should be absent
    And the "cells" reading of trellis plot viewer should be 8
    And no errors should have been logged
