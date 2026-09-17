@journey @viewers @realizes:viewers.trellis-plot
Feature: Trellis plot global scale, axes and range sliders
  Global Scale puts every cell on one pair of axes: the trellis then draws the axis strips its inner
  viewer only asks for under global scale, owns a range slider per cell plus the shared one, and
  every flip of the setting redraws every cell. The axes obey Always, Never and Auto; Show Range
  Sliders gates the sliders without taking the strip away; "Reset Inner Range Sliders" is offered
  only while a slider exists, and it puts the cells back exactly where the shared slider took them
  from. Last, the wheel over a cell zooms nothing until the inner viewer's own Allow Zoom says so.
  One journey on demog-1000 with SEX by RACE and a scatter plot inside.

  Not translated here: "nothing repaints while the viewer is idle" — the trellis paints no canvas
  of its own, so idleness has no pixel to be read from; the cell signatures below carry the same
  claim wherever a change is expected to leave a cell alone.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names | SEX          |
      | Y Column Names | RACE         |
      | Viewer Type    | Scatter plot |
    Then the "cells" reading of trellis plot viewer should be 8
    And the "x axis sliders" reading of trellis plot viewer should be 0

  Scenario: Global Scale redraws every cell on every flip
    Then trellis plot viewer should not have an "x axis" area
    When user sets "Global Scale" property of trellis plot viewer to "true"
    Then the "cell signature F | Caucasian" reading of trellis plot viewer should differ from before
    And the "cell signature M | Asian" reading of trellis plot viewer should differ from before
    And trellis plot viewer should have an "x axis" area
    And trellis plot viewer should have a "y axis" area
    When user sets "Global Scale" property of trellis plot viewer to "false"
    Then the "cell signature F | Caucasian" reading of trellis plot viewer should differ from before
    And the "cell signature M | Asian" reading of trellis plot viewer should differ from before
    And trellis plot viewer should not have an "x axis" area
    When user sets "Global Scale" property of trellis plot viewer to "true"
    Then the "cell signature F | Caucasian" reading of trellis plot viewer should differ from before
    And the "cell signature M | Asian" reading of trellis plot viewer should differ from before
    And no errors should have been logged

  Scenario: Show X and Y Axes obey Always, Never and Auto
    When user sets properties of trellis plot viewer:
      | Show X Axes | Always |
      | Show Y Axes | Always |
    Then trellis plot viewer should have an "x axis" area
    And trellis plot viewer should have an "x axis cell 1" area
    And trellis plot viewer should have a "y axis" area
    And the "x axis sliders" reading of trellis plot viewer should be 3
    And the "y axis sliders" reading of trellis plot viewer should be 5
    When user sets properties of trellis plot viewer:
      | Show X Axes | Never |
      | Show Y Axes | Never |
    Then trellis plot viewer should not have an "x axis" area
    And trellis plot viewer should not have a "y axis" area
    And the "x axis sliders" reading of trellis plot viewer should be 0
    And the "y axis sliders" reading of trellis plot viewer should be 0
    When user sets properties of trellis plot viewer:
      | Show X Axes | Auto |
      | Show Y Axes | Auto |
    Then trellis plot viewer should have an "x axis" area
    And the "x axis sliders" reading of trellis plot viewer should be 3
    And the "y axis sliders" reading of trellis plot viewer should be 5
    And no errors should have been logged

  Scenario: Show Range Sliders gates the sliders and keeps the strip
    When user sets "Show Range Sliders" property of trellis plot viewer to "false"
    Then the "x axis sliders" reading of trellis plot viewer should be 0
    And the "y axis sliders" reading of trellis plot viewer should be 0
    And trellis plot viewer should have an "x axis" area
    When user sets "Show Range Sliders" property of trellis plot viewer to "true"
    Then the "x axis sliders" reading of trellis plot viewer should be 3
    And the "y axis sliders" reading of trellis plot viewer should be 5
    And no errors should have been logged

  Scenario: Reset Inner Range Sliders is offered only while a slider exists
    When user right-clicks on the "view" area of trellis plot viewer
    Then "Properties..." menu item in context menu should be visible
    And "Reset Inner Range Sliders" menu item in context menu should be visible
    When user closes the context menu
    And user sets properties of trellis plot viewer:
      | Show X Axes | Never |
      | Show Y Axes | Never |
    And user right-clicks on the "view" area of trellis plot viewer
    Then "Properties..." menu item in context menu should be visible
    And "Reset Inner Range Sliders" menu item in context menu should be absent
    When user closes the context menu
    And user sets properties of trellis plot viewer:
      | Show X Axes | Always |
      | Show Y Axes | Always |
    And user right-clicks on the "view" area of trellis plot viewer
    Then "Reset Inner Range Sliders" menu item in context menu should be visible
    When user closes the context menu
    Then no errors should have been logged

  Scenario: The shared slider re-bounds every cell and Reset puts them back exactly
    When user hovers over the "cell body F | Caucasian" area of trellis plot viewer
    Then trellis plot viewer should have an "x range slider" area
    And trellis plot viewer should have an "x range slider max handle" area
    When user remembers the "cell signature F | Caucasian" reading of trellis plot viewer
    And user remembers the "cell signature M | Asian" reading of trellis plot viewer
    And user drags the "x range slider max handle" area of trellis plot viewer by 60 pixels to the left
    Then the "cell signature F | Caucasian" reading of trellis plot viewer should not be as remembered
    And the "cell signature M | Asian" reading of trellis plot viewer should not be as remembered
    When user picks "Reset Inner Range Sliders" from the context menu of the "view" area of trellis plot viewer
    Then the "cell signature F | Caucasian" reading of trellis plot viewer should be as remembered
    And the "cell signature M | Asian" reading of trellis plot viewer should be as remembered
    And no errors should have been logged

  Scenario: The wheel zooms a cell only once Allow Zoom says so
    When user sets properties of trellis plot viewer:
      | Global Scale | false |
    And user remembers the "cell signature F | Caucasian" reading of trellis plot viewer
    And user scrolls the mouse wheel down over the "cell body F | Caucasian" area of trellis plot viewer
    Then the "cell signature F | Caucasian" reading of trellis plot viewer should be as remembered
    When user sets "allowZoom" inner property of trellis plot viewer to "true"
    And user scrolls the mouse wheel down over the "cell body F | Caucasian" area of trellis plot viewer
    Then the "cell signature F | Caucasian" reading of trellis plot viewer should not be as remembered
    When user sets "allowZoom" inner property of trellis plot viewer to "false"
    And user remembers the "cell signature F | Caucasian" reading of trellis plot viewer
    And user scrolls the mouse wheel down over the "cell body F | Caucasian" area of trellis plot viewer
    Then the "cell signature F | Caucasian" reading of trellis plot viewer should be as remembered
    And no errors should have been logged

  Scenario: The wheel leaves a bar chart cell alone too
    When user sets "Viewer Type" property of trellis plot viewer to "Bar chart"
    Then the "inner viewer type" reading of trellis plot viewer should be "Bar chart"
    And the "distinct cell signatures" reading of trellis plot viewer should be at least 2
    When user remembers the "cell signature F | Caucasian" reading of trellis plot viewer
    And user scrolls the mouse wheel down over the "cell body F | Caucasian" area of trellis plot viewer
    Then the "cell signature F | Caucasian" reading of trellis plot viewer should be as remembered
    When user sets "Viewer Type" property of trellis plot viewer to "Scatter plot"
    Then the "inner viewer type" reading of trellis plot viewer should be "Scatter plot"
    And no errors should have been logged
