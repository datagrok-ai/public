@journey @viewers @realizes:viewers.trellis-plot
Feature: Trellis plot axes without global scale, the Y slider, paging ends and packing
  The axis settings take effect only under Global Scale: with it off, axes set to Always and range
  sliders on still draw no axis strip and own no slider, and the inner viewer's Allow Zoom starts
  off. Under Global Scale the shared Y range slider re-bounds every cell and Reset Inner Range Sliders
  puts them back, as the X one does. The paging icons go inert at both ends: at entry, when every
  category fits, (+) has nothing to add and (-) pages out; at the far end (+) stops adding. Packing a
  two-column X axis drops the combinations no row holds: RACE by SEVERITY makes 20; three (Critical with Asian,
  Black and Other) hold no row, and Asian Medium holds rows none of which has both the HEIGHT and the
  WEIGHT the scatter plot inside draws, so 16 are left. Translated from TestTrack
  Viewers/TrellisPlot/trellis-plot-global-scale-axes.md (section 1 steps 1-4, section 3 steps 1-2),
  trellis-plot.md "Range sliders with global scale" step 7, trellis-plot-scroll-categories.md
  (section 1 steps 1-6, section 2). One journey on demog-1000 with SEX by RACE and a scatter plot
  inside; every scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names | SEX          |
      | Y Column Names | RACE         |
      | Viewer Type    | Scatter plot |
    Then the "cells" reading of trellis plot viewer should be 8
    And "Global Scale" property of trellis plot viewer should be "false"

  Scenario: Without Global Scale the axis settings draw nothing, and Allow Zoom starts off
    Then "allowZoom" inner property of trellis plot viewer should be "false"
    When user sets properties of trellis plot viewer:
      | Show X Axes        | Always |
      | Show Y Axes        | Always |
      | Show Range Sliders | true   |
    Then trellis plot viewer should not have an "x axis" area
    And trellis plot viewer should not have a "y axis" area
    And the "x axis sliders" reading of trellis plot viewer should be 0
    And the "y axis sliders" reading of trellis plot viewer should be 0
    When user sets "Global Scale" property of trellis plot viewer to "true"
    Then trellis plot viewer should have an "x axis" area
    And trellis plot viewer should have a "y axis" area
    And the "y axis sliders" reading of trellis plot viewer should be 5
    And no errors should have been logged

  Scenario: The shared Y slider re-bounds every cell and Reset puts them back exactly
    When user hovers over the "cell body F | Caucasian" area of trellis plot viewer
    Then trellis plot viewer should have a "y range slider" area
    And trellis plot viewer should have a "y range slider max handle" area
    When user remembers the "cell signature F | Caucasian" reading of trellis plot viewer
    And user remembers the "cell signature M | Asian" reading of trellis plot viewer
    And user drags the "y range slider max handle" area of trellis plot viewer by 40 pixels to the down
    Then the "cell signature F | Caucasian" reading of trellis plot viewer should not be as remembered
    And the "cell signature M | Asian" reading of trellis plot viewer should not be as remembered
    When user picks "Reset Inner Range Sliders" from the context menu of the "view" area of trellis plot viewer
    Then the "cell signature F | Caucasian" reading of trellis plot viewer should be as remembered
    And the "cell signature M | Asian" reading of trellis plot viewer should be as remembered
    When user sets properties of trellis plot viewer:
      | Global Scale       | false |
      | Show X Axes        | Auto  |
      | Show Y Axes        | Auto  |
      | Show Range Sliders | false |
    Then trellis plot viewer should not have an "x axis" area
    And no errors should have been logged

  Scenario: At entry every category fits, so (+) has nothing to add and (-) pages out
    When user sets properties of trellis plot viewer:
      | X Column Names  | DIS_POP |
      | Pack Categories | false   |
    Then the "x categories" reading of trellis plot viewer should be 6
    And the cells of trellis plot viewer should be 6 wide and 4 tall
    And x plus icon should be disabled
    And x minus icon should be enabled
    When user clicks on the "x plus" area of trellis plot viewer
    Then the cells of trellis plot viewer should be 6 wide and 4 tall
    When user clicks on the "x minus" area of trellis plot viewer
    Then the cells of trellis plot viewer should be 5 wide and 4 tall
    And x plus icon should be enabled
    When user clicks on the "x plus" area of trellis plot viewer
    Then the cells of trellis plot viewer should be 6 wide and 4 tall
    And no errors should have been logged

  Scenario: At the far end (+) stops adding and (-) still pages
    When user sets "X Column Names" property of trellis plot viewer to "SEX, DIS_POP"
    Then the "x categories" reading of trellis plot viewer should be 12
    And the cells of trellis plot viewer should be 5 wide and 4 tall
    When user clicks on the "x plus" area of trellis plot viewer
    And user clicks on the "x plus" area of trellis plot viewer
    And user clicks on the "x plus" area of trellis plot viewer
    And user clicks on the "x plus" area of trellis plot viewer
    And user clicks on the "x plus" area of trellis plot viewer
    And user clicks on the "x plus" area of trellis plot viewer
    And user clicks on the "x plus" area of trellis plot viewer
    Then the cells of trellis plot viewer should be 12 wide and 4 tall
    And x plus icon should be disabled
    When user clicks on the "x plus" area of trellis plot viewer
    Then the cells of trellis plot viewer should be 12 wide and 4 tall
    When user clicks on the "x minus" area of trellis plot viewer
    Then the cells of trellis plot viewer should be 11 wide and 4 tall
    And x plus icon should be enabled
    And no errors should have been logged

  Scenario: Packing a two-column X axis drops the combinations no row holds
    When user sets properties of trellis plot viewer:
      | X Column Names  | RACE, SEVERITY |
      | Y Column Names  | SEX            |
      | Pack Categories | true           |
    Then the "x categories" reading of trellis plot viewer should be 20
    And the "x categories packed" reading of trellis plot viewer should be 16
    And the "x scroll handle share" reading of trellis plot viewer should be 0.3125
    When user sets "Pack Categories" property of trellis plot viewer to "false"
    Then the "x categories packed" reading of trellis plot viewer should be 20
    And the "x scroll handle share" reading of trellis plot viewer should be 0.25
    When user sets properties of trellis plot viewer:
      | X Column Names  | SEX  |
      | Y Column Names  | RACE |
      | Pack Categories | true |
    Then the "cells" reading of trellis plot viewer should be 8
    And no errors should have been logged
