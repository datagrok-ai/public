@viewers @realizes:viewers.trellis-plot
Feature: Trellis plot Pick Up / Apply, inner color coding and category scrolling
  Pick Up on one trellis and Apply on another carries its Y split, its inner viewer, its legend and
  its title, and a later change to the first leaves the second as it was. A pie chart or a box plot
  inside, colored by RACE, repaints every cell (the pie chart has no Marker Color, which the case
  names: its slices are colored by its Category; the box plot's is Marker Color). The X column selector's own menu (Reset X columns)
  clears the X split while the axis is paged, leaving the four RACE cells; the category scroll slider dragged along its track, and the mouse wheel
  over the grid, bring other categories into the window. Translated from TestTrack
  Viewers/TrellisPlot/trellis-plot.md "Pick Up / Apply" steps 3-7 and "Scrolling" steps 3 and 5,
  trellis-plot-ui.md "Inner viewer color coding", trellis-plot-scroll-categories.md section 1 step 8
  and trellis-plot-split-and-pick-inner.md section 2 step 6. Each scenario starts on a fresh
  demog-1000 view.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  Scenario: Pick Up / Apply carries the Y split, the inner viewer, the legend and the title
    Given user adds a trellis plot viewer
    And user adds a trellis plot viewer
    Then the open tableview should have 2 trellis plot viewers
    When user sets properties of first trellis plot viewer:
      | Y Column Names    | DIS_POP   |
      | Viewer Type       | Bar chart |
      | Legend Visibility | Always    |
      | Legend Position   | Left      |
      | Show Title        | true      |
      | Title             | First     |
    Then "Y Column Names" property of second trellis plot viewer should not be "DIS_POP"
    And second trellis plot viewer should not have a "y label RA" area
    When user picks "Pick Up / Apply > Pick Up" from the viewer menu of first trellis plot viewer
    And user picks "Pick Up / Apply > Apply" from the viewer menu of second trellis plot viewer
    Then properties of second trellis plot viewer should be:
      | Y Column Names  | DIS_POP   |
      | Viewer Type     | Bar chart |
      | Legend Position | Left      |
      | Title           | First     |
    And the "inner viewer type" reading of second trellis plot viewer should be "Bar chart"
    And second trellis plot viewer should have a "y label RA" area
    And title of second trellis plot viewer should have text "First"
    When user sets "Y Column Names" property of first trellis plot viewer to "RACE"
    Then the "y categories" reading of first trellis plot viewer should be 4
    And the "y categories" reading of second trellis plot viewer should be 6
    And "Y Column Names" property of second trellis plot viewer should be "DIS_POP"
    And no errors should have been logged

  Scenario Outline: A <type> inside, colored by RACE, repaints every cell
    Given user adds a trellis plot viewer with:
      | X Column Names | SEX    |
      | Y Column Names | RACE   |
      | Viewer Type    | <type> |
    Then the "cells drawn" reading of trellis plot viewer should be 8
    When user remembers the "cell signature F | Caucasian" reading of trellis plot viewer
    And user remembers the "cell signature M | Asian" reading of trellis plot viewer
    And user remembers the "cell signature F | Black" reading of trellis plot viewer
    And user remembers the "cell signature M | Other" reading of trellis plot viewer
    And user sets "<color property>" inner property of trellis plot viewer to "RACE"
    Then "<color property>" inner property of trellis plot viewer should be "RACE"
    And the "cell signature F | Caucasian" reading of trellis plot viewer should not be as remembered
    And the "cell signature M | Asian" reading of trellis plot viewer should not be as remembered
    And the "cell signature F | Black" reading of trellis plot viewer should not be as remembered
    And the "cell signature M | Other" reading of trellis plot viewer should not be as remembered
    And the "cells drawn" reading of trellis plot viewer should be 8
    And no errors should have been logged

    Examples:
      | type      | color property        |
      | Pie chart | categoryColumnName    |
      | Box plot  | markerColorColumnName |

  Scenario: The X selector's menu resets a paged X split
    Given user adds a trellis plot viewer with:
      | X Column Names | SEX, DIS_POP |
      | Y Column Names | RACE         |
    Then the cells of trellis plot viewer should be 5 wide and 4 tall
    When user clicks on the "x plus" area of trellis plot viewer
    Then the cells of trellis plot viewer should be 6 wide and 4 tall
    When user right-clicks on the "x selector 1" area of trellis plot viewer
    Then the open menu should list "Reset X columns"
    When user picks "Reset X columns" from the open menu
    Then "X Column Names" property of trellis plot viewer should be ""
    And the "cells drawn" reading of trellis plot viewer should be 4
    And trellis plot viewer should report no error
    And no errors should have been logged

  Scenario: Dragging the category scroll handle and the wheel bring other categories in
    Given user adds a trellis plot viewer with:
      | X Column Names | SEX, DIS_POP  |
      | Y Column Names | DIS_POP, RACE |
    Then the cells of trellis plot viewer should be 5 wide and 5 tall
    And trellis plot viewer should have an "x label F" area
    And trellis plot viewer should not have an "x label M" area
    When user drags the "x scroll handle" area of trellis plot viewer by 200 pixels to the right
    Then trellis plot viewer should have an "x label M" area
    And the cells of trellis plot viewer should be 5 wide and 5 tall
    And trellis plot viewer should not have a "y label PsA" area
    When user scrolls the mouse wheel down 3 times over the "view" area of trellis plot viewer
    Then trellis plot viewer should have a "y label PsA" area
    And no errors should have been logged
