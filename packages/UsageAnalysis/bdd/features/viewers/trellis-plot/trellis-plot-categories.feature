@journey @viewers @realizes:viewers.trellis-plot
Feature: Trellis plot categories, labels and scrolling
  The categories an axis carries and what the viewport shows of them: a second split column
  multiplying the categories and the viewport clamping at five columns, the label strips and the
  angle they are drawn at, the scroll handle that fills its whole track while everything fits and
  shortens when it does not, the (+)/(-) icons paging one category row in and out and going inert
  at the ends, and packing, which drops the categories a filter leaves empty. One journey on
  demog-1000; every scenario puts back what it changed.

  Not translated here: "Reset X columns" from the axis selector's own context menu, and the drag of
  the category scroll slider and the wheel over the grid — the two gestures the old spec left as
  manual gaps; the handle's extent, asserted here, is the state they would leave.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names  | SEX          |
      | Y Column Names  | RACE         |
      | Viewer Type     | Scatter plot |
      | Pack Categories | false        |
    Then the cells of trellis plot viewer should be 2 wide and 4 tall

  Scenario: A second split column grows the grid to the clamped product
    Then the "x categories" reading of trellis plot viewer should be 2
    When user sets "X Column Names" property of trellis plot viewer to "SEX, DIS_POP"
    Then the "x categories" reading of trellis plot viewer should be 12
    And the cells of trellis plot viewer should be 5 wide and 4 tall
    And the "cells" reading of trellis plot viewer should be 20
    And the "cells drawn" reading of trellis plot viewer should be 18
    When user sets "X Column Names" property of trellis plot viewer to "SEX"
    Then the "x categories" reading of trellis plot viewer should be 2
    And the "cells" reading of trellis plot viewer should be 8
    And no errors should have been logged

  Scenario: Show X and Y Labels remove and restore the label strips
    Then the "x labels shown" reading of trellis plot viewer should be 2
    And the "y labels shown" reading of trellis plot viewer should be 4
    And trellis plot viewer should have an "x label F" area
    And trellis plot viewer should have a "y label Caucasian" area
    When user sets "Show X Labels" property of trellis plot viewer to "false"
    Then the "x labels shown" reading of trellis plot viewer should be 0
    And trellis plot viewer should not have an "x label F" area
    And the "y labels shown" reading of trellis plot viewer should be 4
    When user sets "Show Y Labels" property of trellis plot viewer to "false"
    Then the "y labels shown" reading of trellis plot viewer should be 0
    And trellis plot viewer should not have a "y label Caucasian" area
    When user sets properties of trellis plot viewer:
      | Show X Labels | true |
      | Show Y Labels | true |
    Then the "x labels shown" reading of trellis plot viewer should be 2
    And the "y labels shown" reading of trellis plot viewer should be 4
    And trellis plot viewer should have an "x label F" area
    And no errors should have been logged

  Scenario: Label orientation is horizontal, vertical, or one of each
    When user sets properties of trellis plot viewer:
      | X Labels Orientation | Horz |
      | Y Labels Orientation | Horz |
    Then the "x label angle" reading of trellis plot viewer should be 0
    And the "y label angle" reading of trellis plot viewer should be 0
    When user sets properties of trellis plot viewer:
      | X Labels Orientation | Vert |
      | Y Labels Orientation | Vert |
    Then the "x label angle" reading of trellis plot viewer should be -90
    And the "y label angle" reading of trellis plot viewer should be -90
    When user sets properties of trellis plot viewer:
      | X Labels Orientation | Auto |
      | Y Labels Orientation | Auto |
    Then the "x label angle" reading of trellis plot viewer should be 0
    And the "y label angle" reading of trellis plot viewer should be -90
    And no errors should have been logged

  Scenario: Every category fits, so both scroll handles fill their track
    Then the "x scroll handle share" reading of trellis plot viewer should be 1
    And the "y scroll handle share" reading of trellis plot viewer should be 1
    And no errors should have been logged

  Scenario: An overflowing X axis shortens its handle
    When user sets "X Column Names" property of trellis plot viewer to "SEX, DIS_POP"
    Then the "x categories" reading of trellis plot viewer should be 12
    And the cells of trellis plot viewer should be 5 wide and 4 tall
    And the "x scroll handle share" reading of trellis plot viewer should be lower than before
    And the "y scroll handle share" reading of trellis plot viewer should be 1
    When user sets "X Column Names" property of trellis plot viewer to "SEX"
    Then the "x scroll handle share" reading of trellis plot viewer should be 1
    And no errors should have been logged

  Scenario: An overflowing Y axis does the same
    When user sets "Y Column Names" property of trellis plot viewer to "DIS_POP, RACE"
    Then the "y categories" reading of trellis plot viewer should be 24
    And the cells of trellis plot viewer should be 2 wide and 5 tall
    And the "y scroll handle share" reading of trellis plot viewer should be lower than before
    And the "x scroll handle share" reading of trellis plot viewer should be 1
    When user sets "Y Column Names" property of trellis plot viewer to "RACE"
    Then the "y scroll handle share" reading of trellis plot viewer should be 1
    And the cells of trellis plot viewer should be 2 wide and 4 tall
    And no errors should have been logged

  Scenario: The plus and minus icons page one category row in and out
    When user sets "X Column Names" property of trellis plot viewer to "SEX, DIS_POP"
    Then the cells of trellis plot viewer should be 5 wide and 4 tall
    And x plus icon should be enabled
    When user clicks on the "x plus" area of trellis plot viewer
    Then the cells of trellis plot viewer should be 6 wide and 4 tall
    And the "cells" reading of trellis plot viewer should be 24
    When user clicks on the "x minus" area of trellis plot viewer
    Then the cells of trellis plot viewer should be 5 wide and 4 tall
    And the "cells" reading of trellis plot viewer should be 20
    And no errors should have been logged
    When user sets "X Column Names" property of trellis plot viewer to "SEX"

  Scenario: The icons go inert at the ends
    When user sets "X Column Names" property of trellis plot viewer to "SEX, DIS_POP"
    Then the cells of trellis plot viewer should be 5 wide and 4 tall
    When user clicks on the "x minus" area of trellis plot viewer
    And user clicks on the "x minus" area of trellis plot viewer
    And user clicks on the "x minus" area of trellis plot viewer
    And user clicks on the "x minus" area of trellis plot viewer
    Then the cells of trellis plot viewer should be 1 wide and 4 tall
    And x minus icon should be disabled
    And x plus icon should be enabled
    When user clicks on the "x minus" area of trellis plot viewer
    Then the cells of trellis plot viewer should be 1 wide and 4 tall
    And y plus icon should be disabled
    When user clicks on the "y plus" area of trellis plot viewer
    Then the cells of trellis plot viewer should be 1 wide and 4 tall
    When user sets "X Column Names" property of trellis plot viewer to "SEX"
    Then the cells of trellis plot viewer should be 2 wide and 4 tall
    And no errors should have been logged

  Scenario: Packing drops the categories a filter leaves empty
    When user sets properties of trellis plot viewer:
      | X Column Names  | RACE |
      | Y Column Names  | SEX  |
      | Pack Categories | true |
    Then the "x categories packed" reading of trellis plot viewer should be 4
    And the cells of trellis plot viewer should be 4 wide and 2 tall
    When user adds a categorical filter on "RACE" keeping "Caucasian"
    Then 896 rows should pass the filter
    And the "x categories" reading of trellis plot viewer should be 4
    And the "x categories packed" reading of trellis plot viewer should be 1
    And the cells of trellis plot viewer should be 1 wide and 2 tall
    And trellis plot viewer should show 896 rows
    When user sets "Pack Categories" property of trellis plot viewer to "false"
    Then the "x categories packed" reading of trellis plot viewer should be 4
    And the cells of trellis plot viewer should be 4 wide and 2 tall
    When user sets "Pack Categories" property of trellis plot viewer to "true"
    Then the cells of trellis plot viewer should be 1 wide and 2 tall
    When user hovers over "RACE" filter card
    And user clicks on close of "RACE" filter card
    Then "RACE" filter card should be absent
    And all rows should pass the filter
    And the cells of trellis plot viewer should be 4 wide and 2 tall
    And no errors should have been logged
    When user sets properties of trellis plot viewer:
      | X Column Names  | SEX   |
      | Y Column Names  | RACE  |
      | Pack Categories | false |
