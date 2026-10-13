@journey @viewers @realizes:viewers.legend
Feature: Legend interaction
  The legend a viewer hosts: what it lists, what a click on a category does to the viewer and to
  the table, how a selection is added to and taken back, and what a click must never do to the
  legend's own placement, and that a selection outlives a row filter but not a change of the
  column it was made on — the item at the same place in the new column's list is not selected and
  the plot draws every row again; the docked legend is measured to its labels, so the same four
  categories with long labels widen it and short ones narrow it again; and a legend placed by Auto
  in the one free corner (CORNERS keeps the points on the bottom row and, past AGE 60, on the top
  row too, so only the top left is empty) leaves it once the regression equation is drawn there
  (GROK-19041). One journey on demog-1000 with a scatter plot colored by RACE; every
  scenario puts back what it changed. The legend publishes its mode, slot and item total itself,
  so nothing here counts rendered rows or reads pixels to find out where it sits. Under Auto the
  scatter plot parks a four-item legend in a free corner, where its list shows only what fits, so
  the Background docks it on the right and the position scenario is the one that moves it. A
  category's picked color is the column's palette and outlives the coding being switched off:
  removing the coding leaves the items (and the markers' default palette) in the picked colors.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X                 | WEIGHT |
      | Y                 | HEIGHT |
      | Color             | RACE   |
      | Legend Visibility | Always |
      | Legend Position   | Right  |
    Then legend of scatter plot viewer should be visible
    And the legend of scatter plot viewer should list 4 items
    And the legend of scatter plot viewer should be docked

  Scenario: The legend lists the categories of the color column
    Then "Caucasian" legend item in legend of scatter plot viewer should be visible
    And "Asian" legend item in legend of scatter plot viewer should be visible
    And "Caucasian" legend item in legend of scatter plot viewer should not be selected
    And no errors should have been logged

  Scenario: A click on a category filters the viewer, not the table
    When user clicks on "Black" item in the legend of scatter plot viewer
    Then "Black" legend item in legend of scatter plot viewer should be visible
    And "Black" legend item in legend of scatter plot viewer should be selected
    And scatter plot viewer should show fewer rows than before
    And scatter plot viewer should have less ink than before
    And all rows should pass the filter
    And no rows should be selected
    When user clicks on "Black" item in the legend of scatter plot viewer
    Then "Black" legend item in legend of scatter plot viewer should be visible
    And "Black" legend item in legend of scatter plot viewer should not be selected
    And scatter plot viewer should show more rows than before
    And scatter plot viewer should have more ink than before

  Scenario: Control-click adds a category and the cross takes one back
    When user clicks on "Black" item in the legend of scatter plot viewer
    And user clicks on "Asian" item in the legend of scatter plot viewer holding Control
    Then "Black" legend item in legend of scatter plot viewer should be selected
    And "Asian" legend item in legend of scatter plot viewer should be selected
    And scatter plot viewer should show more rows than before
    When user clicks on the cross of "Asian" item in the legend of scatter plot viewer
    Then "Asian" legend item in legend of scatter plot viewer should be visible
    And "Asian" legend item in legend of scatter plot viewer should not be selected
    And "Black" legend item in legend of scatter plot viewer should be selected
    And scatter plot viewer should show fewer rows than before
    When user clicks on "Black" item in the legend of scatter plot viewer
    Then scatter plot viewer should show 872 rows

  Scenario: A click never moves the legend, and a filter never moves it either
    When user takes a snapshot of scatter plot viewer
    And user clicks on "Caucasian" item in the legend of scatter plot viewer
    Then the legend of scatter plot viewer should be placed as before
    And the legend of scatter plot viewer should list the same items as before
    When user clicks on "Caucasian" item in the legend of scatter plot viewer
    And user takes a snapshot of scatter plot viewer
    And user filters rows where "SEX" is "F"
    Then the legend of scatter plot viewer should be placed as before
    And the legend of scatter plot viewer should list the same items as before
    When user resets the filter

  Scenario: A selection stays through a row filter and goes with the color column it was made on
    When user clicks on "Black" item in the legend of scatter plot viewer
    Then "Black" legend item in legend of scatter plot viewer should be selected
    And scatter plot viewer should show 25 rows
    When user filters rows where "SEX" is "F"
    Then "Black" legend item in legend of scatter plot viewer should be selected
    And scatter plot viewer should show 17 rows
    When user resets the filter
    Then scatter plot viewer should show 25 rows
    When user sets "Color" property of scatter plot viewer to "DIS_POP"
    Then the legend of scatter plot viewer should list 5 items
    And "Indigestion" legend item in legend of scatter plot viewer should not be selected
    And scatter plot viewer should show 872 rows
    When user sets "Color" property of scatter plot viewer to "RACE"
    Then the legend of scatter plot viewer should list 4 items
    And "Black" legend item in legend of scatter plot viewer should not be selected
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  Scenario: A filter that empties a category drops it from the legend
    When user takes a snapshot of scatter plot viewer
    And user filters rows where "RACE" is one of "Caucasian, Other"
    Then the legend of scatter plot viewer should list fewer items than before
    And "Asian" legend item in legend of scatter plot viewer should be absent
    And "Caucasian" legend item in legend of scatter plot viewer should be visible
    When user resets the filter
    Then the legend of scatter plot viewer should list 4 items

  Scenario: The category colors are the column's, and a recoloring reaches the items
    Then the "Caucasian" and "Asian" items in the legend of scatter plot viewer should be colored differently
    When user colors "RACE" column categorically:
      | Caucasian | #FF0000 |
      | Asian     | #0000FF |
    Then the "Caucasian" item in the legend of scatter plot viewer should be colored "#FF0000"
    And the "Asian" item in the legend of scatter plot viewer should be colored "#0000FF"
    And the categorical color of "Caucasian" in "RACE" column should be "#FF0000"
    And the "view" area of scatter plot viewer should contain the color "#FF0000"
    When user removes the coloring of "RACE" column
    Then "RACE" column should have no color coding
    And the "Caucasian" item in the legend of scatter plot viewer should be colored "#FF0000"

  Scenario: Position and visibility move and hide the legend
    When user sets "Legend Position" property of scatter plot viewer to "Left"
    Then the legend of scatter plot viewer should be on the left
    And the legend of scatter plot viewer should be in the "left" slot
    And the legend of scatter plot viewer should be docked
    When user sets "Legend Position" property of scatter plot viewer to "Bottom"
    Then the legend of scatter plot viewer should be on the bottom
    When user sets "Legend Visibility" property of scatter plot viewer to "Never"
    Then legend of scatter plot viewer should be hidden
    When user sets properties of scatter plot viewer:
      | Legend Visibility | Always |
      | Legend Position   | Right  |
    Then the legend of scatter plot viewer should be on the right
    And the legend of scatter plot viewer should list 4 items
    And no errors should have been logged

  Scenario: The docked legend is as wide as its labels need
    When user adds a calculated column "RACE described" with formula "\"race of the subject in the study: \" + ${RACE}"
    And user sets "Color" property of scatter plot viewer to "RACE described"
    Then the legend of scatter plot viewer should list 4 items
    And the legend of scatter plot viewer should be wider than before
    And the legend of scatter plot viewer should be docked
    When user sets "Color" property of scatter plot viewer to "RACE"
    Then the legend of scatter plot viewer should be narrower than before
    And the legend of scatter plot viewer should list 4 items
    When user removes "RACE described" column
    Then no errors should have been logged

  Scenario: A numerical color column shows a scale instead of items
    When user sets "Color" property of scatter plot viewer to "AGE"
    Then the legend of scatter plot viewer should list 0 items
    And scatter plot viewer should have a "color scale" area
    And scatter plot viewer should have repainted
    When user sets "Color" property of scatter plot viewer to "RACE"
    Then the legend of scatter plot viewer should list 4 items
    And scatter plot viewer should not have a "color scale" area
    And no errors should have been logged

  @GROK-19041
  Scenario: An auto-placed legend leaves the corner the regression equation is drawn in
    When user adds a calculated column "CORNERS" with formula "if(${AGE} > 60, Mod(${AGE}, 2) * 100, 0)"
    And user sets properties of scatter plot viewer:
      | X               | AGE     |
      | Y               | CORNERS |
      | Legend Position | Auto    |
    Then the legend of scatter plot viewer should be in a corner
    And the legend of scatter plot viewer should be in the "leftTop" slot
    When user sets properties of scatter plot viewer:
      | showRegressionLine         | true |
      | showRegressionLineEquation | true |
    Then scatter plot viewer should have a "regression stats" area
    And the legend of scatter plot viewer should not cover the "regression stats" area
    And the legend of scatter plot viewer should list 4 items
    When user sets properties of scatter plot viewer:
      | showRegressionLine | false   |
      | X                  | WEIGHT  |
      | Y                  | HEIGHT  |
      | Legend Position    | Right   |
    And user removes "CORNERS" column
    Then the legend of scatter plot viewer should be docked
    And the legend of scatter plot viewer should be on the right
    And no errors should have been logged
