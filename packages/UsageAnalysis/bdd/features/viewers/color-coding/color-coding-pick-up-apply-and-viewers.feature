@journey @viewers @realizes:viewers.grid
Feature: Pick Up / Apply Coloring from the header menu, and a column's coloring kept by the viewers around it
  On demog-1000. The header menu's Apply Coloring is disabled until a coloring has been picked up;
  Pick Up Coloring on RACE, coloured with colours of its own, and Apply Coloring on a copy of it give
  the copy those colours, not the default palette. A heat map added to the view leaves the grid's colour coding as it was (public#3043).
  Categorical colours set on the grid column are the colours a scatter plot's legend paints its
  categories in (public#1929). Translated from the closed tickets public#3043, public#1929,
  GROK-17058 and the help's "Pick Up / Apply Coloring"; `color-coding.feature` copies a coloring by
  its tags, this feature does it through the menu a person uses.

  Not translated: the legend's own picker, its Cancel and a second pick (github-3132) are
  `legend/legend-across-viewers.feature`; a conditional coloring's legend is
  `guides/scatter-plot-conditional-colors.feature`.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a calculated column "Race_copy" with formula "${RACE}"
    Then "RACE" column should have no color coding
    And "Race_copy" column should have no color coding

  Scenario: Apply Coloring is disabled until a coloring has been picked up
    When user right-clicks on the "header Race_copy" area of grid
    And user hovers over "Color Coding" menu item
    Then "Color-Coding---Apply-Coloring" menu item should be disabled
    And "Color-Coding---Pick-Up-Coloring" menu item should be enabled
    When user closes the context menu
    Then no errors should have been logged

  Scenario: Pick Up Coloring on one column and Apply Coloring on another give it the same colours
    When user colors "RACE" column categorically:
      | Caucasian | #3366CC |
      | Other     | #CC6699 |
    Then "RACE" column should be color-coded categorically
    When user picks "Color Coding > Pick Up Coloring" from the context menu of the "header RACE" area of grid
    And user right-clicks on the "header Race_copy" area of grid
    And user hovers over "Color Coding" menu item
    Then "Color-Coding---Apply-Coloring" menu item should be enabled
    When user closes the context menu
    And user picks "Color Coding > Apply Coloring" from the context menu of the "header Race_copy" area of grid
    Then "Race_copy" column should be color-coded categorically
    And the categorical color of "Caucasian" in "Race_copy" column should be "3366CC"
    And the categorical color of "Other" in "Race_copy" column should be "CC6699"
    And the "color of cell 1 of Race_copy" and "color of cell 1 of RACE" readings of grid should be the same
    And the "color of cell 2 of Race_copy" and "color of cell 2 of RACE" readings of grid should be the same
    And the "color of cell 1 of Race_copy" and "color of cell 2 of Race_copy" readings of grid should differ
    When user picks "Color Coding > Off" from the context menu of the "header Race_copy" area of grid
    And user removes the coloring of "RACE" column
    Then "Race_copy" column should have no color coding
    And "RACE" column should have no color coding
    And no errors should have been logged

  Scenario: A heat map added to the view leaves the grid's colour coding as it was
    When user picks "Color Coding > Conditional" from the context menu of the "header AGE" area of grid
    And user picks "Color Coding > Categorical" from the context menu of the "header SEX" area of grid
    Then "AGE" column should be color-coded conditionally
    And "SEX" column should be color-coded categorically
    When user remembers the "color of cell 1 of AGE" reading of grid
    And user adds a heat map viewer
    Then heat map viewer should be visible
    And "AGE" column should be color-coded conditionally
    And "SEX" column should be color-coded categorically
    And the "color of cell 1 of AGE" reading of grid should be as remembered
    When user clicks on close icon of heat map viewer
    Then heat map viewer should be hidden
    And "AGE" column should be color-coded conditionally
    And "SEX" column should be color-coded categorically
    And the "color of cell 1 of AGE" reading of grid should be as remembered
    When user picks "Color Coding > Off" from the context menu of the "header AGE" area of grid
    And user picks "Color Coding > Off" from the context menu of the "header SEX" area of grid
    Then "AGE" column should have no color coding
    And no errors should have been logged

  Scenario: The categorical colours set on the grid column are the scatter plot legend's
    When user colors "SEX" column categorically:
      | M | #3366CC |
      | F | #CC6699 |
    And user adds a scatter plot viewer with:
      | X     | WEIGHT |
      | Y     | HEIGHT |
      | Color | SEX    |
    Then the "M" item in the legend of scatter plot viewer should be colored "#3366CC"
    And the "F" item in the legend of scatter plot viewer should be colored "#CC6699"
    When user clicks on close icon of scatter plot viewer
    And user removes the coloring of "SEX" column
    Then "SEX" column should have no color coding
    And no errors should have been logged
