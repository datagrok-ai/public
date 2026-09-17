@journey @viewers @realizes:viewers.grid
Feature: Colour coding across columns
  What a colouring keeps when it is switched off, copied to another column, linked to another
  column's, or inverted. The types themselves (Linear, Conditional, Categorical from the header
  menu) and the grid-wide Color Coding are covered by `grid/grid-appearance.feature`; this journey
  is about the state that survives an edit.
  The claim is always the colour the grid renderer resolved for a cell (`color of cell <r> of <c>`),
  never the tag alone — a tag a renderer ignores would otherwise pass.
  One journey on demog-1000 (row 1: AGE 26, SEX F, RACE Caucasian; row 2: AGE 30, SEX F, RACE Other;
  row 3: AGE 58, SEX F; row 4: AGE 45, SEX M). Each scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then grid should show 1000 rows
    And "SEX" column should have no color coding
    And "AGE" column should have no color coding

  Scenario: A colouring switched off and back on keeps the colours it was given
    When user colors "SEX" column categorically:
      | M | #3366CC |
      | F | #CC6699 |
    Then "SEX" column should be color-coded categorically
    And the categorical color of "M" in "SEX" column should be "3366CC"
    And the "color of cell 1 of SEX" and "color of cell 4 of SEX" readings of grid should differ
    When user removes the coloring of "SEX" column
    Then "SEX" column should have no color coding
    And the "color of cell 1 of SEX" and "color of cell 4 of SEX" readings of grid should be the same
    When user colors "SEX" column categorically again
    Then "SEX" column should be color-coded categorically
    And the categorical color of "M" in "SEX" column should be "3366CC"
    And the categorical color of "F" in "SEX" column should be "CC6699"
    And the "color of cell 1 of SEX" and "color of cell 4 of SEX" readings of grid should differ
    When user removes the coloring of "SEX" column
    Then no errors should have been logged

  Scenario: Applying one column's colouring to another gives it the same colours
    When user colors "AGE" column conditionally:
      | <30   | #00FF00 |
      | 30-60 | #FFFF00 |
      | >60   | #FF0000 |
    Then "AGE" column should be color-coded conditionally
    And the "color of cell 1 of AGE" reading of grid should be "#00ff00"
    And the "color of cell 3 of AGE" reading of grid should be "#ffff00"
    When user adds a calculated column "Age_copy" with formula "${AGE}"
    Then "Age_copy" column should have no color coding
    And the "color of cell 1 of AGE" and "color of cell 1 of Age_copy" readings of grid should differ
    When user applies the coloring of "AGE" column to "Age_copy" column
    Then "Age_copy" column should be color-coded conditionally
    And the "color of cell 1 of AGE" and "color of cell 1 of Age_copy" readings of grid should be the same
    And the "color of cell 3 of AGE" and "color of cell 3 of Age_copy" readings of grid should be the same
    When user removes "Age_copy" column
    Then the table should not have a column "Age_copy"
    And no errors should have been logged

  Scenario: A linked column paints its cells with the source column's colours
    When user colors "RACE" column linked to "AGE" column
    Then "RACE" column should be color-coded linked
    And the coloring of "RACE" column should be linked to "AGE" column
    And the "color of cell 1 of RACE" reading of grid should be "#00ff00"
    And the "color of cell 3 of RACE" reading of grid should be "#ffff00"
    And the "color of cell 1 of RACE" and "color of cell 3 of RACE" readings of grid should differ
    And no errors should have been logged

  Scenario: A change of the source's type leaves the link in place
    When user colors "AGE" column linearly from "#0000FF" to "#FF0000"
    Then "AGE" column should be color-coded linearly
    And "RACE" column should be color-coded linked
    And the coloring of "RACE" column should be linked to "AGE" column
    And the "color of cell 1 of RACE" and "color of cell 1 of AGE" readings of grid should be the same
    And the "color of cell 1 of RACE" reading of grid should not be "#00ff00"
    When user colors "AGE" column conditionally:
      | <40 | #00FF00 |
      | >40 | #FF0000 |
    Then "AGE" column should be color-coded conditionally
    And "RACE" column should be color-coded linked
    And the "color of cell 1 of RACE" reading of grid should be "#00ff00"
    And the "color of cell 3 of RACE" reading of grid should be "#ff0000"
    And no errors should have been logged

  Scenario: A linked colouring applied to the text paints the letters, not the cell
    When user colors the text of "HEIGHT" column linked to "AGE" column
    Then "HEIGHT" column should be color-coded linked
    And the text of "HEIGHT" column should be color-coded
    And the "color of cell 1 of HEIGHT" and "color of cell 1 of AGE" readings of grid should be the same
    And the "color of cell 1 of HEIGHT" and "color of cell 3 of HEIGHT" readings of grid should differ
    And the "cell 1 of HEIGHT" area of grid should contain the color "#00FF00"
    When user removes the coloring of "HEIGHT" column
    And user removes the coloring of "RACE" column
    Then "RACE" column should have no color coding
    And no errors should have been logged

  Scenario: A five-level chain of links all reports Linked
    When user colors "SEX" column linked to "AGE" column
    And user colors "DIS_POP" column linked to "SEX" column
    And user colors "CONTROL" column linked to "DIS_POP" column
    And user colors "STARTED" column linked to "CONTROL" column
    Then "SEX" column should be color-coded linked
    And "DIS_POP" column should be color-coded linked
    And "CONTROL" column should be color-coded linked
    And "STARTED" column should be color-coded linked
    And the coloring of "STARTED" column should be linked to "CONTROL" column
    And the "color of cell 1 of SEX" reading of grid should be "#00ff00"
    And the "color of cell 3 of SEX" reading of grid should be "#ff0000"
    When user removes the coloring of "SEX" column
    And user removes the coloring of "DIS_POP" column
    And user removes the coloring of "CONTROL" column
    And user removes the coloring of "STARTED" column
    Then no errors should have been logged

  Scenario: Inverting a linear scheme swaps the ends of the gradient
    When user colors "AGE" column linearly from "#0000FF" to "#FF0000"
    Then the color scheme of "AGE" column should be "#0000FF, #FF0000"
    And the "color of cell 1 of AGE" and "color of cell 3 of AGE" readings of grid should differ
    When user remembers the "color of cell 1 of AGE" reading of grid
    And user inverts the color scheme of "AGE" column
    Then the color scheme of "AGE" column should be "#FF0000, #0000FF"
    And grid should have repainted
    And the "color of cell 1 of AGE" reading of grid should not be as remembered
    When user removes the coloring of "AGE" column
    Then "AGE" column should have no color coding
    And no errors should have been logged
