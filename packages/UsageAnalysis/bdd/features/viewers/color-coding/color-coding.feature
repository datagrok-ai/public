@journey @viewers @realizes:viewers.grid
Feature: Colour coding across columns
  What a colouring keeps when it is switched off, copied to another column, linked to another
  column's, or inverted. The types themselves (Linear, Conditional, Categorical from the header
  menu) and the grid-wide Color Coding are covered by `grid/grid-appearance.feature`; this journey
  is about the state that survives an edit — the md's boolean CONTROL on default categorical colours,
  STARTED on a three-stop linear scheme, AGE and STARTED switched off and back on with their schemes
  kept, and a categorical colouring (RACE to a copy of it) and a linear scheme (STARTED to HEIGHT)
  applied to other columns (through the colouring's own copy, the step the md's Pick Up / Apply
  performs, not the header menu) — and a layout and a project: a categorical colouring
  comes back when the saved layout is applied to the table opened anew, and a linked colouring (a
  link to a link, RACE to AGE, HEIGHT to RACE) comes back from a layout and from a project.
  Not translated: the md's "Grid Color Coding" header-menu ladder (All / scheme / None / Auto),
  which is `grid/grid-appearance.feature`'s, and the SPGI_v2 Edit Color Scheme dialog, which the
  plan for this round does not include.
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

  Scenario: A boolean column takes the default categorical colours, a date column a three-stop scheme
    When user colors "CONTROL" column categorically again
    Then "CONTROL" column should be color-coded categorically
    And the "color of cell 1 of CONTROL" and "color of cell 3 of CONTROL" readings of grid should differ
    When user colors "STARTED" column linearly through "#0000FF, #FFFFFF, #FF0000"
    Then "STARTED" column should be color-coded linearly
    And the color scheme of "STARTED" column should be "#0000FF, #FFFFFF, #FF0000"
    And the "color of cell 1 of STARTED" and "color of cell 2 of STARTED" readings of grid should differ
    When user removes the coloring of "CONTROL" column
    Then "CONTROL" column should have no color coding
    And no errors should have been logged

  Scenario: AGE and STARTED switched off and back on keep the schemes they were given
    When user colors "AGE" column conditionally:
      | <30   | #00FF00 |
      | 30-60 | #FFFF00 |
      | >60   | #FF0000 |
    Then the "color of cell 1 of AGE" reading of grid should be "#00ff00"
    When user removes the coloring of "AGE" column
    And user removes the coloring of "STARTED" column
    Then "AGE" column should have no color coding
    And "STARTED" column should have no color coding
    And the "color of cell 1 of STARTED" and "color of cell 2 of STARTED" readings of grid should be the same
    When user colors "AGE" column conditionally again
    And user colors "STARTED" column linearly again
    Then "AGE" column should be color-coded conditionally
    And the "color of cell 1 of AGE" reading of grid should be "#00ff00"
    And the "color of cell 3 of AGE" reading of grid should be "#ffff00"
    And "STARTED" column should be color-coded linearly
    And the color scheme of "STARTED" column should be "#0000FF, #FFFFFF, #FF0000"
    And the "color of cell 1 of STARTED" and "color of cell 2 of STARTED" readings of grid should differ
    And no errors should have been logged

  Scenario: Pick Up / Apply copies a categorical colouring and a linear scheme to other columns
    When user adds a calculated column "Race_copy" with formula "${RACE}"
    Then "Race_copy" column should have no color coding
    And "HEIGHT" column should have no color coding
    When user colors "RACE" column categorically:
      | Caucasian | #3366CC |
      | Other     | #CC6699 |
    And user applies the coloring of "RACE" column to "Race_copy" column
    Then "Race_copy" column should be color-coded categorically
    And the categorical color of "Caucasian" in "Race_copy" column should be "3366CC"
    And the categorical color of "Other" in "Race_copy" column should be "CC6699"
    And the "color of cell 1 of Race_copy" and "color of cell 1 of RACE" readings of grid should be the same
    When user applies the coloring of "STARTED" column to "HEIGHT" column
    Then "HEIGHT" column should be color-coded linearly
    And the color scheme of "HEIGHT" column should be "#0000FF, #FFFFFF, #FF0000"
    And the "color of cell 1 of HEIGHT" and "color of cell 2 of HEIGHT" readings of grid should differ
    When user removes the coloring of "HEIGHT" column
    And user removes the coloring of "RACE" column
    And user removes the coloring of "STARTED" column
    And user removes the coloring of "AGE" column
    And user removes "Race_copy" column
    Then the table should not have a column "Race_copy"
    And no errors should have been logged

  Scenario: A column's colouring comes back from a saved layout on the reopened table
    When user colors "SEX" column categorically:
      | M | #3366CC |
      | F | #CC6699 |
    And user adds a scatter plot viewer
    And user saves the layout of the current table view to the server
    And user closes all views
    And user opens demog-1000 dataset
    Then "SEX" column should have no color coding
    When user loads the saved layout
    Then scatter plot viewer should be visible
    And "SEX" column should be color-coded categorically
    And the categorical color of "M" in "SEX" column should be "3366CC"
    And the categorical color of "F" in "SEX" column should be "CC6699"
    And the "color of cell 1 of SEX" and "color of cell 4 of SEX" readings of grid should differ
    When user clicks on close icon of scatter plot viewer
    Then "SEX" column should be color-coded categorically
    When user removes the coloring of "SEX" column
    Then no errors should have been logged

  Scenario: Linked colourings come back from a saved layout and from a saved project
    When user colors "AGE" column conditionally:
      | <40 | #00FF00 |
      | >40 | #FF0000 |
    And user colors "RACE" column linked to "AGE" column
    And user colors "HEIGHT" column linked to "RACE" column
    Then the coloring of "RACE" column should be linked to "AGE" column
    And the coloring of "HEIGHT" column should be linked to "RACE" column
    When user saves the layout of the current table view to the server
    And user removes the coloring of "HEIGHT" column
    And user removes the coloring of "RACE" column
    Then "RACE" column should have no color coding
    And "HEIGHT" column should have no color coding
    When user loads the saved layout
    Then "RACE" column should be color-coded linked
    And the coloring of "RACE" column should be linked to "AGE" column
    And "HEIGHT" column should be color-coded linked
    And the coloring of "HEIGHT" column should be linked to "RACE" column
    And the "color of cell 1 of RACE" reading of grid should be "#00ff00"
    And the "color of cell 3 of RACE" reading of grid should be "#ff0000"
    When user saves the current view as project "bdd color coding links"
    And user closes all views
    And user opens the "bdd color coding links" project
    Then "RACE" column should be color-coded linked
    And the coloring of "RACE" column should be linked to "AGE" column
    And "HEIGHT" column should be color-coded linked
    And the coloring of "HEIGHT" column should be linked to "RACE" column
    And the "color of cell 1 of RACE" reading of grid should be "#00ff00"
    And the "color of cell 3 of RACE" reading of grid should be "#ff0000"
    When user removes the coloring of "HEIGHT" column
    And user removes the coloring of "RACE" column
    And user removes the coloring of "AGE" column
    Then "RACE" column should have no color coding
    And no errors should have been logged
