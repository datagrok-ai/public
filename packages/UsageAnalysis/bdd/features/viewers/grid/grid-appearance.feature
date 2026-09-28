@journey @viewers @realizes:viewers.grid
Feature: Grid cell appearance
  What decides the colour and the text of a cell: the per-column colour coding picked from the
  header menu (Linear, Conditional, Categorical, Linked), the grid-wide Color Coding that overrides
  it, a custom number format (shown in the cells and in the `format` row of the Context Panel's
  Details), the missing-value colour, the row height, the default cell font, the Selected Rows
  Color, and a column's Style > Content background, which its colour coding overrides
  (GROK-18638). The grid reports the renderer-resolved colour and text of every visible cell as the
  readings `color of cell <r> of <c>` and `text of cell <r> of <c>`, so a claim about a cell reads
  the cell - except the selection, which the grid paints over the cells on its overlay: the
  Selected Rows Color is claimed as a colour inside a cell's area (the library composites the
  overlay into the picture), with the same area checked without it before and after. The style
  scenario uses SEVERITY and Categorical where the TestTrack spec uses AGE and Linear, since AGE is
  already Linear in this journey. One journey on demog-1000 (AGE 26, 30, 58, 45 in rows 1-4; HEIGHT
  174.705, 150.288, 174.066, 183.83; the first missing HEIGHT is row 298).

  Not translated, and why: from `grid-cell-appearance.md` scenario 4, that a narrowed column keeps
  the full formatted value - `text of cell` is the value whatever the width, so the claim could not
  fail - and the format shown in the HEIGHT column's panel - the tag is AGE's, and AGE's own Details
  show it.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then grid should show 1000 rows
    And "Color Coding" property of grid should be "Auto"
    And "HEIGHT" column should have missing values

  Scenario: Linear coding from the header menu paints the column by value
    When user picks "Color Coding > Linear" from the context menu of the "header AGE" area of grid
    Then "AGE" column should have tag ".color-coding-type" equal to "Linear"
    And grid should have repainted
    When user moves the pointer away from grid
    Then the "cell 2 of AGE" and "cell 3 of AGE" areas of grid should be painted in different colors
    And the "cell 2 of AGE" and "cell 2 of USUBJID" areas of grid should be painted in different colors
    And no errors should have been logged

  Scenario: The grid-wide Color Coding overrides the column's
    When user picks "Grid Color Coding > None" from the context menu of the "cell 8 of AGE" area of grid
    And user moves the pointer away from grid
    Then "Color Coding" property of grid should be "None"
    And the "cell 2 of AGE" and "cell 3 of AGE" areas of grid should be painted in the same colors
    When user picks "Grid Color Coding > All" from the context menu of the "cell 8 of AGE" area of grid
    And user moves the pointer away from grid
    Then "Color Coding" property of grid should be "All"
    And the "cell 2 of HEIGHT" and "cell 4 of HEIGHT" areas of grid should be painted in different colors
    And the "cell 2 of AGE" and "cell 3 of AGE" areas of grid should be painted in different colors
    When user picks "Grid Color Coding > Auto" from the context menu of the "cell 8 of AGE" area of grid
    And user moves the pointer away from grid
    Then "Color Coding" property of grid should be "Auto"
    And the "cell 2 of AGE" and "cell 3 of AGE" areas of grid should be painted in different colors
    And the "cell 2 of HEIGHT" and "cell 4 of HEIGHT" areas of grid should be painted in the same colors
    And no errors should have been logged

  Scenario: Narrowing a coloured column keeps its colour
    When user picks "Grid Color Coding > All" from the context menu of the "cell 1 of AGE" area of grid
    And user remembers the "color of cell 1 of AGE" reading of grid
    And user drags the "column resizer AGE" area of grid by 30 pixels to the left
    Then the "column width of AGE" reading of grid should be lower than before
    And the "color of cell 1 of AGE" reading of grid should be as remembered
    When user picks "Column Sizing > Optimal" from the context menu of the "cell 1 of AGE" area of grid
    And user picks "Grid Color Coding > Auto" from the context menu of the "cell 1 of AGE" area of grid
    Then no errors should have been logged

  Scenario: A custom format shows in the cell
    When user picks "Format > Custom..." from the context menu of the "header AGE" area of grid
    Then Format AGE dialog should be visible
    When user enters "0.00" into Custom input in Format AGE dialog
    And user clicks on OK button in Format AGE dialog
    Then Format AGE dialog should be hidden
    And "AGE" column should have tag "format" equal to "0.00"
    And the "text of cell 1 of AGE" reading of grid should be "26.00"
    When user clicks on the "header AGE" area of grid
    Given the context panel is open
    Then the context panel should show "AGE"
    Given Details accordion header in context panel is expanded
    Then "format" table row in context panel should contain text "0.00"
    When user picks "Format > Custom..." from the context menu of the "header AGE" area of grid
    And user clears Custom input in Format AGE dialog
    And user clicks on OK button in Format AGE dialog
    Then the "text of cell 1 of AGE" reading of grid should be "26"
    And no errors should have been logged

  Scenario: Conditional coding paints the cells its ranges name
    When user picks "Color Coding > Conditional" from the context menu of the "header HEIGHT" area of grid
    Then "HEIGHT" column should have tag ".color-coding-type" equal to "Conditional"
    When user colors "HEIGHT" column conditionally:
      | <160 | #0000FF |
      | >180 | #FF0000 |
    Then the "cell 2 of HEIGHT" area of grid should contain the color "#0000FF"
    And the "cell 4 of HEIGHT" area of grid should contain the color "#FF0000"
    And the "cell 1 of HEIGHT" area of grid should not contain the color "#0000FF"
    And the "cell 1 of HEIGHT" area of grid should not contain the color "#FF0000"
    When user removes the coloring of "HEIGHT" column
    Then "HEIGHT" column should have no color coding
    And no errors should have been logged

  Scenario: Categorical coding tells the categories apart
    When user picks "Color Coding > Categorical" from the context menu of the "header SEX" area of grid
    Then "SEX" column should be color-coded categorically
    And the "cell 1 of SEX" and "cell 4 of SEX" areas of grid should be painted in different colors
    When user removes the coloring of "SEX" column
    Then "SEX" column should have no color coding
    And no errors should have been logged

  Scenario: Row height and the missing-value colour paint through the grid's own properties
    When user sets "Row Height" property of grid to "48"
    Then the "cell 1 of AGE" area of grid should be at least 40 pixels tall
    And grid should have repainted
    When user sets "Row Height" property of grid to "28"
    And user makes row 298 current
    And user makes row 300 current
    And user sets "Missing Value Color" property of grid to "#FFAAAA"
    Then the "cell 298 of HEIGHT" area of grid should contain the color "#FFAAAA"
    When user sets "Missing Value Color" property of grid to "#00FF00"
    Then the "cell 298 of HEIGHT" area of grid should contain the color "#00FF00"
    And the "cell 298 of HEIGHT" area of grid should not contain the color "#FFAAAA"
    When user sets "Missing Value Color" property of grid to "#FFFFFF"
    And user makes row 1 current
    Then no errors should have been logged

  Scenario: A bigger cell font repaints the grid
    When user sets "Default Cell Font" property of grid to "20px Roboto"
    Then grid should have repainted by at least 3000 pixels
    When user sets "Default Cell Font" property of grid to "12px Roboto"
    Then grid should have repainted
    And no errors should have been logged

  Scenario: Selected Rows Color paints the selected rows and Escape takes it away
    When user presses Escape
    Then no rows should be selected
    When user sets "Selected Rows Color" property of grid to "#00FF00"
    Then the "cell 6 of USUBJID" area of grid should not contain the color "#00FF00"
    When user drags a selection box from the "row header 5" area to the "row header 7" area of grid
    Then rows 5 to 7 should be selected
    And the "cell 6 of USUBJID" area of grid should contain the color "#00FF00"
    And the "cell 9 of USUBJID" area of grid should not contain the color "#00FF00"
    When user presses Escape
    Then no rows should be selected
    And the "cell 6 of USUBJID" area of grid should not contain the color "#00FF00"
    When user sets "Selected Rows Color" property of grid to "819780688"
    Then no errors should have been logged

  Scenario: A column's colour coding wins over the background its Style sets (GROK-18638)
    When user clicks on the "header SEVERITY" area of grid
    Given the context panel is open
    Then the context panel should show "SEVERITY"
    Given Style accordion header in context panel is expanded
    And Content accordion header in context panel is expanded
    When user clicks on editor of "Back Color" property in context panel
    And user picks the color "#FFA500" in the color picker
    Then the "color of cell 1 of SEVERITY" reading of grid should be "#ffa500"
    And the "color of cell 2 of SEVERITY" reading of grid should be "#ffa500"
    When user picks "Color Coding > Categorical" from the context menu of the "header SEVERITY" area of grid
    Then the "color of cell 1 of SEVERITY" reading of grid should not be "#ffa500"
    And the "color of cell 1 of SEVERITY" and "color of cell 3 of SEVERITY" readings of grid should differ
    When user removes the coloring of "SEVERITY" column
    Then the "color of cell 1 of SEVERITY" reading of grid should be "#ffa500"
    And no errors should have been logged
