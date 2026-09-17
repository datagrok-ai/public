@journey @viewers @realizes:viewers.grid
Feature: Grid cell appearance
  What decides the colour and the text of a cell: the per-column colour coding picked from the
  header menu (Linear, Conditional, Categorical, Linked), the grid-wide Color Coding that overrides
  it, a custom number format, the missing-value colour, the row height and the default cell font.
  The grid reports the renderer-resolved colour and text of every visible cell as the readings
  `color of cell <r> of <c>` and `text of cell <r> of <c>`, so a claim about a cell reads the cell.
  One journey on demog-1000 (AGE 26, 30, 58, 45 in rows 1-4; HEIGHT 174.705, 150.288, 174.066,
  183.83; the first missing HEIGHT is row 298).

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

  Scenario: Narrowing a coloured column keeps its colour and its text
    When user picks "Grid Color Coding > All" from the context menu of the "cell 1 of AGE" area of grid
    And user remembers the "color of cell 1 of AGE" reading of grid
    And user drags the "column resizer AGE" area of grid by 30 pixels to the left
    Then the "column width of AGE" reading of grid should be lower than before
    And the "color of cell 1 of AGE" reading of grid should be as remembered
    And the "text of cell 1 of AGE" reading of grid should be "26"
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

  Scenario: A bigger cell font repaints the grid and an idle grid does not
    When user takes a snapshot of grid
    Then grid should not have repainted
    When user sets "Default Cell Font" property of grid to "20px Roboto"
    Then grid should have repainted by at least 3000 pixels
    When user sets "Default Cell Font" property of grid to "12px Roboto"
    Then grid should have repainted
    And no errors should have been logged
