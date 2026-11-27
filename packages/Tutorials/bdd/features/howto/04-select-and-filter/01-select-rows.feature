@guide @tutorial
Feature: How-to · Select and filter
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). Selecting rows: a click on a cell makes its row current,
  Shift+click selects the rows from there to the clicked one, Ctrl+click adds one more, the status bar counts
  them, Select > Invert swaps selected and unselected, Esc on the grid clears the selection, and a Shift+drag on
  a scatter plot selects the points inside the box, which the grid shows too. demog (5850 rows) and a scatter
  plot of AGE and WEIGHT are opened off camera; nothing is made on the server.

  Scenario: Select rows
    Click, Shift+click and Ctrl+click in the grid, use the Select menu, draw a box on a plot, clear with Esc.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And user adds a scatter plot viewer with:
      | X | AGE    |
      | Y | WEIGHT |
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    Then no rows should be selected
    # caption: A click makes a row current, it does not select it
    When user clicks on the "cell 3 of AGE" area of grid
    Then row 3 should be current
    And no rows should be selected
    # caption: Shift+click selects every row in between
    When user clicks on the "cell 8 of AGE" area of grid holding Shift
    Then 6 rows should be selected
    And rows 3 to 8 should all be selected
    # caption: Ctrl+click adds a row, or takes it out
    When user clicks on the "cell 11 of AGE" area of grid holding Control
    Then 7 rows should be selected
    When user hovers over "Selected: 7" text in status bar
    And the video says "The status bar counts the selected rows"
    # caption: The Select menu: all, none, invert and more
    When user picks "Select > Invert" from the top menu
    Then 5843 rows should be selected
    # caption: Esc clears the selection
    When user presses Escape in grid
    Then no rows should be selected
    # caption: On a plot, Shift+drag a box around the points
    When user remembers the "rows selected" reading of scatter plot viewer
    And user drags a selection box from the "marker of row 11" area to the "marker of row 13" area of scatter plot viewer
    Then the "rows selected" reading of scatter plot viewer should be higher than remembered
    And some rows should be selected
    When the video says "The same rows are selected in the grid and in every viewer"
    Then every selected row should pass the filter
    And no errors should have been logged
