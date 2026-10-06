@journey @viewers @realizes:viewers.scatter-plot
Feature: Scatter plot selection and viewport navigation
  Selecting markers and moving the viewport: a Shift-drag selects the rows inside the rectangle, a
  Control+Shift-drag takes them out again, a click on marker-free space clears the selection, and a
  jitter change leaves it alone. The viewport answers to an Alt-drag, a plain drag, the wheel and
  the axis range slider, and Reset View, a double-click on empty space and the H key all bring it
  home. Ctrl+A selects what the filter passes, Ctrl+Shift+A and Escape clear it, and L turns on the
  lasso. A zoom is animated over a few frames and must not touch the selection, which is what the
  `rows selected` claim after each zoom says. Read from what the plot reports: `rows selected`,
  `x axis min` / `x axis span`, `hovered row`, and the `empty space` and `marker of row <n>` hit
  areas. One journey on demog-1000, X = WEIGHT, Y = HEIGHT (872 of the 1000 rows have a HEIGHT);
  every scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | WEIGHT |
      | Y | HEIGHT |
    Then scatter plot viewer should show 872 rows
    And no rows should be selected

  Scenario: A drag selects the markers inside it, a second drag adds, Control+Shift takes them out
    Then "Reset Selection On Background Click" property of scatter plot viewer should be "true"
    When user drags a selection box from the "marker of row 11" area to the "marker of row 3" area of scatter plot viewer
    Then some rows should be selected
    And the "rows selected" reading of scatter plot viewer should be higher than before
    When user drags a selection box over the "view" area of scatter plot viewer
    Then the "rows selected" reading of scatter plot viewer should be higher than before
    And every selected row should pass the filter
    When user drags a deselection box over the "view" area of scatter plot viewer
    Then the "rows selected" reading of scatter plot viewer should be lower than before
    When user drags a selection box over the "view" area of scatter plot viewer
    Then some rows should be selected
    When user clicks on the "empty space" area of scatter plot viewer
    Then no rows should be selected
    And no errors should have been logged

  Scenario: The selection survives a jitter change
    When user sets properties of scatter plot viewer:
      | Jitter Size   | 20 |
      | Jitter Size Y | 15 |
    And user drags a selection box over the "view" area of scatter plot viewer
    Then some rows should be selected
    When user remembers the "rows selected" reading of scatter plot viewer
    And user sets "Jitter Size" property of scatter plot viewer to "30"
    Then scatter plot viewer should have repainted
    And the "rows selected" reading of scatter plot viewer should be as remembered
    When user drags a deselection box over the "view" area of scatter plot viewer
    Then the "rows selected" reading of scatter plot viewer should be lower than before
    When user clears the row selection
    And user sets properties of scatter plot viewer:
      | Jitter Size   | 0 |
      | Jitter Size Y | 0 |
    Then no rows should be selected
    And no errors should have been logged

  Scenario: An Alt-drag zooms, a plain drag pans, the wheel zooms and Reset View comes home
    When user selects the first 20 rows
    Then 20 rows should be selected
    When user sets "Zoom and Filter" property of scatter plot viewer to "no action"
    And user picks "Reset View" from the context menu of scatter plot viewer
    And user remembers the value range of scatter plot viewer
    And user remembers the "x axis span" reading of scatter plot viewer
    And user drags a zoom box over the "view" area of scatter plot viewer
    Then the "rows selected" reading of scatter plot viewer should be the same as before
    And scatter plot viewer should show a narrower value range than before
    And the "x axis span" reading of scatter plot viewer should be lower than before
    And all rows should pass the filter
    When user picks "Reset View" from the context menu of scatter plot viewer
    Then scatter plot viewer should show the remembered value range
    And the "x axis span" reading of scatter plot viewer should be as remembered
    When user drags across the "view" area of scatter plot viewer
    Then the "x axis min" reading of scatter plot viewer should differ from before
    And the "x axis span" reading of scatter plot viewer should be the same as before
    When user picks "Reset View" from the context menu of scatter plot viewer
    Then scatter plot viewer should show the remembered value range
    When user scrolls the mouse wheel up over the "view" area of scatter plot viewer
    Then the "rows selected" reading of scatter plot viewer should be the same as before
    And scatter plot viewer should show a narrower value range than before
    And the "x axis span" reading of scatter plot viewer should be lower than before
    When user picks "Reset View" from the context menu of scatter plot viewer
    Then scatter plot viewer should show the remembered value range
    And all rows should pass the filter
    When user sets "Zoom and Filter" property of scatter plot viewer to "filter by zoom"
    And user clears the row selection
    Then no rows should be selected
    And no errors should have been logged

  Scenario: The X range slider narrows the axis and a double-click on empty space resets the view
    When user sets "Zoom and Filter" property of scatter plot viewer to "no action"
    And user picks "Reset View" from the context menu of scatter plot viewer
    And user remembers the value range of scatter plot viewer
    And user takes a snapshot of scatter plot viewer
    And user drags the min handle of the "x" range slider of scatter plot viewer by 30 pixels
    Then the "x axis min" reading of scatter plot viewer should be higher than before
    And the "x axis span" reading of scatter plot viewer should be lower than before
    When user picks "Reset View" from the context menu of scatter plot viewer
    Then scatter plot viewer should show the remembered value range
    When user scrolls the mouse wheel up over the "view" area of scatter plot viewer
    Then the "rows selected" reading of scatter plot viewer should be the same as before
    And scatter plot viewer should show a narrower value range than before
    When user listens for "d4-scatterplot-reset-view" event on scatter plot viewer
    And user hovers over the "empty space" area of scatter plot viewer
    Then the "hovered row" reading of scatter plot viewer should be 0
    When user double-clicks on the "empty space" area of scatter plot viewer
    Then "d4-scatterplot-reset-view" event should have fired on scatter plot viewer
    And scatter plot viewer should show the remembered value range
    When user sets "Zoom and Filter" property of scatter plot viewer to "filter by zoom"
    Then no errors should have been logged

  Scenario: Ctrl+A selects the rows the filter passes, Ctrl+Shift+A and Escape clear them
    When user filters rows where "SEX" is "F"
    Then 553 rows should pass the filter
    When user clicks on the "empty space" area of scatter plot viewer
    And user presses Control+a
    Then 553 rows should be selected
    And only rows where "SEX" is "F" should be selected
    When user presses Control+Shift+a
    Then no rows should be selected
    When user drags a selection box over the "view" area of scatter plot viewer
    Then some rows should be selected
    When user presses Escape
    Then no rows should be selected
    When user resets the filter
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  Scenario: H restores the viewport after a wheel zoom
    When user picks "Reset View" from the context menu of scatter plot viewer
    And user remembers the value range of scatter plot viewer
    And user scrolls the mouse wheel up over the "view" area of scatter plot viewer
    Then the "rows selected" reading of scatter plot viewer should be the same as before
    And scatter plot viewer should show a narrower value range than before
    When user clicks on the "empty space" area of scatter plot viewer
    And user presses h
    Then scatter plot viewer should show the remembered value range
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: L turns on the Lasso Tool and a lasso selects the markers it encloses
    Then "Lasso Tool" property of scatter plot viewer should be "false"
    When user clicks on the "empty space" area of scatter plot viewer
    And user presses l
    Then "Lasso Tool" property of scatter plot viewer should be "true"
    When user drags a lasso over the "view" area of scatter plot viewer
    Then some rows should be selected
    And every selected row should pass the filter
    When user presses l
    Then "Lasso Tool" property of scatter plot viewer should be "false"
    When user clicks on the "empty space" area of scatter plot viewer
    Then no rows should be selected
    And no errors should have been logged
