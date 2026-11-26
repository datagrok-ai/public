@guide @tutorial
Feature: How-to · Viewers
  For the maintainer: the grid and every viewer of a table share its rows. A hovered histogram bin lights
  its rows on the scatter plot (the mouse-over group); a row clicked in the grid becomes the current row
  in every viewer; a Shift-drag across histogram bins selects their rows, and the scatter plot and the
  bar chart show the selection; a click on a bar selects that bar's rows instead (DIS_POP = RA, 2,550 rows); and
  the filter panel's RA category filters the grid and every viewer at once (one bar left). demog is
  opened and the three viewers are added off camera; nothing is made on the server.

  Scenario: Viewers work together
    Hover, select, pick the current row or filter in one place, and the grid and every other viewer follow.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user adds a scatter plot viewer
    And user adds a histogram viewer with:
      | Value | AGE |
    And user adds a bar chart viewer
    When the video says "The grid and all viewers of a table show the same rows"
    # caption: Point at a bin: its rows light up on the other viewers
    And user takes a snapshot of scatter plot viewer
    And user hovers over the "bin 8" area of histogram viewer
    Then scatter plot viewer should have repainted
    # caption: Click a row to make it current everywhere
    When user clicks on the "cell 10 of AGE" area of grid
    Then row 10 should be current
    And histogram viewer should have a "current row marker" area
    # caption: Shift+drag across bins selects their rows, shown everywhere
    When user drags a selection box from the "bin 4" area to the "bin 7" area of histogram viewer
    Then some rows should be selected
    And scatter plot viewer should show a selection highlight
    And bar chart viewer should show a selection highlight
    # caption: Click a bar to select its rows
    When user clicks on the "bar RA" area of bar chart viewer
    Then only rows where "DIS_POP" is "RA" should be selected
    And 2550 rows should be selected
    # caption: A filter applies to the grid and every viewer
    When user clicks on filter icon in toolbar
    And user clicks on the "category Psoriasis of DIS_POP" area of filter panel
    Then 1209 rows should pass the filter
    And histogram viewer should show 1209 rows
    And the "bars" reading of bar chart viewer should be 1
    And no errors should have been logged
