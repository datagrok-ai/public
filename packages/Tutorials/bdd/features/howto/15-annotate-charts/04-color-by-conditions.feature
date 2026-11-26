@guide @tutorial
Feature: How-to · Annotations
  For the maintainer: conditional color coding, on demog-1000 with a scatter plot of WEIGHT against HEIGHT colored
  by AGE (a linear scale to begin with). The AGE header's context menu, Color Coding > Conditional, cuts AGE into
  four equal ranges (18 to 89), and the plot's legend lists them; the item is a radio choice, which leaves the menu
  open, so Escape closes it. Color Coding > Edit... opens the "Color-coding: AGE" editor: Bin... asks for the
  number of ranges (3, cut at 41.67 and 65.33); the middle row's color bar ("Color 2") opens a color picker titled
  with that row's range, where grey is picked (before the patterns are typed: the picker keeps the title it was
  made with), and each row's pattern is typed over: <40, 40-65, >65 (of the 872 plotted rows: 280, 523
  and 69). CLOSE keeps it; the grid's AGE cells and the plot take the same three colors and the legend lists the
  three patterns. Coloring a column by thresholds is not in the Basics clip of color coding (linear and
  categorical only). demog-1000 is opened and the plot added off camera; the coloring stays in the browser,
  nothing is made on the server.

  Scenario: Color points by conditions
    Give each range of a column its own color; every chart colored by that column follows.

    Given user is logged in
    And simple mode is off
    And user opens demog-1000 dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user adds a scatter plot viewer with:
      | X     | WEIGHT |
      | Y     | HEIGHT |
      | Color | AGE    |
    And scatter plot viewer takes 60% of the space it shares with the grid
    When the video says "Conditional colors: one color for each range of values"
    # caption: Right-click the column's header: Color Coding > Conditional
    And user picks "Color Coding > Conditional" from the context menu of the "header AGE" area of grid
    And user presses Escape
    # caption: The plot takes the column's colors, one per range
    Then "AGE" column should be color-coded conditionally
    And the legend of scatter plot viewer should list 4 items
    # caption: Color Coding > Edit... sets the ranges
    When user picks "Color Coding > Edit..." from the context menu of the "header AGE" area of grid
    Then "Color-coding: AGE" dialog should be visible
    # caption: Bin...: how many ranges
    When user clicks on "Bin..." link in "Color-coding: AGE" dialog
    Then "Edit" dialog should be visible
    When user enters "3" into "Number of bins" input in "Edit" dialog
    And user clicks on OK button in "Edit" dialog
    Then "Edit" dialog should be hidden
    # caption: Click a color to change it
    When user clicks on "Color 2" button in "Color-coding: AGE" dialog
    Then "41.67 - 65.33" dialog should be visible
    When user picks the color "#c7c7c7" in the color picker dialog
    And user clicks on OK button in "41.67 - 65.33" dialog
    # caption: Type each range: under 40, 40 to 65, over 65
    When user enters "<40" into "Pattern 1" text input in "Color-coding: AGE" dialog
    And user enters "40-65" into "Pattern 2" text input in "Color-coding: AGE" dialog
    And user enters ">65" into "Pattern 3" text input in "Color-coding: AGE" dialog
    And user clicks on CLOSE button in "Color-coding: AGE" dialog
    # caption: The grid and the plot follow the new ranges
    Then the legend of scatter plot viewer should list 3 items
    And the "40-65" item in the legend of scatter plot viewer should be colored "#c7c7c7"
    And "AGE" column should be color-coded conditionally
    And the "color of cell 1 of AGE" and "color of cell 3 of AGE" readings of grid should differ
    And no errors should have been logged
