@guide @tutorial
Feature: How-to · Viewers
  For the maintainer: what the legend of a viewer is for, on a scatter plot of demog-1000 colored by
  RACE (872 rows have a WEIGHT and a HEIGHT: Caucasian 775, Other 59, Black 25, Asian 13). A click on
  "Other" shows only those 59 on this plot while every row still passes the table's filter;
  Ctrl+click on "Black" adds its 25 (84); the cross on "Black" takes it back (59); a second click on
  "Other" shows all 872 again. A hovered item shows a palette icon: it opens the category's color
  dialog, and the picked purple recolors Caucasian on the legend and the plot (the column's color
  coding, in the browser). The plot is added with its legend docked on the right and given most of the
  view off camera; demog-1000 is opened off camera; nothing is made on the server.

  Scenario: Use the legend
    Click legend items to show only some categories on a viewer, and change a category's color.

    Given user is logged in
    And simple mode is off
    And user opens demog-1000 dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user adds a scatter plot viewer with:
      | X                 | WEIGHT |
      | Y                 | HEIGHT |
      | Color             | RACE   |
      | Legend Visibility | Always |
      | Legend Position   | Right  |
    And scatter plot viewer takes 65% of the space it shares with the grid
    Then the legend of scatter plot viewer should list 4 items
    When the video says "The legend names the colors on a viewer"
    # caption: Click a category to show only it on this plot
    And user clicks on "Other" legend item in legend of scatter plot viewer
    Then scatter plot viewer should show 59 rows
    And "Other" legend item in legend of scatter plot viewer should be selected
    # caption: The table keeps all its rows
    And all rows should pass the filter
    # caption: Ctrl+click adds another category
    When user clicks on "Black" legend item in legend of scatter plot viewer holding Control
    Then scatter plot viewer should show 84 rows
    # caption: The x takes one back; a click on the last one shows all
    When user clicks on cross of "Black" legend item in legend of scatter plot viewer
    Then scatter plot viewer should show 59 rows
    When user clicks on "Other" legend item in legend of scatter plot viewer
    Then scatter plot viewer should show 872 rows
    # caption: Point at an item and click the palette to recolor it
    When user hovers over "Caucasian" legend item in legend of scatter plot viewer
    And user clicks on color picker icon
    Then "Caucasian" dialog should be visible
    When user picks the color "#9467BD" in the color picker dialog
    And user clicks on OK button in "Caucasian" dialog
    Then the "Caucasian" item in the legend of scatter plot viewer should be colored "#9467BD"
    And the "view" area of scatter plot viewer should contain the color "#9467BD"
    And no errors should have been logged
