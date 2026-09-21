@guide @help:visualize/viewers
Feature: Color a scatter plot by conditions, then change one condition's color
  A guide: the answer to "how do I get conditional color coding on a scatter plot, and then change
  the color of one of the conditions?". A scatter plot takes its colors from the column set as its
  Color, so the coding is switched on the column: the column header's context menu in the grid,
  Color Coding > Conditional, cuts a numeric column into ranges (four by default) and the plot's
  legend lists one entry per range. A range's color is the palette icon a hovered legend entry
  shows: it opens the color dialog for that range, and OK keeps the pick — on the plot, in the
  legend and in the grid column alike. Demo: demog-1000, WEIGHT by HEIGHT, colored by AGE (18 to 89).

  Scenario: Switch the color column to conditional coding and recolor one range from the legend
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X     | WEIGHT |
      | Y     | HEIGHT |
      | Color | AGE    |
    Then scatter plot viewer should have a "color scale" area
    When user picks "Color Coding > Conditional" from the context menu of the "header AGE" area of grid
    And user presses Escape
    Then "AGE" column should be color-coded conditionally
    And context menu should be hidden
    And the legend of scatter plot viewer should list 4 items
    When user hovers over "18 - 35.75" legend item in legend of scatter plot viewer
    And user clicks on color picker icon
    Then "18 - 35.75" dialog should be visible
    When user picks the color "#9467BD" in the color picker dialog
    And user clicks on OK button in "18 - 35.75" dialog
    Then the "18 - 35.75" item in the legend of scatter plot viewer should be colored "#9467BD"
    And the "view" area of scatter plot viewer should contain the color "#9467BD"
