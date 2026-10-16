@guide
Feature: Show Legend of curves in grid cells is set for the table, a column or one cell
  A guide: the answer to "how is the curve legend switched off and on in the grid?" (Curves 1.13.1+),
  filmed as a help-mode GIF for a community post. A curve cell names its curves in a legend once it is
  at least 325 x 275 px. Clicking it shows Context Panel > Options: Level (Dataframe, Column or Cell,
  Column by default) says where the options below it are set, and Chart options > Show Legend turns
  the legend off or on. Set for the column, it changes every cell of the column; set for one cell, only
  that cell. Demo: the Assay Curves table of the Curves package, one named curve per cell.

  Scenario: Hide the legend for the column, then show it in one cell
    Given user is logged in
    And simple mode is off
    And user opens assay-curves dataset
    And the grid of the current view does not show the "Molecule" column
    And the grid of the current view has rows 290 px high and the "Fitted Curve" column 340 px wide
    And the context panel is open
    When user clicks on the "cell 1 of Fitted Curve" area of grid
    And user expands "Options" pane in context panel
    Then "Level" input in context panel should have value "Column"
    And every curve cell in view of "Fitted Curve" column should draw a legend
    When user scrolls "Show Legend" input in context panel into view
    And user unchecks "Show Legend" input in context panel
    Then no curve cell in view of "Fitted Curve" column should draw a legend
    When user scrolls "Level" input in context panel into view
    And user selects "Cell" in "Level" input in context panel
    And user scrolls "Show Legend" input in context panel into view
    And user checks "Show Legend" input in context panel
    Then only the current curve cell of "Fitted Curve" column should draw a legend
