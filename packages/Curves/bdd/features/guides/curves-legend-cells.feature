@guide
Feature: Show Legend of curves in grid cells is set for the table, a column or one cell
  A guide: the answer to "how is the curve legend switched off and on in the grid?" (Curves 1.13.1+),
  filmed as a help-mode GIF for a community post. A curve cell names its curves in a legend once it is
  at least 325 x 275 px. Clicking it shows Context Panel > Options: Level (Dataframe, Column or Cell,
  Column by default) says where the options below it are set, and Chart options > Show Legend turns
  the legend off or on. Set for the column, it changes every cell of the column; set for one cell, only
  that cell. Demo: the Assay Curves table of the Curves package, one named curve per cell.

  What each cell then draws is the renderer's, and the package's legend-tests claim it ("Show Legend
  set for the column names no curve, set for one cell names that one again"); the guide claims what
  the panel shows: another cell of the column keeps the column's setting, the one set for itself its own.

  Scenario: Hide the legend for the column, then show it in one cell
    Given user is logged in
    And simple mode is off
    And user opens assay-curves dataset
    And user removes "Molecule" column
    And user sets properties of grid:
      | Row Height        | 290 |
      | Col Header Height | 30  |
    And the "Fitted Curve" column of the grid is 340 px wide
    And the context panel is open
    When user clicks on the "cell 1 of Fitted Curve" area of grid
    And user expands "Options" pane in context panel
    Then "Level" input in context panel should have value "Column"
    When user scrolls "Show Legend" input in context panel into view
    And user unchecks "Show Legend" input in context panel
    And user scrolls "Level" input in context panel into view
    And user selects "Cell" in "Level" input in context panel
    And user scrolls "Show Legend" input in context panel into view
    And user checks "Show Legend" input in context panel
    And user clicks on the "cell 2 of Fitted Curve" area of grid
    Then "Show Legend" input in context panel should be unchecked
    When user clicks on the "cell 1 of Fitted Curve" area of grid
    Then "Show Legend" input in context panel should be checked
    And no errors should have been logged
