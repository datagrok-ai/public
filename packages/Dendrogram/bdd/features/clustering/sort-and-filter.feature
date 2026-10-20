@journey @realizes:GROK-13041
Feature: The tree next to the grid under a filter and under a sort
  A filter set in the filter panel on a table whose grid carries a tree hides rows without reordering
  them: the tree keeps a leaf for every row the grid shows and shows no "Revert columns sort order
  to see Dendrogram Tree" overlay (GROK-13041). A sort from the column header does show the overlay
  with a Revert sort button, and Revert sort puts the grid's rows back in the order they had.

  Background:
    Given user is logged in
    And user opens mol1K dataset
    And user watches the task bar
    When user picks "Chem > Analyze > Hierarchical Clustering..." from the top menu
    And user clicks on OK button in "Hierarchical Clustering" dialog
    Then the task bar should have finished "Creating dendrogram"
    And the "tree leaves" reading of grid should be 1000

  Scenario: A filter keeps the tree in step with the grid and shows no overlay
    When user remembers the "tree leaves" reading of grid
    And user clicks on first "Toggle filters" icon
    And user clicks on the "category Active_Integrase of Activity_Integrase" area of filter panel
    Then fewer than 1000 rows should pass the filter
    And the "tree leaves" reading of grid should not be as remembered
    And the "tree leaves" and "rows shown" readings of grid should be the same
    And grid should not contain text "Revert columns sort order to see Dendrogram Tree"
    And "Revert sort" button should be absent
    And "Assign Clusters" icon should be visible
    When user clicks on the "checkbox Inactive_Integrase of Activity_Integrase" area of filter panel
    Then all rows should pass the filter
    And the "tree leaves" reading of grid should be 1000
    And grid should not contain text "Revert columns sort order to see Dendrogram Tree"
    And no errors should have been logged

  Scenario: A sort shows the overlay, and Revert sort puts the rows back in the tree's order
    When user remembers the "row order" reading of grid
    And user double-clicks on the "header pIC50_HIV_Integrase" area of grid
    Then the "sort column" reading of grid should be "pIC50_HIV_Integrase"
    And the "row order" reading of grid should not be as remembered
    And "Revert sort" button should be visible
    And grid should contain text "Revert columns sort order to see Dendrogram Tree"
    When user clicks on "Revert sort" button
    Then the "row order" reading of grid should be as remembered
    And "Revert sort" button should be absent
    And grid should not contain text "Revert columns sort order to see Dendrogram Tree"
    And no errors should have been logged
