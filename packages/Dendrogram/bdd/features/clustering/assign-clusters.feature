@journey @realizes:dendrogram.cp.assign-clusters-column-creation
Feature: Assigning clusters from the tree next to the grid
  mol1K, whose molecule column the grid draws as structures, clustered on that column (euclidean,
  ward): the task bar shows "Creating dendrogram ...", and a tree with a leaf per row appears beside
  the grid. The tree follows the grid — the current row, the selection, the filter — and the grid
  follows the tree: a click on a leaf makes its row current, hovering it puts its row under the
  mouse, a click on an inner node selects its leaves. The wheel scrolls the tree with the grid;
  Control with the wheel zooms it between 1 and 100 times, and Reset Zoom puts it back.
  Assign Clusters opens from the tree's context menu and from the magic wand, draws its cut line at
  the threshold, and its Threshold and Clusters inputs drive each other. The tree is 639.09 high and
  the dialog opens at half of it, which cuts it into 6 clusters; a threshold of 200 gives 4 clusters
  and 5 clusters sit at 239.66. Assign adds a string column named after the threshold with a cluster
  number in every row, next to the Medoid Rank and Avg Distance to Cluster columns, and a second
  Assign adds its own columns beside the first. Running the clustering again warns and replaces the
  tree.

  A double-click on the tree's empty margin leaves the zoom as it was, where the case expects it to
  reset. That scenario is last and a known failure; the tag goes when the double-click resets the zoom.

  Background:
    Given user is logged in
    And user opens mol1K dataset
    And user watches the task bar
    When user picks "Chem > Analyze > Hierarchical Clustering..." from the top menu
    And user clicks on OK button in "Hierarchical Clustering" dialog
    Then the "tree leaves" reading of grid should be 1000

  Scenario: The grid draws structures and the clustering shows its progress
    Then the "cell type of molecule" reading of grid should be "Molecule"
    And the task bar should have shown "Creating dendrogram"
    And no errors should have been logged

  Scenario: The tree follows the grid's current row, selection and filter
    When user clicks on the "cell 3 of prID" area of grid
    Then row 3 should be current
    And the "tree current node" reading of grid should be "2"
    When user selects the first 10 rows
    Then the "tree selected leaves" reading of grid should be 10
    When user filters rows where "pIC50_HIV_Integrase" is between 6 and 7
    Then fewer than 1000 rows should pass the filter
    And the "tree leaves" reading of grid should be lower than before
    And the "tree leaves" and "rows shown" readings of grid should be the same
    When user resets the filter
    And user selects no rows
    Then the "tree leaves" reading of grid should be 1000
    And the "tree selected leaves" reading of grid should be 0
    And no errors should have been logged

  Scenario: The grid follows the tree: a leaf click, a hover, an inner node click
    When user clicks on the "cell 266 of prID" area of grid
    Then the "tree current node" reading of grid should be "265"
    When user clicks on the "leaf 14" area of grid
    Then row 15 should be current
    And the "tree current node" reading of grid should be "14"
    When user hovers over the "leaf 104" area of grid
    Then the "tree mouse over node" reading of grid should be "104"
    And row 105 should be under the mouse
    When user clicks on the "node 265-393" area of grid
    Then 41 rows should be selected
    And the "tree selected leaves" reading of grid should be 41
    When user selects no rows
    Then the "tree selected leaves" reading of grid should be 0
    And no errors should have been logged

  Scenario: The tree's context menu opens Assign Clusters at half the tree's height, with its cut line there
    When user picks "Assign Clusters" from the context menu of the "tree" area of grid
    Then "Assign Clusters" dialog should be visible
    And Threshold input in "Assign Clusters" dialog should have a value between 319.5 and 319.6
    And Clusters input in "Assign Clusters" dialog should have value "6"
    And "Medoid columns" input in "Assign Clusters" dialog should be checked
    And the "cut threshold" reading of grid should be 319.54
    And grid should have a "cut line" area
    When user clicks on CANCEL button in "Assign Clusters" dialog
    Then "Assign Clusters" dialog should be hidden
    And grid should not have a "cut line" area
    And no errors should have been logged

  Scenario: The magic wand opens the same dialog
    When user hovers over "Assign Clusters" icon
    Then tooltip should contain text "Assign Clusters"
    When user clicks on "Assign Clusters" icon
    Then "Assign Clusters" dialog should be visible
    And Threshold input in "Assign Clusters" dialog should have a value between 319.5 and 319.6
    When user clicks on CANCEL button in "Assign Clusters" dialog
    Then "Assign Clusters" dialog should be hidden
    And no errors should have been logged

  Scenario: Threshold sets Clusters and moves the cut line, Clusters sets Threshold, and Assign adds the columns of that cut
    When user clicks on "Assign Clusters" icon
    And user enters "200" into Threshold input in "Assign Clusters" dialog
    Then Clusters input in "Assign Clusters" dialog should have value "4"
    And the "cut threshold" reading of grid should be 200
    When user enters "5" into Clusters input in "Assign Clusters" dialog
    Then Threshold input in "Assign Clusters" dialog should have a value between 239.65 and 239.67
    And the "cut threshold" reading of grid should be 239.66
    When user clicks on Assign button in "Assign Clusters" dialog
    Then the "Assign Clusters" dialog should close
    And 3 new columns matching "\(239\.66\)$" should have been added
    And the table should have a column "Cluster (239.66)"
    And the table should have a column "Medoid Rank (239.66)"
    And the table should have a column "Avg Distance to Cluster (239.66)"
    And "Cluster (239.66)" column should have type "string"
    And "Cluster (239.66)" column should have no missing values
    And the newest column matching "^Cluster \(" should have 5 distinct values
    And "Medoid Rank (239.66)" column should have no missing values
    And "Avg Distance to Cluster (239.66)" column should have no missing values
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A second Assign adds its own columns and keeps the first
    When user clicks on "Assign Clusters" icon
    And user enters "3" into Clusters input in "Assign Clusters" dialog
    And user clicks on Assign button in "Assign Clusters" dialog
    Then the "Assign Clusters" dialog should close
    And 3 new columns matching "\(159\.77\)$" should have been added
    And 2 new columns matching "^Cluster \(" should have been added
    And the table should have a column "Cluster (239.66)"
    And the newest column matching "^Cluster \(" should have 3 distinct values
    And the newest column matching "^Cluster \(" should have no missing values
    And no errors should have been logged

  Scenario: The wheel scrolls the tree with the grid, and Control with the wheel zooms it between 1 and 100
    When user scrolls the mouse wheel down over the "tree" area of grid
    Then the "tree top row" reading of grid should be higher than before
    When user scrolls the mouse wheel up over the "tree" area of grid
    Then the "tree top row" reading of grid should be 0
    When user scrolls the mouse wheel up over the "tree" area of grid holding Control
    Then the "tree zoom" reading of grid should be higher than before
    When user scrolls the mouse wheel up 150 times over the "tree" area of grid holding Control
    Then the "tree zoom" reading of grid should be 100
    When user scrolls the mouse wheel down 150 times over the "tree" area of grid holding Control
    Then the "tree zoom" reading of grid should be 1
    And no errors should have been logged

  Scenario: Reset Zoom from the tree's context menu puts the zoom back
    When user scrolls the mouse wheel up 3 times over the "tree" area of grid holding Control
    Then the "tree zoom" reading of grid should be higher than before
    When user picks "Reset Zoom" from the context menu of the "tree" area of grid
    Then the "tree zoom" reading of grid should be 1
    And no errors should have been logged

  Scenario: Running the clustering again warns and replaces the tree
    When user picks "Chem > Analyze > Hierarchical Clustering..." from the top menu
    And user clicks on OK button in "Hierarchical Clustering" dialog
    Then a warning balloon containing "Closing existing dendrogram" should have been shown
    And the "tree leaves" reading of grid should be 1000
    And there should be 1 visible "Assign Clusters" icon
    And there should be 1 visible "Remove Dendrogram" icon
    And no errors should have been logged

  Scenario: Control with the wheel zooms the new tree in
    When user scrolls the mouse wheel up 3 times over the "tree" area of grid holding Control
    Then the "tree zoom" reading of grid should be higher than before

  @known-failure
  Scenario: A double-click on the tree's empty margin resets the zoom
    When user double-clicks on the "tree margin" area of grid
    Then the "tree zoom" reading of grid should be 1
