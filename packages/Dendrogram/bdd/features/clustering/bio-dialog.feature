@journey @realizes:dendrogram.cp.hier-clustering-bio-sequence-path
Feature: Hierarchical clustering from the Bio menu
  Bio | Analyze | Hierarchical Clustering... over FASTA_PT_activity opens the clustering dialog on
  the sequence column, a macromolecule the grid draws as sequences, with the two distances and the
  seven linkages. OK shows "Creating dendrogram ..." in the task bar and attaches a tree to the grid
  with a leaf for every sequence; the tree follows the grid's current row and puts a hovered leaf's
  row under the mouse; another distance and linkage attach a tree of another height; Assign
  Clusters works on the sequence tree, where 5 clusters sit at a threshold of 11.06.

  A threshold above the tree's height leaves Clusters at 0, below the input's minimum of 1. That
  scenario is last and a known failure; the tag goes when Clusters stays at 1 or more.

  Background:
    Given user is logged in
    And user opens FASTA_PT_activity dataset

  Scenario: The dialog opens on the sequence column with every distance and linkage
    Then "sequence" column should have semantic type "Macromolecule"
    And the "cell type of sequence" reading of grid should be "sequence"
    When user picks "Bio > Analyze > Hierarchical Clustering..." from the top menu
    Then "Hierarchical Clustering" dialog should be visible
    And Table input in "Hierarchical Clustering" dialog should have value "FASTA_PT_activity"
    And editor of Features input in "Hierarchical Clustering" dialog should have text "(1) sequence"
    And Distance input in "Hierarchical Clustering" dialog should have value "euclidean"
    And Linkage input in "Hierarchical Clustering" dialog should have value "ward"
    And Distance input in "Hierarchical Clustering" dialog should offer "euclidean, manhattan"
    And Linkage input in "Hierarchical Clustering" dialog should offer "single, complete, average, weighted, centroid, median, ward"
    When user clicks on CANCEL button in "Hierarchical Clustering" dialog
    Then "Hierarchical Clustering" dialog should be hidden
    And no errors should have been logged

  Scenario: Euclidean distance and ward linkage attach a tree with a leaf for every sequence
    Given user watches the task bar
    When user picks "Bio > Analyze > Hierarchical Clustering..." from the top menu
    And user clicks on OK button in "Hierarchical Clustering" dialog
    Then "Hierarchical Clustering" dialog should be hidden
    And the "tree leaves" reading of grid should be 99
    And the task bar should have shown "Creating dendrogram"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The sequence tree follows the grid, and the grid follows the tree
    When user clicks on the "cell 3 of sequence_id" area of grid
    Then row 3 should be current
    And the "tree current node" reading of grid should be "2"
    When user hovers over the "leaf 10" area of grid
    Then the "tree mouse over node" reading of grid should be "10"
    And row 11 should be under the mouse
    And no errors should have been logged

  Scenario: Removing the tree, then manhattan distance and complete linkage attach a tree of another height
    When user remembers the "tree height" reading of grid
    And user clicks on "Remove Dendrogram" icon
    Then grid should not report a "tree leaves" reading
    When user picks "Bio > Analyze > Hierarchical Clustering..." from the top menu
    And user selects "manhattan" in Distance input in "Hierarchical Clustering" dialog
    And user selects "complete" in Linkage input in "Hierarchical Clustering" dialog
    Then Distance input in "Hierarchical Clustering" dialog should have value "manhattan"
    And Linkage input in "Hierarchical Clustering" dialog should have value "complete"
    When user clicks on OK button in "Hierarchical Clustering" dialog
    Then the "tree leaves" reading of grid should be 99
    And the "tree height" reading of grid should not be as remembered
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Assign Clusters on the sequence tree adds a cluster number to every sequence
    When user clicks on "Remove Dendrogram" icon
    And user picks "Bio > Analyze > Hierarchical Clustering..." from the top menu
    And user clicks on OK button in "Hierarchical Clustering" dialog
    Then the "tree leaves" reading of grid should be 99
    When user clicks on "Assign Clusters" icon
    Then "Assign Clusters" dialog should be visible
    When user enters "5" into Clusters input in "Assign Clusters" dialog
    Then Threshold input in "Assign Clusters" dialog should have a value between 11.05 and 11.07
    When user clicks on Assign button in "Assign Clusters" dialog
    Then the "Assign Clusters" dialog should close
    And the table should have a column "Cluster (11.06)"
    And "Cluster (11.06)" column should have no missing values
    And the newest column matching "^Cluster \(" should have 5 distinct values
    And no errors should have been logged

  @known-failure
  Scenario: A threshold above the tree's height keeps at least one cluster
    When user clicks on "Assign Clusters" icon
    And user enters "100" into Threshold input in "Assign Clusters" dialog
    Then Clusters input in "Assign Clusters" dialog should have a value between 1 and 99
