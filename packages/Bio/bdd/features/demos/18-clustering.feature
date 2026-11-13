@guide @demo
Feature: Bioinformatics · Clustering
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag
  leaves the set-up and the checks off camera). On FASTA_PT_activity (99 peptides of 16 residues in
  three families, with an activity): Bio > Analyze > Hierarchical Clustering (the Dendrogram package)
  clusters the sequences by their distance and attaches the tree to the grid, the rows in the tree's
  order; a click on a branch selects the peptides under it, and Assign Clusters cuts the tree into as
  many clusters as asked, adding a cluster column to the table (nothing is saved on the server). Off
  camera: the rows are set 7 px tall so the whole tree is in view, the activity is colored from blue
  (2.0 and below) to red (4.5), and after the cut the grid shows the new cluster column (its medoid
  columns hidden) colored by cluster, since a 7 px row draws no text.

  Scenario: Cluster the peptides into a tree and cut it
    A tree of the peptides by sequence distance grows beside the grid; a click on a branch selects its family, and a cut assigns clusters.

    Given user is logged in
    And simple mode is off
    And user opens FASTA_PT_activity dataset
    And the Bio package is initialized
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "sequence, activity"
    And the "sequence" column of the grid is 220 px wide
    And the rows of the grid are 7 px tall
    And user colors "activity" column linearly from "#4575b4" to "#d73027" over 2.0 to 4.5
    And user watches the task bar
    When user picks "Bio > Analyze > Hierarchical Clustering..." from the top menu
    Then "Hierarchical Clustering" dialog should be visible
    When user clicks on OK button in "Hierarchical Clustering" dialog
    Then the task bar should have finished "Creating dendrogram"
    And the "tree leaves" reading of grid should be 99
    When the video says "A tree of the peptides by sequence distance"
    And user clicks on the first "node" area of grid
    Then the "tree selected leaves" reading of grid should be at least 20
    When user clicks on "Assign Clusters" icon
    Then "Assign Clusters" dialog should be visible
    When user enters "3" into Clusters input in "Assign Clusters" dialog
    And user clicks on Assign button in "Assign Clusters" dialog
    Then the "Assign Clusters" dialog should close
    And the newest column matching "^Cluster \(" should have 3 distinct values
    And every value of "Cluster (6.60)" column should match "^[123]$"
    Given the grid shows only the columns "sequence, activity, Cluster (6.60)"
    And user colors "Cluster (6.60)" column categorically:
      | 1 | #8e6cc0 |
      | 2 | #2a9d8f |
      | 3 | #e9c46a |
    And the "Cluster (6.60)" column of the grid is 110 px wide
    And user scrolls the grid to the "sequence" column
    When the video says "Three clusters, three families: one stands out in activity"
    Then no error or warning balloon should have been shown
    And no errors should have been logged
