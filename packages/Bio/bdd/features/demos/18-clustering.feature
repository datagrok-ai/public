@guide @demo
Feature: 4 · A peptide series · Clustering
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the
  set-up and the checks off camera). On FASTA_PT_activity (99 peptides of 16 residues in three families, with
  an activity): Bio > Analyze > Hierarchical Clustering (the Dendrogram package) clusters the sequences by
  their distance and attaches the tree to the grid, the rows in the tree's order; a click on a branch selects
  the peptides under it, and Assign Clusters cuts the tree into clusters, adding a cluster column to the table
  (nothing is saved on the server). In the Assign Clusters dialog, opened from the tree's icon, the threshold
  is typed (6.5, three clusters; the dialog opens on the tree's half height written to sixteen decimals,
  replaced as the typing starts) and Medoid columns is unchecked, so the cut adds the cluster column alone.
  After the cut, off camera, the grid shows the new column colored by cluster, and a box plot of the activity
  per cluster moves to the right of the grid. The rows are 14 px tall, so half of the tree is in view with
  every residue readable; the activity is colored from blue (2.0 and below) to red (4.5).

  Scenario: Cluster the peptides into a tree and cut it
    The peptides as a tree by sequence, cut into clusters compared by activity.

    Given user is logged in
    And simple mode is off
    And user opens FASTA_PT_activity dataset
    And the Bio package is initialized
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "sequence, activity"
    And the "sequence" column of the grid is 270 px wide
    And the rows of the grid are 14 px tall
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
    When the video says "Cut it where it splits into three"
    And user enters "6.5" into Threshold input in "Assign Clusters" dialog
    Then Clusters input in "Assign Clusters" dialog should have value "3"
    When user unchecks "Medoid columns" checkbox in "Assign Clusters" dialog
    And user clicks on Assign button in "Assign Clusters" dialog
    Then the "Assign Clusters" dialog should close
    And the newest column matching "^Cluster \(" should have 3 distinct values
    And every value of "Cluster (6.50)" column should match "^[123]$"
    And the table should not have a column "Medoid Rank (6.50)"
    Given the grid shows only the columns "sequence, activity, Cluster (6.50)"
    And user colors "Cluster (6.50)" column categorically:
      | 1 | #8e6cc0 |
      | 2 | #2a9d8f |
      | 3 | #e9c46a |
    And the "Cluster (6.50)" column of the grid is 125 px wide
    And user adds a box plot viewer with:
      | valueColumnName     | activity       |
      | category1ColumnName | Cluster (6.50) |
    And box plot viewer is docked to the right of the grid, 45% of the view wide
    And user scrolls the grid to the "sequence" column
    Then box plot viewer should be painted
    When the video says "Three clusters, three families: one stands out in activity"
    Then no error or warning balloon should have been shown
    And no errors should have been logged
