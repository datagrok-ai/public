@guide @tutorial
Feature: How-to · ML
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). ML > Cluster > Hierarchical Clustering... is the Dendrogram
  package's dialog; for numeric features the distance matrix and the clustering run in browser workers
  (DistanceMatrixService, getClusterMatrixWorker), no script. On iris the Features take All and leave out the row
  number "col 1"; Euclidean distance, Ward linkage; OK draws the tree beside the grid, a leaf per row. Assign Clusters
  from the tree's context menu cuts it into 3 clusters and adds the cluster column (with Medoid Rank and Avg Distance
  to Cluster); with 3 clusters the cut is at 38.58, which names the columns. Off camera, Species is colored by
  category and the grid's rows are 6 px tall, so all 150 rows and the whole tree fit on screen and the species show
  as colored bands beside the branches (at the default height only 28 setosa rows and a sliver of the tree show).
  Text is not drawn in 6 px rows, so the clip colors the new Cluster column on camera (Color Coding > Categorical
  from its header's menu: the default palette, 1 blue, 2 light orange, 3 green; the pick is a radio item that leaves
  the menu open, so the clip closes it). Species gets the same three colors off
  camera, matched to the clusters Ward linkage gives (setosa is cluster 1, cluster 2 is virginica only, cluster 3 is
  versicolor and the rest of virginica), so where the colors agree the tree found the species. Colored off camera,
  the Cluster column was a jump: it came out of Assign already painted, which Assign does not do.
  iris is opened off camera; nothing is made on the server. k-means (Cluster Data) is in the "Work with data"
  playlist.

  Scenario: Cluster rows into a tree
    Hierarchical clustering joins similar rows step by step into a tree; cut the tree to get clusters.

    Given user is logged in
    And simple mode is off
    And user opens iris dataset
    And user colors "Species" column categorically:
      | setosa     | #1f77b4 |
      | versicolor | #2ca02c |
      | virginica  | #ffbb78 |
    And the rows of the grid are 6 px tall
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When the video says "Group the flowers by their four measurements" at the bottom
    # caption: Hierarchical Clustering is under ML > Cluster
    And user picks "ML > Cluster > Hierarchical Clustering..." from the top menu
    Then "Hierarchical Clustering" dialog should be visible
    # caption: Features: the columns that make rows alike
    When user clicks on Features input in "Hierarchical Clustering" dialog
    Then "Select columns..." dialog should be visible
    When user clicks on "All" link in "Select columns..." dialog
    # caption: Leave out the row number
    And user toggles the "col 1" column in the column list of "Select columns..." dialog
    Then "Select columns..." dialog should contain text "4 checked"
    When user clicks on OK button in "Select columns..." dialog
    # caption: Distance and Linkage: how closeness is measured and groups are joined
    And user hovers over Linkage input in "Hierarchical Clustering" dialog
    Given user watches the task bar
    When user clicks on OK button in "Hierarchical Clustering" dialog
    Then the task bar should have finished "Creating dendrogram"
    And the "tree leaves" reading of grid should be 150
    When the video says "The tree beside the grid: one leaf per row, branches follow the species" at the bottom
    And user hovers over the "tree" area of grid
    # caption: Right-click the tree to cut it into clusters
    And user picks "Assign Clusters" from the context menu of the "tree" area of grid
    Then "Assign Clusters" dialog should be visible
    # caption: Clusters: how many groups to cut the tree into
    When user enters "3" into Clusters input in "Assign Clusters" dialog
    Then the "cut threshold" reading of grid should be a finite number
    When user clicks on Assign button in "Assign Clusters" dialog
    Then the "Assign Clusters" dialog should close
    And a new column matching "^Cluster \(" should have been added
    And the newest column matching "^Cluster \(" should have 3 distinct values
    # caption: Assign adds a Cluster column; color it to see the groups
    When user picks "Color Coding > Categorical" from the context menu of the "header Cluster (38.58)" area of grid
    Then "Cluster (38.58)" column should be color-coded categorically
    When user closes the context menu
    # caption: Setosa gets a cluster of its own; the other two species overlap
    And user moves the pointer away from grid
    Then no error or warning balloon should have been shown
    And no errors should have been logged
