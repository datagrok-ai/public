@journey @realizes:dendrogram.cp.hier-clustering-chem-dialog-end-to-end
Feature: Hierarchical clustering from the Chem menu
  Chem | Analyze | Hierarchical Clustering... over mol1K opens the clustering dialog on the molecule
  column, with the two distances and the seven linkages; OK attaches a tree to the grid with a leaf
  for every row. The Remove Dendrogram icon takes the tree away, and a run with another distance and
  linkage attaches a tree of another height. The task bar shows "Creating dendrogram ..." while the
  tree is built.

  Centroid linkage on the numeric columns attaches no tree: turning the cluster matrix into a tree
  fails with "Cannot read properties of undefined (reading 'children')" (GROK-19595). The last
  scenario states the tree it should attach and is a known failure; the tag goes when the tree has
  a leaf for every row.

  Background:
    Given user is logged in
    And user opens mol1K dataset

  Scenario: The dialog opens on the molecule column with every distance and linkage
    When user picks "Chem > Analyze > Hierarchical Clustering..." from the top menu
    Then "Hierarchical Clustering" dialog should be visible
    And Table input in "Hierarchical Clustering" dialog should have value "mol1K"
    And editor of Features input in "Hierarchical Clustering" dialog should have text "(1) molecule"
    And Distance input in "Hierarchical Clustering" dialog should have value "euclidean"
    And Linkage input in "Hierarchical Clustering" dialog should have value "ward"
    And Distance input in "Hierarchical Clustering" dialog should offer "euclidean, manhattan"
    And Linkage input in "Hierarchical Clustering" dialog should offer "single, complete, average, weighted, centroid, median, ward"
    When user clicks on CANCEL button in "Hierarchical Clustering" dialog
    Then "Hierarchical Clustering" dialog should be hidden
    And no errors should have been logged

  Scenario: Euclidean distance and ward linkage attach a tree with a leaf for every molecule
    Given user watches the task bar
    When user picks "Chem > Analyze > Hierarchical Clustering..." from the top menu
    And user clicks on OK button in "Hierarchical Clustering" dialog
    Then "Hierarchical Clustering" dialog should be hidden
    And the "tree leaves" reading of grid should be 1000
    And the task bar should have shown "Creating dendrogram"
    And "Assign Clusters" icon should be visible
    And no errors should have been logged

  Scenario: Removing the tree, then manhattan distance and single linkage attach a tree of another height
    When user remembers the "tree height" reading of grid
    And user clicks on "Remove Dendrogram" icon
    Then "Assign Clusters" icon should be absent
    And grid should not report a "tree leaves" reading
    When user picks "Chem > Analyze > Hierarchical Clustering..." from the top menu
    And user selects "manhattan" in Distance input in "Hierarchical Clustering" dialog
    And user selects "single" in Linkage input in "Hierarchical Clustering" dialog
    Then Distance input in "Hierarchical Clustering" dialog should have value "manhattan"
    And Linkage input in "Hierarchical Clustering" dialog should have value "single"
    When user clicks on OK button in "Hierarchical Clustering" dialog
    Then the "tree leaves" reading of grid should be 1000
    And the "tree height" reading of grid should not be as remembered
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Numeric columns with median linkage attach a tree with a leaf for every row
    When user remembers the "tree height" reading of grid
    And user clicks on "Remove Dendrogram" icon
    And user picks "Chem > Analyze > Hierarchical Clustering..." from the top menu
    And user clicks on editor of Features input in "Hierarchical Clustering" dialog
    Then "Select columns..." dialog should be visible
    When user clicks on None label in "Select columns..." dialog
    And user clicks on the "cell 4 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 5 of x" area of grid viewer in "Select columns..." dialog
    Then the "text of cell 4 of __name" reading of grid viewer in "Select columns..." dialog should be "pIC50_HIV_Integrase"
    And the "text of cell 5 of __name" reading of grid viewer in "Select columns..." dialog should be "Q"
    And "Select columns..." dialog should contain text "2 checked"
    When user clicks on OK button in "Select columns..." dialog
    Then editor of Features input in "Hierarchical Clustering" dialog should contain text "(2)"
    When user selects "median" in Linkage input in "Hierarchical Clustering" dialog
    Then Linkage input in "Hierarchical Clustering" dialog should have value "median"
    When user clicks on OK button in "Hierarchical Clustering" dialog
    Then the "tree leaves" reading of grid should be 1000
    And the "tree height" reading of grid should not be as remembered
    And no errors should have been logged

  Scenario: Numeric columns with centroid linkage are run
    When user clicks on "Remove Dendrogram" icon
    And user picks "Chem > Analyze > Hierarchical Clustering..." from the top menu
    And user clicks on editor of Features input in "Hierarchical Clustering" dialog
    And user clicks on None label in "Select columns..." dialog
    And user clicks on the "cell 4 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 5 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user selects "centroid" in Linkage input in "Hierarchical Clustering" dialog
    Then Linkage input in "Hierarchical Clustering" dialog should have value "centroid"
    When user clicks on OK button in "Hierarchical Clustering" dialog
    Then "Hierarchical Clustering" dialog should be hidden

  @known-failure @GROK-19595
  Scenario: Numeric columns with centroid linkage attach a tree with a leaf for every row
    Then the "tree leaves" reading of grid should be 1000
    And no errors should have been logged
