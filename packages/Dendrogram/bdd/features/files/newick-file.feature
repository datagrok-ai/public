@journey @realizes:dendrogram.cp.newick-file-open-via-files-browser
Feature: A newick file opened from Browse and through its file handler
  nwk1.nwk ships with the package: leaf1, and an inner node holding leaf2 and leaf3. Opened from the
  Files tree — with a click or with a double-click on the file — it shows as a Phylocanvas GL viewer
  that fills its view, over the table the file parses into, a row per node. The handler registered
  for .nwk makes a table view of that table instead — the file's text on its .newick tag — with a
  Dendrogram viewer drawing the same leaves; Browse opens the preview, so the handler is called
  directly.

  Background:
    Given user is logged in
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And "Files > App Data" tree node inside browse tree is expanded
    And "Files > App Data > Dendrogram" tree node inside browse tree is expanded

  Scenario: A double-click on the file shows its tree
    When user clicks on "Files > App Data > Dendrogram > data" tree node inside browse tree
    And user double-clicks on nwk1.nwk link in gallery
    Then PhylocanvasGL viewer should be visible
    And PhylocanvasGL viewer should fill its parent
    And the table of PhylocanvasGL viewer should hold a tree with leaves "leaf1, leaf2, leaf3"
    And no errors should have been logged

  Scenario: A click on the file previews the same tree
    When user clicks on "Files > App Data > Dendrogram > data" tree node inside browse tree
    And user clicks on nwk1.nwk link in gallery
    Then PhylocanvasGL viewer should be visible
    And PhylocanvasGL viewer should fill its parent
    And the table of PhylocanvasGL viewer should hold a tree with leaves "leaf1, leaf2, leaf3"
    And no errors should have been logged

  Scenario: The file handler makes a tree table with a Dendrogram of the same leaves
    When user opens the newick file "System:AppData/Dendrogram/data/nwk1.nwk" with its file handler
    Then Dendrogram viewer should be visible
    And table "Table" should have 5 rows
    And table "Table" should have columns "node, parent, leaf, distance"
    And the table should have tag ".newick" equal to the text of "System:AppData/Dendrogram/data/nwk1.nwk" file
    And the "leaves" reading of Dendrogram viewer should be "leaf1, leaf2, leaf3"
    And no errors should have been logged
