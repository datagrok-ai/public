@journey @realizes:chem.cp.scaffold-tree-add-filter
Feature: The Scaffold Tree viewer — building, checking, editing and filtering
  Chem | Analyze | Scaffold Tree adds the viewer to spgi-100 (100 rows, molecule column Structure) in
  its empty state. The magic wand builds a tree from Structure. Checking its first node keeps exactly
  the molecules of Structure that contain that node's own scaffold, and unchecking it lets every row
  through again; the toolbar offers generate, sketch, upload, save, expand/collapse, clear filter and
  drop, and Clear filter unchecks a checked node. Edit scaffold on the node replaces its structure
  with quinoline and the filter follows it. A scaffold sketched by hand on the plus icon becomes a node that keeps the molecules containing it, and a
  clone of the view carries the same scaffold and the same filtered rows.

  Background:
    Given user is logged in
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And user opens spgi dataset

  Scenario: The menu adds the viewer in its empty state
    When user picks "Chem > Analyze > Scaffold Tree" from the top menu
    Then the top menu command should have completed
    And Scaffold Tree viewer should be visible
    And Scaffold Tree viewer should be bound to table "spgi-100"
    And the "nodes" reading of Scaffold Tree viewer should be 0
    And the "message" reading of Scaffold Tree viewer should include the text "Scaffold Tree is empty"
    And the "generate blocked reason" reading of Scaffold Tree viewer should be ""
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: The magic wand builds a tree from the molecule column
    When user hovers over Scaffold Tree viewer
    And user clicks on "Generate" icon inside Scaffold Tree viewer
    Then Scaffold Tree viewer should have finished building its tree
    And the "nodes" reading of Scaffold Tree viewer should be at least 4
    And the "root nodes" reading of Scaffold Tree viewer should be at least 1
    And the "checked nodes" reading of Scaffold Tree viewer should be 0
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: Checking a node keeps the molecules that contain its scaffold
    When user clicks on the "checkbox of node 1" area of Scaffold Tree viewer
    Then the "checked nodes" reading of Scaffold Tree viewer should be 1
    And the "checked of node 1" reading of Scaffold Tree viewer should be "true"
    And fewer than 100 rows should pass the filter
    And the "hits of node 1" reading of Scaffold Tree viewer should be at least 1
    And the filter should pass exactly the molecules of "Structure" column containing the "scaffold of node 1" reading of Scaffold Tree viewer
    When user clicks on the "checkbox of node 1" area of Scaffold Tree viewer
    Then the "checked nodes" reading of Scaffold Tree viewer should be 0
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: The viewer's own toolbar offers its actions, and Clear filter unchecks the tree
    When user clicks on the "checkbox of node 1" area of Scaffold Tree viewer
    Then the "checked nodes" reading of Scaffold Tree viewer should be 1
    And the filter should pass exactly the molecules of "Structure" column containing the "scaffold of node 1" reading of Scaffold Tree viewer
    When user hovers over Scaffold Tree viewer
    Then the following elements should be visible:
      | "Generate" icon inside Scaffold Tree viewer                   |
      | "Sketch scaffolds manually" icon inside Scaffold Tree viewer  |
      | "Upload saved tree file" icon inside Scaffold Tree viewer     |
      | "Save this tree to disk" icon inside Scaffold Tree viewer     |
      | "Expand / collapse all" icon inside Scaffold Tree viewer      |
      | "Clear filter" icon inside Scaffold Tree viewer               |
      | "Drop all trees" icon inside Scaffold Tree viewer             |
    When user clicks on "Clear filter" icon inside Scaffold Tree viewer
    Then the "checked nodes" reading of Scaffold Tree viewer should be 0
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: Edit scaffold replaces the structure and the filter follows
    When user hovers over the "node 1" area of Scaffold Tree viewer
    And user clicks on the "edit icon of node 1" area of Scaffold Tree viewer
    Then sketcher dialog should be visible
    When user clears molecule input of sketcher dialog
    And user types "c1ccc2ncccc2c1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on "Save" button in sketcher dialog
    Then the "scaffold of node 1" reading of Scaffold Tree viewer should be the molecule "c1ccc2ncccc2c1"
    When user clicks on the "checkbox of node 1" area of Scaffold Tree viewer
    Then the "checked nodes" reading of Scaffold Tree viewer should be 1
    And the filter should pass exactly the molecules of "Structure" column containing "c1ccc2ncccc2c1"
    And no errors should have been logged


  Scenario: A scaffold sketched by hand keeps the molecules that contain it
    When user hovers over Scaffold Tree viewer
    And user clicks on "Drop all trees" icon inside Scaffold Tree viewer
    And user clicks on "Yes" button in "Delete Tree" dialog
    Then the "nodes" reading of Scaffold Tree viewer should be 0
    When user hovers over Scaffold Tree viewer
    And user clicks on "Sketch scaffolds manually" icon inside Scaffold Tree viewer
    And user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on "Add" button in sketcher dialog
    Then the "nodes" reading of Scaffold Tree viewer should be 1
    And the "scaffold of node 1" reading of Scaffold Tree viewer should be the molecule "c1ccncc1"
    And the "hits of node 1" reading of Scaffold Tree viewer should be 17
    When user clicks on the "checkbox of node 1" area of Scaffold Tree viewer
    Then 17 rows should pass the filter
    And the filter should pass exactly the molecules of "Structure" column containing "c1ccncc1"

  Scenario: Cloning the view brings the same scaffold and the same filtered rows
    When user picks "View > Layout > Clone View" from the top menu
    Then the open tableview should have 1 Scaffold Tree viewer
    And the "nodes" reading of Scaffold Tree viewer should be 1
    And the "scaffold of node 1" reading of Scaffold Tree viewer should be the molecule "c1ccncc1"
    And 17 rows should pass the filter
    And no errors should have been logged
