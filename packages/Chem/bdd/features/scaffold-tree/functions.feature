@journey @realizes:chem.cp.scaffold-tree-add-filter
Feature: Scaffold Tree — add, generate, filter, inspect
  The canonical walk on smiles-50 (50 rows, molecule column canonical_smiles). Chem | Analyze |
  Scaffold Tree adds the viewer to the table view in its empty state, reading "Scaffold Tree is
  empty". The magic wand, whose tooltip is "Generate from molecular column", builds the scaffold
  hierarchy. Checking the first scaffold node filters the table to the molecules that contain that
  scaffold and the viewer reports the rows it keeps. Its toolbar offers generate, sketch, upload,
  save, expand/collapse, clear filter and drop, and it reports the column and size it holds.
  Add Filter | Scaffold Tree Filter... offers the molecule column alone and puts a scaffold tree card
  on the filter panel.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens smiles-50 dataset

  Scenario: The dataset opens with a molecule column
    Then "canonical_smiles" column should have semantic type "Molecule"
    And the table should have 50 rows
    And no errors should have been logged

  Scenario: The menu adds the viewer in its empty state
    When user picks "Chem > Analyze > Scaffold Tree" from the top menu
    Then the top menu command should have completed
    And Scaffold Tree viewer should be visible
    And the open tableview should have 1 Scaffold Tree viewer
    And Scaffold Tree viewer should be bound to table "smiles-50"
    And the "molecule column" reading of Scaffold Tree viewer should be "canonical_smiles"
    And the "nodes" reading of Scaffold Tree viewer should be 0
    And the "message" reading of Scaffold Tree viewer should include the text "Scaffold Tree is empty"
    And no errors should have been logged

  Scenario: The magic wand generates the scaffold hierarchy
    When user hovers over Scaffold Tree viewer
    Then "Generate" icon inside Scaffold Tree viewer should be enabled
    When user clicks on "Generate" icon inside Scaffold Tree viewer
    Then Scaffold Tree viewer should have finished building its tree
    Then the "nodes" reading of Scaffold Tree viewer should be at least 4
    And the "root nodes" reading of Scaffold Tree viewer should be at least 1
    And the "message" reading of Scaffold Tree viewer should be ""
    And the "checked nodes" reading of Scaffold Tree viewer should be 0
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: A click on the first scaffold node filters the table
    When user clicks on the "checkbox of node 1" area of Scaffold Tree viewer
    Then the "checked nodes" reading of Scaffold Tree viewer should be 1
    And the "checked of node 1" reading of Scaffold Tree viewer should be "true"
    And fewer than 50 rows should pass the filter
    And the "hits of node 1" reading of Scaffold Tree viewer should be at least 1
    And the "rows kept" reading of Scaffold Tree viewer should be at least 1
    And no errors should have been logged

  Scenario: The viewer's own toolbar offers its actions
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

  Scenario: The viewer reports the column and the size it was built with
    Then "size" property of Scaffold Tree viewer should be "large"
    And "molecule" property of Scaffold Tree viewer should be "canonical_smiles"
    And no errors should have been logged

  Scenario: Add Filter | Scaffold Tree Filter... offers the molecule column alone
    When user clicks on filter icon in toolbar
    And user picks "Add Filter | Scaffold Tree Filter..." from the viewer menu of filter panel
    Then "Select columns..." dialog should be visible
    And the "text of cell 1 of __name" reading of grid viewer in "Select columns..." dialog should be "canonical_smiles"
    When user clicks on All label in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then there should be 2 visible "canonical_smiles" filter card
    And all rows should pass the filter
    And no errors should have been logged
