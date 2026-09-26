@journey @realizes:chem.cp.scaffold-tree-add-filter
Feature: Scaffold Tree — added from the menu, and its filter card
  On smiles-50 (50 rows, one molecule column canonical_smiles), Chem | Analyze | Scaffold Tree adds
  the viewer to the table view in its empty state, reading "Scaffold Tree is empty", bound to that
  column, and it reports the column and the size it holds. Add Filter | Scaffold Tree Filter... offers
  the molecule column alone and puts a scaffold tree card on the filter panel. Generating a tree,
  checking a node, the viewer's toolbar and its Clear filter are claimed on spgi-100 in
  scaffold-tree.feature: a second generation here would claim nothing new.

  Background:
    Given user is logged in
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And user opens smiles-50 dataset

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

  Scenario: The viewer reports the column and the size it holds
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
