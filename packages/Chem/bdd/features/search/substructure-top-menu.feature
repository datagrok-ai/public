@journey @realizes:chem.cp.substructure-search-top-menu
Feature: Substructure search from the Search menu
  Chem | Search | Substructure Search... on smiles puts an empty substructure card for
  canonical_smiles on the filter panel and opens its sketcher; nothing is filtered yet. Benzene
  keeps exactly the 924 molecules that contain it, a gold atom keeps none, and the table keeps its
  1000 rows. A second invocation asks for the molecule column (the first search leaves a hidden
  canonical SMILES column behind), then empties the card again; carboxylic acid keeps exactly the
  314 molecules that contain it, on the same single card.

  Background:
    Given user is logged in
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And user opens smiles dataset

  Scenario: The command prepares an empty card and opens its sketcher
    When user picks "Chem > Search > Substructure Search..." from the top menu
    Then sketcher dialog should be visible
    And the "cards" reading of filter panel should contain "canonical_smiles"
    And the "structure of canonical_smiles" reading of filter panel should be ""
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: Benzene keeps the molecules that contain it
    When user types "c1ccccc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then 924 rows should pass the filter
    And the "structure of canonical_smiles" reading of filter panel should be "c1ccccc1"
    And the filter should pass exactly the molecules of "canonical_smiles" column containing "c1ccccc1"
    And no errors should have been logged

  Scenario: A gold atom keeps no molecule and deletes no row
    When user clicks on the "card canonical_smiles" area of filter panel
    Then sketcher dialog should be visible
    When user types "[Au]" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then 0 rows should pass the filter
    And the filter should pass exactly the molecules of "canonical_smiles" column containing "[Au]"
    And the table should have 1000 rows
    And no errors should have been logged

  Scenario: A second invocation asks for the column and starts from an empty card
    When user picks "Chem > Search > Substructure Search..." from the top menu
    Then "Substructure search" dialog should be visible
    When user clicks on OK button in "Substructure search" dialog
    Then sketcher dialog should be visible
    And all rows should pass the filter
    When user types "C(=O)O" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then 314 rows should pass the filter
    And the filter should pass exactly the molecules of "canonical_smiles" column containing "C(=O)O"
    And the "filters" reading of filter panel should be 1
    And no errors should have been logged
