@journey @realizes:bio.menu.manage.monomers @realizes:bio.menu.manage.match-with-library @realizes:bio.op.standardise_monomer_library
Feature: Monomers, and matching molecules against the library
  Bio | Manage | Monomers opens a table view of every monomer the selected libraries hold, one row
  per monomer with its symbol, name and polymer type. Match with Monomer Library asks for a
  molecule column and the polymer type to match against — peptide, RNA or chemical monomers. The
  library standardization behind it keeps every monomer's symbol and polymer type. The Manage
  Monomer Libraries app is the same view as the menu command.

  Not translated: the Manage Monomer Libraries app's tree browser in Browse (its children load
  with no signal the tree gives, a core gap); creating or editing a monomer in the Monomers view
  (its editor dialog has no feature yet); running Match on a molecule column — the features' Bio
  datasets hold no Molecule column, and Chem's is Chem's to test.

  Background:
    Given user is logged in
    And user opens filter_HELM dataset
    And the Bio package is initialized

  Scenario: Bio > Manage > Monomers lists every monomer of the libraries
    When user picks "Bio > Manage > Monomers" from the top menu
    Then the "Manage Monomers" view should be current
    And the current view should hold at least 1 viewer
    And the table should have a column "Symbol"
    And the table should have a column "Polymer Type"
    And "Symbol" column should have no missing values
    And "Polymer Type" column should have no missing values
    And the monomer sketcher of the Manage Monomers view should be ready
    When user closes the current view

  Scenario: The Manage Monomer Libraries app opens the manager view
    Given user opens the "Manage Monomer Libraries" app
    Then the "Manage Monomer Libraries" view should be current
    And "HELMCoreLibrary.json" checkbox should be visible
    When user closes the current view

  Scenario: Match with Monomer Library offers the three polymer types
    When user switches to the "filter_HELM" table view
    And user picks "Bio > Manage > Match with Monomer Library..." from the top menu
    Then "Match with Monomer Library" dialog should be visible
    And "Polymer Type" input in "Match with Monomer Library" dialog should offer "PEPTIDE, RNA, CHEM"
    And "Polymer Type" input in "Match with Monomer Library" dialog should have value "PEPTIDE"
    When user clicks on CANCEL button in "Match with Monomer Library" dialog
    Then "Match with Monomer Library" dialog should be hidden

  Scenario: Standardizing the core library keeps its monomers and their polymer types
    When user standardises the "HELMCoreLibrary.json" monomer library
    Then the result should be a list of 500 or more items
    And the result should hold "PEPTIDE" monomer "A"
    And the result should hold "PEPTIDE" monomer "meI"
    And the result should hold "RNA" monomer "A"
    And no errors should have been logged
