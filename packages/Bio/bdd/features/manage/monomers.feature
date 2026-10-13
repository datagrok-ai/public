@journey @realizes:bio.menu.manage.monomers @realizes:bio.menu.manage.match-with-library
Feature: Monomers, and matching molecules against the library
  Bio | Manage | Monomers opens a table view of every monomer the selected libraries hold, one row
  per monomer with its symbol, name and polymer type. Match with Monomer Library asks for a
  molecule column and the polymer type to match against — peptide, RNA or chemical monomers. The
  Manage Monomer Libraries app is the same view as the menu command. Match run on a column of
  amino acids drawn as molecules (a table the feature makes; Chem reads the molecules) opens a copy
  of the table that names the monomer of each: alanine is A, glycine G. Which molecules match which
  monomers, and how, is the matcher's to test (match-molecules.ts), not this feature's.

  Not translated: the Manage Monomer Libraries app's tree browser in Browse (its children load
  with no signal the tree gives, a core gap); creating and deleting a monomer in the Monomers
  view are monomer-editor.feature's, editing one is not written; the library standardization
  function called directly (no UI) — see the bdd library's CLAUDE.md, "What never becomes a
  feature".

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

  # last: the gate skips the rest of the journey on a stand without Chem
  Scenario: Match with Monomer Library names the monomer of each molecule
    Given the "Chem" package is installed
    And user opens a table "amino-acids" with:
      | molecule        |
      | C[C@H](N)C(=O)O |
      | NCC(=O)O        |
    And the semantic types of the current table are detected
    Then "molecule" column should have semantic type "Molecule"
    When user picks "Bio > Manage > Match with Monomer Library..." from the top menu
    Then "Match with Monomer Library" dialog should be visible
    When user clicks on OK button in "Match with Monomer Library" dialog
    Then the top menu command should have completed
    And the table should have a column "Matched monomer symbol"
    And the value of "Matched monomer symbol" column in row 1 should be "A"
    And the value of "Matched monomer symbol" column in row 2 should be "G"
    And no error or warning balloon should have been shown
    And no errors should have been logged
