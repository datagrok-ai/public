@journey @realizes:chem.cp.names-to-smiles
Feature: Names To Smiles over a column of compound names
  Chem | Transform | Names To Smiles... looks the names up in ChEMBL and writes the structures into a
  canonical_smiles column beside them — aspirin and caffeine, compared with their structures through
  RDKit — so the lookup needs the ChEMBL database of the stand. On a
  table that already carries a canonical_smiles column the run ends in "Column named
  'canonical_smiles' already exists" and writes nothing (GROK-20955); that scenario is last and a
  known failure.

  Background:
    Given user is logged in
    And the package autostarts have completed

  Scenario: The names become molecules in a canonical_smiles column
    Given user opens a table "compound_names" with:
      | name     | mol |
      | aspirin  | CCO |
      | caffeine | CCC |
    And user sets the semantic type of "mol" column to "Molecule"
    When user picks "Chem > Transform > Names To Smiles..." from the top menu
    Then "Names To Smiles" dialog should be visible
    When user clicks on OK button in "Names To Smiles" dialog
    Then the top menu command should have completed
    And a new column "canonical_smiles" should have been added
    And "canonical_smiles" column should have semantic type "Molecule"
    And the molecule in row 1 of "canonical_smiles" column should be "CC(=O)Oc1ccccc1C(=O)O"
    And the molecule in row 2 of "canonical_smiles" column should be "Cn1c(=O)c2c(ncn2C)n(C)c1=O"
    And the table should have 2 rows
    And no errors should have been logged

  @known-failure @GROK-20955
  Scenario: The names go into a column of their own beside an existing canonical_smiles
    Given user opens a table "named_molecules" with:
      | name     | canonical_smiles |
      | aspirin  | CCO              |
      | caffeine | CCC              |
    And user sets the semantic type of "canonical_smiles" column to "Molecule"
    When user picks "Chem > Transform > Names To Smiles..." from the top menu
    And user clicks on OK button in "Names To Smiles" dialog
    Then the top menu command should have completed
    And 1 new column should have been added
    And no error or warning balloon should have been shown
