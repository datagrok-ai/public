@journey @realizes:chem.cp.elemental-analysis
Feature: Chem > Analyze commands that append columns to SMILES molecules
  Chem | Analyze | Elemental Analysis... counts the atoms of every element the molecules hold and
  appends a column per element plus Molecule Charge, leaving the rows as they were. The Radar Viewer
  option adds a viewer. Structural Alerts... appends a boolean column per rule set ticked (PAINS by
  default), named "<rule set> (<molecule column>)"; Pharmacophores... appends the molecules with their
  features highlighted and a boolean column per feature type ticked (Donor, Acceptor, Hydrophobic and
  Aromatic by default).

  V2000 and V3000 molecules are tested in Chem src/tests/screening-tools.ts (elementalAnalysis.molV2000,
  elementalAnalysis.molV3000).

  SMARTS patterns are tested in Chem src/tests/screening-tools.ts 'elementalAnalysis.smarts'. Which molecule
  raises which alert or has which feature, and the rule sets and feature types ticked by hand, are Chem's
  package tests ('calculate: structural alerts', 'pharmacophores').

  Background:
    Given user is logged in
    And the package autostarts have completed

  Scenario: The counts of every element are appended for SMILES molecules
    Given user opens smiles-50 dataset
    When user picks "Chem > Analyze > Elemental Analysis..." from the top menu
    Then "Elemental Analysis" dialog should be visible
    And Molecules input in "Elemental Analysis" dialog should contain text "canonical_smiles"
    And "Radar Viewer" input in "Elemental Analysis" dialog should not be checked
    And "Radar Grid" input in "Elemental Analysis" dialog should not be checked
    When user clicks on OK button in "Elemental Analysis" dialog
    Then the top menu command should have completed
    And a new column "C" should have been added
    And a new column "N" should have been added
    And a new column "O" should have been added
    And a new column "Molecule Charge" should have been added
    And "C" column should have no missing values
    And the table should have 50 rows
    And no errors should have been logged

  Scenario: The Radar Viewer option adds a viewer
    Given user opens smiles-50 dataset
    When user picks "Chem > Analyze > Elemental Analysis..." from the top menu
    And user checks "Radar Viewer" input in "Elemental Analysis" dialog
    And user clicks on OK button in "Elemental Analysis" dialog
    Then the top menu command should have completed
    And a new column "C" should have been added
    And the current view should hold at least 2 viewers
    And no errors should have been logged

  Scenario: Structural Alerts appends a boolean column for the rule set ticked
    Given user opens smiles-50 dataset
    When user picks "Chem > Analyze > Structural Alerts..." from the top menu
    Then "Structural Alerts" dialog should be visible
    And Molecules input in "Structural Alerts" dialog should contain text "canonical_smiles"
    When user clicks on OK button in "Structural Alerts" dialog
    Then a new column "PAINS (canonical_smiles)" should have been added
    And "PAINS (canonical_smiles)" column should have type "bool"

  Scenario: Pharmacophores appends the highlighted molecules and a boolean column per feature type ticked
    Given user opens smiles-50 dataset
    When user picks "Chem > Analyze > Pharmacophores..." from the top menu
    Then "Pharmacophores" dialog should be visible
    And Molecules input in "Pharmacophores" dialog should contain text "canonical_smiles"
    When user clicks on OK button in "Pharmacophores" dialog
    Then 5 new columns should have been added
    And a new column "Pharmacophore (canonical_smiles)" should have been added
    And "Pharmacophore (canonical_smiles)" column should have semantic type "Molecule"
