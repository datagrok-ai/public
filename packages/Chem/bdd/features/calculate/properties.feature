@journey @realizes:chem.cp.calculate-properties
Feature: Chemical properties, toxicity risks and InChI from the Calculate menu
  On smiles (1000 molecules in canonical_smiles) the Calculate menu adds columns to the table.
  Chemical Properties opens on the molecule column with MW the only property ticked and no calculator
  ticked, and OK with no calculator asks for one and adds nothing; with the OCL calculator ticked it adds MW;
  run again with every property ticked it adds nine columns, MW under a new name. Toxicity Risks
  opens with Mutagenicity the only risk ticked; each risk column holds Unknown, None, Low or High,
  with at least one Low or High. To InchI adds an inchi column of InChI strings, To InchI Keys an
  inchi_key column of 27-character keys. No command changes the number of rows. A molecule
  OpenChemLib cannot parse keeps its property cells empty.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens smiles dataset

  Scenario: Chemical Properties with no calculator ticked asks for one
    When user picks "Chem > Calculate > Chemical Properties..." from the top menu
    Then "Chemical Properties" dialog should be visible
    And "Chemical Properties (OCL)" calculator in "Chemical Properties" dialog should not be checked
    When user clicks on OK button in "Chemical Properties" dialog
    Then a warning balloon containing "Please select at least one calculation" should have been shown
    And no new column should have been added
    And the table should have 1000 rows

  Scenario: Chemical Properties with MW alone adds an MW column
    When user picks "Chem > Calculate > Chemical Properties..." from the top menu
    Then "Chemical Properties" dialog should be visible
    And Molecules input in "Chemical Properties" dialog should contain text "canonical_smiles"
    And "Chemical Properties (OCL)" calculator in "Chemical Properties" dialog should not be checked
    And MW input in "Chemical Properties" dialog should be checked
    And HBA input in "Chemical Properties" dialog should not be checked
    And "Molecule charge" input in "Chemical Properties" dialog should not be checked
    When user checks "Chemical Properties (OCL)" calculator in "Chemical Properties" dialog
    And user clicks on OK button in "Chemical Properties" dialog
    Then "Chemical Properties" dialog should be hidden
    And 1 new column should have been added
    And a new column "MW" should have been added
    And every value of "MW" column should lie between 1 and 5000
    And the table should have 1000 rows
    And no errors should have been logged

  Scenario: Chemical Properties with every property adds nine columns
    When user picks "Chem > Calculate > Chemical Properties..." from the top menu
    And user checks "Chemical Properties (OCL)" calculator in "Chemical Properties" dialog
    And user checks HBA input in "Chemical Properties" dialog
    And user checks HBD input in "Chemical Properties" dialog
    And user checks "Log P" input in "Chemical Properties" dialog
    And user checks "Log S" input in "Chemical Properties" dialog
    And user checks PSA input in "Chemical Properties" dialog
    And user checks "Rotatable bonds" input in "Chemical Properties" dialog
    And user checks "Stereo centers" input in "Chemical Properties" dialog
    And user checks "Molecule charge" input in "Chemical Properties" dialog
    And user clicks on OK button in "Chemical Properties" dialog
    Then 9 new columns should have been added
    And a new column "MW (2)" should have been added
    And every value of "HBA" column should lie between 0 and 1000
    And every value of "HBD" column should lie between 0 and 1000
    And every value of "LogP" column should lie between -50 and 50
    And every value of "LogS" column should lie between -50 and 50
    And every value of "PSA" column should lie between 0 and 1000
    And every value of "Rotatable bonds" column should lie between 0 and 1000
    And every value of "Stereo centers" column should lie between 0 and 1000
    And every value of "Molecule charge" column should lie between -20 and 20
    And the table should have 1000 rows
    And no errors should have been logged

  Scenario: Toxicity Risks with Mutagenicity alone adds a column of risk levels
    When user picks "Chem > Calculate > Toxicity Risks..." from the top menu
    Then Mutagenicity input in "Toxicity Risks" dialog should be checked
    And Tumorigenicity input in "Toxicity Risks" dialog should not be checked
    And "Irritating effects" input in "Toxicity Risks" dialog should not be checked
    And "Reproductive effects" input in "Toxicity Risks" dialog should not be checked
    When user clicks on OK button in "Toxicity Risks" dialog
    Then 1 new column should have been added
    And a new column "Mutagenicity" should have been added
    And "Mutagenicity" column should have type "string"
    And every value of "Mutagenicity" column should match "^(Unknown|None|Low|High)$"
    And "Mutagenicity" column should have no missing values
    And some value of "Mutagenicity" column should match "^(Low|High)$"
    And the table should have 1000 rows
    And no errors should have been logged

  Scenario: Toxicity Risks with every risk adds four columns of risk levels
    When user picks "Chem > Calculate > Toxicity Risks..." from the top menu
    And user checks Tumorigenicity input in "Toxicity Risks" dialog
    And user checks "Irritating effects" input in "Toxicity Risks" dialog
    And user checks "Reproductive effects" input in "Toxicity Risks" dialog
    And user clicks on OK button in "Toxicity Risks" dialog
    Then 4 new columns should have been added
    And a new column "Mutagenicity (2)" should have been added
    And every value of "Tumorigenicity" column should match "^(Unknown|None|Low|High)$"
    And some value of "Tumorigenicity" column should match "^(Low|High)$"
    And every value of "Irritating effects" column should match "^(Unknown|None|Low|High)$"
    And some value of "Irritating effects" column should match "^(Low|High)$"
    And every value of "Reproductive effects" column should match "^(Unknown|None|Low|High)$"
    And some value of "Reproductive effects" column should match "^(Low|High)$"
    And the table should have 1000 rows
    And no errors should have been logged

  Scenario: To InchI adds a column of InChI strings
    When user picks "Chem > Calculate > To InchI..." from the top menu
    Then Molecules input in "To InchI" dialog should contain text "canonical_smiles"
    When user clicks on OK button in "To InchI" dialog
    Then a new column "inchi" should have been added
    And "inchi" column should have no missing values
    And every value of "inchi" column should match "^InChI=1S/"
    And the table should have 1000 rows
    And no errors should have been logged

  Scenario: To InchI Keys adds a column of InChI keys
    When user picks "Chem > Calculate > To InchI Keys..." from the top menu
    And user clicks on OK button in "To InchI Keys" dialog
    Then a new column "inchi_key" should have been added
    And "inchi_key" column should have no missing values
    And every value of "inchi_key" column should match "^[A-Z]{14}-[A-Z]{10}-[A-Z]$"
    And the table should have 1000 rows
    And no errors should have been logged
