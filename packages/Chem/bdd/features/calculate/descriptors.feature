@journey @realizes:chem.cp.calculate-descriptors-docker
Feature: Descriptors from the chem-chem service
  With the chem-chem container running, Chem | Calculate | Descriptors (RDKit)... offers the
  descriptor tree of the service. None clears the selection, MolWt and MolLogP tick two descriptors,
  and OK appends those two columns to smiles-50, filled on every row, with the rows left alone.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens smiles-50 dataset

  Scenario: The dialog offers the descriptor tree and appends the two chosen descriptors
    Given the "chem-chem" container is running
    When user picks "Chem > Calculate > Descriptors (RDKit)..." from the top menu
    Then "Chemical Descriptors" dialog should be visible
    And Molecules input in "Chemical Descriptors" dialog should contain text "canonical_smiles"
    When user clicks on None label in "Chemical Descriptors" dialog
    Then "Chemical Descriptors" dialog should contain text "0 checked"
    When user expands Descriptors tree node in "Chemical Descriptors" dialog
    And user checks MolWt tree node in "Chemical Descriptors" dialog
    And user expands Crippen tree node in "Chemical Descriptors" dialog
    And user checks MolLogP tree node in "Chemical Descriptors" dialog
    Then "Chemical Descriptors" dialog should contain text "2 checked"
    When user clicks on OK button in "Chemical Descriptors" dialog
    Then the top menu command should have completed
    And 2 new columns should have been added
    And a new column "MolWt" should have been added
    And a new column "MolLogP" should have been added
    And "MolWt" column should have no missing values
    And every value of "MolWt" column should lie between 1 and 5000
    And "MolLogP" column should have no missing values
    And the table should have 50 rows
    And no errors should have been logged
