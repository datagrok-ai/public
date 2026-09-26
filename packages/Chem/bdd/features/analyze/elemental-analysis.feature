@journey @realizes:chem.cp.elemental-analysis
Feature: Elemental Analysis over SMILES, V2000 and V3000 molecules
  Chem | Analyze | Elemental Analysis... counts the atoms of every element the molecules hold and
  appends a column per element plus Molecule Charge, leaving the rows as they were. It does the same
  for molecules read from a V2000 SDF and from a V3000 SDF. The Radar Viewer option adds a viewer.

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
    And every value of "C" column should lie between 1 and 200
    And "C" column should have type "int"
    And the table should have 50 rows
    And no errors should have been logged

  Scenario: The same counts come out of V2000 molecules
    Given user opens mol1K.sdf dataset
    When user picks "Chem > Analyze > Elemental Analysis..." from the top menu
    And user clicks on OK button in "Elemental Analysis" dialog
    Then the top menu command should have completed
    And a new column "C" should have been added
    And a new column "Molecule Charge" should have been added
    And "C" column should have no missing values
    And every value of "C" column should lie between 1 and 200
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The same counts come out of V3000 molecules
    Given user opens ApprovedDrugs2015 dataset
    When user picks "Chem > Analyze > Elemental Analysis..." from the top menu
    And user clicks on OK button in "Elemental Analysis" dialog
    Then the top menu command should have completed
    And a new column "C" should have been added
    And a new column "Molecule Charge" should have been added
    And "C" column should have no missing values
    And no error or warning balloon should have been shown
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

  Scenario: SMARTS patterns are counted the same way
    Given user opens ex-smarts dataset
    When user picks "Chem > Analyze > Elemental Analysis..." from the top menu
    Then "Elemental Analysis" dialog should be visible
    When user clicks on OK button in "Elemental Analysis" dialog
    Then the top menu command should have completed
    And a new column "Molecule Charge" should have been added
    And no error or warning balloon should have been shown
    And no errors should have been logged
