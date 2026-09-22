@journey @realizes:chem.cp.calculate-identifiers
Feature: Map Identifiers and Generate Conformers from the Calculate menu
  Map Identifiers opens on the open table with its molecule column in Ids and "smiles" in From
  Source. With the chem-chem service running, To Source "inchi_key" appends an inchi_key column of
  27-character keys and To Source "smiles" appends a canonical smiles column holding, row by row,
  the same molecule as the column it was made from.

  Generate Conformers runs the RDKit ETKDGv3 script on a single molecule — the dialog's own default,
  butane — with Num conformers 50, Optimize on, RMS threshold 0.1, Max attempts 5000 and Random seed
  42. It opens a conformers table whose smiles column repeats the input molecule, whose conformer
  numbers start at 1, whose molblocks carry coordinates and end in "M  END", whose MMFF94 energies
  are all present, and whose rmsd reads 0 on the reference row and something else further down.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens smiles-50 dataset

  Scenario: Map Identifiers appends InChI keys for the chosen source
    Given the "chem-chem" container is running
    When user picks "Chem > Calculate > Map Identifiers..." from the top menu
    Then "Map Identifiers" dialog should be visible
    And Ids input in "Map Identifiers" dialog should contain text "canonical_smiles"
    And "From Source" input in "Map Identifiers" dialog should have value "smiles"
    When user selects "inchi_key" in "To Source" input in "Map Identifiers" dialog
    And user clicks on OK button in "Map Identifiers" dialog
    Then "Map Identifiers" dialog should be hidden
    And the top menu command should have completed
    And 1 new column should have been added
    And a new column "inchi_key" should have been added
    And "inchi_key" column should have no missing values
    And every value of "inchi_key" column should match "^[A-Z]{14}-[A-Z]{10}-[A-Z]$"
    And every value of "inchi_key" column should have the same length
    And the table should have 50 rows
    And no errors should have been logged

  Scenario: Map Identifiers to smiles returns the same molecules
    Given the "chem-chem" container is running
    When user picks "Chem > Calculate > Map Identifiers..." from the top menu
    And user selects "smiles" in "To Source" input in "Map Identifiers" dialog
    And user clicks on OK button in "Map Identifiers" dialog
    Then the top menu command should have completed
    And a new column "smiles" should have been added
    And "smiles" column should have no missing values
    And every molecule of "smiles" column should be the same as in "canonical_smiles" column
    And the table should have 50 rows
    And no errors should have been logged

  Scenario: Generate Conformers builds a conformer table for the dialog's own molecule
    When user picks "Chem > Calculate > Generate Conformers..." from the top menu
    Then "Generate Conformers" dialog should be visible
    And "Num conformers" input in "Generate Conformers" dialog should have value "50"
    And Optimize input in "Generate Conformers" dialog should be checked
    And "RMS threshold" input in "Generate Conformers" dialog should have value "0.1"
    And "Max attempts" input in "Generate Conformers" dialog should have value "5000"
    And "Random seed" input in "Generate Conformers" dialog should have value "42"
    When user clicks on OK button in "Generate Conformers" dialog
    Then "Generate Conformers" dialog should be hidden
    And the top menu command should have completed
    And table "conformers" should be open
    And table "conformers" should have columns "smiles, molblock, conformer, energy, rmsd"
    When user switches to the "conformers" table view
    Then every value of "smiles" column should match "^CCCC$"
    And the value of "conformer" column in row 1 should be "1"
    And "conformer" column should have no missing values
    And "conformer" column should have at least 2 distinct values
    And every value of "molblock" column should contain "M  END"
    And "molblock" column should have no missing values
    And "energy" column should have no missing values
    And "rmsd" column should have no missing values
    And "rmsd" column should have at least 2 distinct values
    And every value of "rmsd" column should lie between 0 and 100
    And no errors should have been logged
