@journey @realizes:chem.cp.sketcher-cell-editor
Feature: The sketcher opened from a molecule cell
  A double-click on a canonical_smiles cell of smiles-50 opens the sketcher on that cell's molecule.
  Cyclohexane typed into its molecule field and confirmed with OK replaces the molecule in the cell,
  and reopening the cell brings the sketcher up on the new molecule. What the sketcher holds is read
  the way a user takes it out: Copy as SMILES from its options menu, and the clipboard read by RDKit
  against the cell (the molecule field itself shows only what was typed into it).

  Background:
    Given user is logged in
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And user opens smiles-50 dataset

  Scenario: A double-click opens the sketcher on the cell's molecule
    When user double-clicks on the "cell 1 of canonical_smiles" area of grid
    Then sketcher dialog should be visible
    And molecule input of sketcher dialog should be visible
    When user clicks on "Options" icon in sketcher dialog
    And user picks "Copy as SMILES" from the open menu
    Then the clipboard should hold the molecule of row 1 of "canonical_smiles" column
    And no errors should have been logged

  Scenario: A molecule typed into the sketcher replaces the one in the cell
    When user types "C1CCCCC1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    And the value of "canonical_smiles" column in row 1 should be "C1CCCCC1"
    And the table should have 50 rows
    And no errors should have been logged

  Scenario: Reopening the cell brings the sketcher up on the new molecule
    When user double-clicks on the "cell 1 of canonical_smiles" area of grid
    Then sketcher dialog should be visible
    When user clicks on "Options" icon in sketcher dialog
    And user picks "Copy as SMILES" from the open menu
    Then the clipboard should hold the molecule of row 1 of "canonical_smiles" column
    When user clicks on CANCEL button in sketcher dialog
    Then sketcher dialog should be absent
    And the value of "canonical_smiles" column in row 1 should be "C1CCCCC1"
    And no errors should have been logged
