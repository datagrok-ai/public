@journey @realizes:chem.cp.transform-notation-roundtrip
Feature: Convert Notation and Recalculate Coordinates keep the molecules
  On smiles, Convert Notation to molblock adds a column of V2000 molblocks beside canonical_smiles
  (Overwrite off, Join on), and converting that column back to smiles gives the same molecules row
  for row. Recalculate Coordinates with CoordGen adds a column of molblocks whose coordinates differ
  from the converted ones while the atoms and bonds stay the same.

  The stereochemistry survives in every row but four: converting to molblock inverts the vinyl
  stereocentre of the quinuclidines in rows 478, 487 and 488 (GROK-20956), and CoordGen flips the
  C=N geometry of row 833. The inversion is RDKit's own: its molblock round trip of those SMILES
  flips the centre in RDKit 2024.09 and keeps it in 2026.03, so the tag goes when Chem's
  RDKit_minimal (1.2.23) is upgraded. The main scenarios compare molecules without
  stereochemistry; the last scenario compares them with it and is a known failure.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens smiles dataset

  Scenario: Converting to molblock adds a column of V2000 molblocks
    When user picks "Chem > Transform > Convert Notation..." from the top menu
    Then Molecules input in "Convert Notation" dialog should contain text "canonical_smiles"
    And "Target Notation" input in "Convert Notation" dialog should have value "smiles"
    And Overwrite input in "Convert Notation" dialog should not be checked
    And Join input in "Convert Notation" dialog should be checked
    When user selects "molblock" in "Target Notation" input in "Convert Notation" dialog
    And user clicks on OK button in "Convert Notation" dialog
    Then 1 new column should have been added
    And a new column "canonical_smiles_molblock" should have been added
    And "canonical_smiles_molblock" column should have no missing values
    And every value of "canonical_smiles_molblock" column should contain "V2000"
    And every value of "canonical_smiles_molblock" column should contain "M  END"
    And "canonical_smiles_molblock" column should have semantic type "Molecule"
    And every molecule of "canonical_smiles_molblock" column should be the same as in "canonical_smiles" column ignoring stereochemistry
    And the table should have 1000 rows
    And no errors should have been logged

  Scenario: Converting the molblocks back to smiles gives the same molecules
    When user picks "Chem > Transform > Convert Notation..." from the top menu
    And user selects "canonical_smiles_molblock" in Molecules input in "Convert Notation" dialog
    And user selects "smiles" in "Target Notation" input in "Convert Notation" dialog
    And user clicks on OK button in "Convert Notation" dialog
    Then a new column "canonical_smiles_molblock_smiles" should have been added
    And every molecule of "canonical_smiles_molblock_smiles" column should be the same as in "canonical_smiles" column ignoring stereochemistry
    And no errors should have been logged

  Scenario: Recalculating coordinates with CoordGen moves the atoms and keeps the molecules
    When user picks "Chem > Transform > Recalculate Coordinates..." from the top menu
    Then Molecules input in "Recalculate Coordinates" dialog should contain text "canonical_smiles"
    And Method input in "Recalculate Coordinates" dialog should have value "OCL"
    And Join input in "Recalculate Coordinates" dialog should be checked
    When user selects "CoordGen" in Method input in "Recalculate Coordinates" dialog
    And user clicks on OK button in "Recalculate Coordinates" dialog
    Then a new column "canonical_smiles_recalcCoords" should have been added
    And every value of "canonical_smiles_recalcCoords" column should contain "M  END"
    And some value of "canonical_smiles_recalcCoords" column should differ from "canonical_smiles_molblock" column in the same row
    And every molecule of "canonical_smiles_recalcCoords" column should be the same as in "canonical_smiles" column ignoring stereochemistry
    And the table should have 1000 rows
    And no errors should have been logged

  @known-failure @GROK-20956
  Scenario: Converted molecules keep their stereochemistry
    Then every molecule of "canonical_smiles_molblock" column should be the same as in "canonical_smiles" column
    And every molecule of "canonical_smiles_recalcCoords" column should be the same as in "canonical_smiles" column
