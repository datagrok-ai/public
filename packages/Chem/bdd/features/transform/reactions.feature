@journey @realizes:chem.cp.transform-reactions
Feature: Reactions from the Transform menu
  Remove Water and Salts adds Desalted(canonical_smiles) on smiles. Transformation with the Ester
  Hydrolysis card adds a product column in which some molecules, not all, differ from their desalted
  input. Two-Component Reaction with the Amide Coupling card on ten acid/amine pairs and two pairs
  that cannot couple makes an amide for the ten pairs.

  The molecules each reaction makes (the water and hydrogen halides gone, the halide anions of quaternary ammonium cations kept, a second desalting into the same column, the products of every row) are tested in Chem src/tests/reactions-tests.ts ('Reactions').

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens smiles dataset

  Scenario: Remove Water and Salts adds the desalted molecules
    When user picks "Chem > Transform > Reactions > Remove Water and Salts..." from the top menu
    Then Molecules input in "Remove Water and Salts" dialog should contain text "canonical_smiles"
    When user clicks on OK button in "Remove Water and Salts" dialog
    Then a new column "Desalted(canonical_smiles)" should have been added
    And "Desalted(canonical_smiles)" column should have semantic type "Molecule"
    And no errors should have been logged

  Scenario: Transformation with ester hydrolysis changes some molecules
    When user picks "Chem > Transform > Reactions > Transformation..." from the top menu
    Then Molecules input in "Run Reaction" dialog should contain text "canonical_smiles"
    And "Remove salts and water" input in "Run Reaction" dialog should be checked
    When user clicks on "Ester Hydrolysis (Saponification)" reaction in "Run Reaction" dialog
    And user clicks on OK button in "Run Reaction" dialog
    Then a new column "Reacted(canonical_smiles)" should have been added
    And some but not all molecules of "Reacted(canonical_smiles)" column should differ from "Desalted(canonical_smiles)" column
    And no errors should have been logged

  Scenario: Two-Component Reaction couples the acids with the amines
    Given user opens a table "amide_coupling" with:
      | smiles1               | smiles2             |
      | OC(=O)c1ccccc1        | Nc1ccccc1           |
      | CC(=O)O               | NCc1ccccc1          |
      | CCC(=O)O              | NC1CCCCC1           |
      | CCCC(=O)O             | CCCCN               |
      | OC(=O)C1CCCCC1        | Cc1ccc(N)cc1        |
      | Cc1ccc(cc1)C(=O)O     | Nc1ccc(Cl)cc1       |
      | OC(=O)c1ccc(Cl)cc1    | NCCc1ccccc1         |
      | OC(=O)Cc1ccccc1       | Nc1cccnc1           |
      | OC(=O)c1cccnc1        | CCN                 |
      | OC(=O)c1ccco1         | COc1ccc(N)cc1       |
      | Cc1ccccc1             | Nc1ccccc1           |
      | OC(=O)c1ccccc1        | COc1ccccc1          |
    When user picks "Chem > Transform > Reactions > Two-Component Reaction..." from the top menu
    Then "Reactant 1" input in "Two-Component Reaction" dialog should contain text "smiles1"
    And "Combination Mode" input in "Two-Component Reaction" dialog should have value "pairwise"
    When user selects "smiles2" in "Reactant 2" input in "Two-Component Reaction" dialog
    And user clicks on "Amide Coupling" reaction in "Two-Component Reaction" dialog
    And user clicks on OK button in "Two-Component Reaction" dialog
    Then 1 new column should have been added
    And a new column "Product(smiles1+smiles2)" should have been added
    And 10 molecules of "Product(smiles1+smiles2)" column should contain "[CX3](=O)[NX3;H1]"
    And no errors should have been logged
