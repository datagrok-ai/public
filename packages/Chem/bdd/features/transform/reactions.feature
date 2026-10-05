@journey @realizes:chem.cp.transform-reactions
Feature: Reactions from the Transform menu
  Remove Water and Salts adds Desalted(canonical_smiles) on smiles: no molecule keeps water or a
  hydrogen halide, none gains a fragment, and some lose one. A halide anion paired with a quaternary
  ammonium cation stays (rows 52, 938, 941, 942), and a second run writes into the same column. Transformation with the
  Ester Hydrolysis card adds a product column in which some molecules, not all, differ from their
  desalted input. Two-Component Reaction with the Amide Coupling card on ten acid/amine
  pairs and two pairs that cannot couple makes an amide for the ten pairs only.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens smiles dataset

  Scenario: Remove Water and Salts strips counter-ions and water
    When user picks "Chem > Transform > Reactions > Remove Water and Salts..." from the top menu
    Then Molecules input in "Remove Water and Salts" dialog should contain text "canonical_smiles"
    When user clicks on OK button in "Remove Water and Salts" dialog
    Then a new column "Desalted(canonical_smiles)" should have been added
    And "Desalted(canonical_smiles)" column should have semantic type "Molecule"
    And "Desalted(canonical_smiles)" column should have no missing values
    And no molecule of "Desalted(canonical_smiles)" column should have more fragments than in "canonical_smiles" column
    And some but not all molecules of "Desalted(canonical_smiles)" column should differ from "canonical_smiles" column
    And no molecule of "Desalted(canonical_smiles)" column should contain "[OH2]"
    And no molecule of "Desalted(canonical_smiles)" column should contain "[ClH,BrH,IH]"
    And no errors should have been logged

  Scenario: Running Remove Water and Salts again writes into the same column
    When user picks "Chem > Transform > Reactions > Remove Water and Salts..." from the top menu
    And user clicks on OK button in "Remove Water and Salts" dialog
    Then the top menu command should have completed
    And no new column should have been added
    And the table should have a column "Desalted(canonical_smiles)"
    And no errors should have been logged

  Scenario: Transformation with ester hydrolysis changes the esters only
    When user picks "Chem > Transform > Reactions > Transformation..." from the top menu
    Then Molecules input in "Run Reaction" dialog should contain text "canonical_smiles"
    And "Remove salts and water" input in "Run Reaction" dialog should be checked
    When user clicks on "Ester Hydrolysis (Saponification)" reaction in "Run Reaction" dialog
    And user clicks on OK button in "Run Reaction" dialog
    Then a new column "Reacted(canonical_smiles)" should have been added
    And "Reacted(canonical_smiles)" column should have no missing values
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
    And the value of "Product(smiles1+smiles2)" column in row 11 should be ""
    And the value of "Product(smiles1+smiles2)" column in row 12 should be ""
    And no molecule of "Product(smiles1+smiles2)" column should be the same as in "smiles1" column
    And no molecule of "Product(smiles1+smiles2)" column should be the same as in "smiles2" column
    And no errors should have been logged
