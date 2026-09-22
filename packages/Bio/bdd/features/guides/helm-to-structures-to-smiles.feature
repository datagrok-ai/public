@guide @help:datagrok/solutions/domains/bio
Feature: HELM sequences to structures, then to SMILES
  A guide: the answer to "how do I convert HELM sequences to molecular structures, and then
  convert the notation to SMILES?". Bio | Transform | To Atomic Level assembles a molecule per
  sequence from the monomer library and adds it as a Molecule column (V3000 molfiles); Chem |
  Transform | Convert Notation then rewrites that column in SMILES as a second new column, with
  SMILES already proposed as the target. Demo: Bio's filter_HELM fixture (four HELM peptides:
  a branched one and three cyclic ones). Film it with BDD_GUIDE_VIEWPORT=1920x1080: at the
  default 1600 px the ribbon's current-cell display folds Bio and Chem under the "more" group.

  Scenario: Build molecules from the HELM sequences, then write them out as SMILES
    Given user is logged in
    And simple mode is off
    And user opens filter_HELM dataset
    And the Bio package is initialized
    When user picks "Bio > Transform > To Atomic Level..." from the top menu
    And user clicks on OK button in "To Atomic Level" dialog
    Then the top menu command should have completed
    And a new column matching "^molfile\(HELM string\)" should have been added
    And every value of "molfile(HELM string)" column should contain "V30 BEGIN CTAB"
    When user picks "Chem > Transform > Convert Notation..." from the top menu
    Then "Target Notation" input in "Convert Notation" dialog should have value "smiles"
    When user clicks on OK button in "Convert Notation" dialog
    Then the top menu command should have completed
    And a new column "molfile(HELM string)_smiles" should have been added
    And every value of "molfile(HELM string)_smiles" column should match "^[A-Za-z0-9@+()\[\]\\/%=#$.:-]+$"
    When user picks "Column Properties..." from the context menu of the "header molfile(HELM string)_smiles" area of grid
    And user enters "canonical_smiles" into "New name" input in "molfile(HELM string)_smiles" dialog
    And user clicks on OK button in "molfile(HELM string)_smiles" dialog
    Then the table should have a column "canonical_smiles"
    And no error or warning balloon should have been shown
