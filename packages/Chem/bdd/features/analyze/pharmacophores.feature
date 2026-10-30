@journey
Feature: Chem > Analyze > Pharmacophores adds a highlighted molecule column and a boolean column per feature type
  Chem > Analyze > Pharmacophores... opens the function's dialog over the current table: the molecule column
  and a checkbox per feature type (Donor, Acceptor, Hydrophobic and Aromatic ticked; Positive, Negative and
  Halogen Bond not). OK adds a molecule column "Pharmacophore (<molecule column>)" with the features
  highlighted, and a boolean column per ticked type, "<type> (<molecule column>)", true in the rows whose
  molecule has that feature (help `datagrok/solutions/domains/chem/chem.md`, "Pharmacophores" > How to use).
  The rows are ethanol (donor, acceptor, not aromatic), benzene (aromatic, no donor), acetic acid (negative
  ionizable) and ethylamine (positive ionizable), so each claimed column has a row of each value.

  Left out on purpose:
  - the colours the highlight column draws are its cell renderer's (Chem `src/tests/pharmacophore-features-tests.ts`);
    the feature claims the column, its type and its molecules;
  - which SMARTS of which family matches which molecule is the same test file's matrix, run in both engines;
  - the Pharmacophores info pane is in `panels/info-panels.feature`.
  Nothing is put on the server.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens a table "pharm" with:
      | name       | smiles   |
      | ethanol    | CCO      |
      | benzene    | c1ccccc1 |
      | acid       | CC(=O)O  |
      | ethylamine | CCN      |
    And the semantic types of the current table are detected

  Scenario: The dialog opens with four feature types ticked, and OK adds the highlight column and a column per type
    When user picks "Chem > Analyze > Pharmacophores..." from the top menu
    Then "Pharmacophores" dialog should be visible
    And Molecules input in "Pharmacophores" dialog should contain text "smiles"
    And Donor input in "Pharmacophores" dialog should be checked
    And Acceptor input in "Pharmacophores" dialog should be checked
    And Hydrophobic input in "Pharmacophores" dialog should be checked
    And Aromatic input in "Pharmacophores" dialog should be checked
    And Positive input in "Pharmacophores" dialog should not be checked
    And Negative input in "Pharmacophores" dialog should not be checked
    And "Halogen Bond" input in "Pharmacophores" dialog should not be checked
    When user clicks on OK button in "Pharmacophores" dialog
    Then "Pharmacophores" dialog should be hidden
    And a new column "Pharmacophore (smiles)" should have been added
    And 5 new columns should have been added
    And "Pharmacophore (smiles)" column should have semantic type "Molecule"
    And "Pharmacophore (smiles)" column should hold the same values as "smiles" column
    And "Donor (smiles)" column should have type "bool"
    And the value of "Donor (smiles)" column in row 1 should be "true"
    And the value of "Donor (smiles)" column in row 2 should be "false"
    And the value of "Acceptor (smiles)" column in row 1 should be "true"
    And the value of "Acceptor (smiles)" column in row 2 should be "false"
    And the value of "Aromatic (smiles)" column in row 1 should be "false"
    And the value of "Aromatic (smiles)" column in row 2 should be "true"
    And the value of "Hydrophobic (smiles)" column in row 2 should be "true"
    And the table should have 4 rows
    And no errors should have been logged
    When user removes "Pharmacophore (smiles)" column
    And user removes "Donor (smiles)" column
    And user removes "Acceptor (smiles)" column
    And user removes "Aromatic (smiles)" column
    And user removes "Hydrophobic (smiles)" column
    Then the table should have the columns "name, smiles"

  Scenario: Only the feature types ticked get a column
    When user picks "Chem > Analyze > Pharmacophores..." from the top menu
    And user unchecks Donor input in "Pharmacophores" dialog
    And user unchecks Acceptor input in "Pharmacophores" dialog
    And user unchecks Hydrophobic input in "Pharmacophores" dialog
    And user unchecks Aromatic input in "Pharmacophores" dialog
    And user checks Positive input in "Pharmacophores" dialog
    And user checks Negative input in "Pharmacophores" dialog
    And user clicks on OK button in "Pharmacophores" dialog
    Then a new column "Pharmacophore (smiles)" should have been added
    And a new column "Positive (smiles)" should have been added
    And a new column "Negative (smiles)" should have been added
    And 3 new columns should have been added
    And the value of "Positive (smiles)" column in row 3 should be "false"
    And the value of "Positive (smiles)" column in row 4 should be "true"
    And the value of "Negative (smiles)" column in row 3 should be "true"
    And the value of "Negative (smiles)" column in row 4 should be "false"
    And no errors should have been logged
    When user removes "Pharmacophore (smiles)" column
    And user removes "Positive (smiles)" column
    And user removes "Negative (smiles)" column
    Then the table should have the columns "name, smiles"
