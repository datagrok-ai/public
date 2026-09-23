@journey @realizes:chem.cp.transform-curate-mutate
Feature: Curate and Mutate from the Transform menu
  On chem_standards, 14 molecules written both as salts and as their standardised parents, Chem |
  Transform | Curate... opens with normalization, reionization, neutralization and main fragment on
  and kekulization and tautomerization off. OK joins a curated_molecule column: no row gains a
  fragment, no counter-ion of succinic acid is left, and some but not all rows hold another molecule
  than they started with — an output that repeats every input means nothing was standardised. A
  second run with kekulization on joins its own column of the same molecules.

  Mutate takes a single molecule — the dialog's own default — with steps 1, randomize on and max
  random results 100, and opens a mutations table of 100 molecules, most of them other than the
  parent.

  Both are server-side Python scripts of the Chem compute service.

  Background:
    Given user is logged in
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And user opens chem_standards dataset

  Scenario: Curate with its default options standardises the salts and leaves the parents alone
    When user picks "Chem > Transform > Curate..." from the top menu
    Then "Curate" dialog should be visible
    And Molecules input in "Curate" dialog should contain text "smiles"
    And Normalization input in "Curate" dialog should be checked
    And Reionization input in "Curate" dialog should be checked
    And Neutralization input in "Curate" dialog should be checked
    And "Main fragment" input in "Curate" dialog should be checked
    And Kekulization input in "Curate" dialog should not be checked
    And Tautomerization input in "Curate" dialog should not be checked
    When user clicks on OK button in "Curate" dialog
    Then "Curate" dialog should be hidden
    And the top menu command should have completed
    And 1 new column should have been added
    And a new column "curated_molecule" should have been added
    And "curated_molecule" column should have no missing values
    And "curated_molecule" column should have semantic type "Molecule"
    And no molecule of "curated_molecule" column should have more fragments than in "smiles" column
    And some but not all molecules of "curated_molecule" column should differ from "smiles" column
    And 9 molecules of "curated_molecule" column should differ from "smiles" column
    And no molecule of "curated_molecule" column should contain "OC(=O)CCC(=O)O"
    And the table should have 14 rows
    And no errors should have been logged

  Scenario: A second run with kekulization joins its own column
    When user picks "Chem > Transform > Curate..." from the top menu
    And user checks Kekulization input in "Curate" dialog
    And user clicks on OK button in "Curate" dialog
    Then the top menu command should have completed
    And a new column "curated_molecule (2)" should have been added
    And "curated_molecule (2)" column should have no missing values
    And every molecule of "curated_molecule (2)" column should be the same as in "curated_molecule" column
    And the table should have 14 rows
    And no errors should have been logged

  Scenario: Mutate returns a hundred molecules
    When user picks "Chem > Transform > Mutate..." from the top menu
    Then "Mutate" dialog should be visible
    And Steps input in "Mutate" dialog should have value "1"
    And Randomize input in "Mutate" dialog should be checked
    And "Max random results" input in "Mutate" dialog should have value "100"
    When user clicks on OK button in "Mutate" dialog
    Then "Mutate" dialog should be hidden
    And the top menu command should have completed
    And table "mutations" should be open
    When user switches to the "mutations" table view
    Then the table should have 100 rows
    And "mutations" column should have no missing values
    And "mutations" column should have semantic type "Molecule"
    And "mutations" column should have at least 20 distinct values
    And no errors should have been logged
    # the second run opens a table of the same name: this one must not answer for it
    When user closes the current view
    And user switches to the "chem_standards" table view

  Scenario: Mutate with two steps still returns a hundred molecules
    When user picks "Chem > Transform > Mutate..." from the top menu
    And user enters "2" into Steps input in "Mutate" dialog
    And user clicks on OK button in "Mutate" dialog
    Then the top menu command should have completed
    And table "mutations" should be open
    When user switches to the "mutations" table view
    Then the table should have 100 rows
    And "mutations" column should have no missing values
    And "mutations" column should have at least 20 distinct values
    And no errors should have been logged
