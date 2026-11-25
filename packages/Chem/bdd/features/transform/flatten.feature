@journey
Feature: Flatten Molecules keeps the E/Z of double bonds unless asked
  Chem > Transform > Flatten Molecules... removes the stereocenters and, by default, keeps the E/Z of
  double bonds: the flat SMILES column it adds still holds their / and \ bonds. Its "Flatten E/Z" option,
  off by default, flattens the double bonds too.

  The SMILES both engines write, for each option, are Chem's package tests ('chemical engine: column jobs').

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens a table "alkenes" with:
      | smiles             |
      | C/C=C/[C@H](N)C    |
      | F/C=C\F            |
      | C1CC/C=C/CCC1      |
    And the semantic types of the current table are detected

  Scenario: Flattening keeps the E/Z of double bonds by default
    When user picks "Chem > Transform > Flatten Molecules..." from the top menu
    Then "Flatten E/Z" input in "Flatten Molecules" dialog should not be checked
    When user clicks on OK button in "Flatten Molecules" dialog
    Then a new column "smiles_flat" should have been added
    And every value of "smiles_flat" column should match "^[^@]*[/\\][^@]*$"
    And no errors should have been logged

  Scenario: Flatten E/Z flattens the double bonds too
    When user picks "Chem > Transform > Flatten Molecules..." from the top menu
    And user checks "Flatten E/Z" input in "Flatten Molecules" dialog
    And user clicks on OK button in "Flatten Molecules" dialog
    Then a new column "smiles_flat (2)" should have been added
    And every value of "smiles_flat (2)" column should match "^[^@/\\]*$"
    And no errors should have been logged
