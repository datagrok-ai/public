@journey @realizes:chem.cp.mpo-profile-crud
Feature: MPO Score over a profile's own properties
  Chem | Calculate | MPO Score... on a table that carries Caco2, Lipophilicity and Solubility columns
  offers the package's ADME Test profile, whose properties are exactly those three, marked as fitting
  the table. Picked, it maps the three columns (OK becomes enabled) and aggregates by the average; OK
  appends the "MPO ADME Test" score column, between 0 and 1 on every row.

  The profile is picked rather than claimed as the default: the dialog opens on the first fitting
  profile by name, which depends on the profiles of the stand (with the shipped ones seeded it is
  "ADME Profile with Categorical Property", whose Type property this table has no column for). The
  dialog's text is no evidence of a profile either: every property row lists the table's numeric
  columns in its weight choice. The menu command returns once the dialog is shown, so the run is
  over when its column is there. The profiles live in the mpo domain table, which a stand gets from
  the package's mpo folder through Chem:seedMpoProfiles; a stand that was never seeded is seeded by
  the Background, and keeps the shipped profiles.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the shipped MPO profile "ADME Test" is on the server

  Scenario: The dialog offers the ADME Test profile, which maps the table's three columns
    Given user opens a table "adme" with:
      | smiles     | Caco2 | Lipophilicity | Solubility |
      | c1ccccc1   | -6    | 2             | -4         |
      | CCO        | -7    | 3             | -5         |
      | c1ccncc1   | -5    | 1             | -3         |
      | CC(=O)O    | -4    | 4             | -6         |
    When user picks "Chem > Calculate > MPO Score..." from the top menu
    Then "MPO Score" dialog should be visible
    When user selects "✓ ADME Test" in Profile input in "MPO Score" dialog
    Then Aggregation input in "MPO Score" dialog should have value "Average"
    And OK button in "MPO Score" dialog should be enabled
    And no errors should have been logged

  Scenario: The run appends the profile's score between 0 and 1 for every row
    When user clicks on OK button in "MPO Score" dialog
    Then a new column "MPO ADME Test" should have been added
    And 1 new column should have been added
    And "MPO ADME Test" column should have no missing values
    And every value of "MPO ADME Test" column should lie between 0 and 1
    And the table should have 4 rows
    And no errors should have been logged
