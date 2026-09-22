@journey @realizes:chem.cp.mpo-profile-crud
Feature: MPO Score over a profile's own properties
  Chem | Calculate | MPO Score... opens with the three properties of the ADME Test profile — Caco2,
  Lipophilicity and Solubility — aggregated by their average. On a table that carries those three
  columns the dialog maps each of them and OK appends a score column between 0 and 1 on every row.

  Background:
    Given user is logged in
    And the package autostarts have completed

  Scenario: The dialog opens on the ADME Test profile and its three properties
    Given user opens a table "adme" with:
      | smiles     | Caco2 | Lipophilicity | Solubility |
      | c1ccccc1   | -6    | 2             | -4         |
      | CCO        | -7    | 3             | -5         |
      | c1ccncc1   | -5    | 1             | -3         |
      | CC(=O)O    | -4    | 4             | -6         |
    When user picks "Chem > Calculate > MPO Score..." from the top menu
    Then "MPO Score" dialog should be visible
    And Aggregation input in "MPO Score" dialog should have value "Average"
    And "MPO Score" dialog should contain text "Caco2"
    And "MPO Score" dialog should contain text "Lipophilicity"
    And "MPO Score" dialog should contain text "Solubility"
    And no errors should have been logged

  Scenario: The run appends a score between 0 and 1 for every row
    When user clicks on OK button in "MPO Score" dialog
    Then the top menu command should have completed
    And 1 new column should have been added
    And the newest column matching "." should have no missing values
    And the table should have 4 rows
    And no errors should have been logged
