@journey
Feature: Chem > Analyze > Structural Alerts adds a boolean column per rule set ticked
  Chem > Analyze > Structural Alerts... opens the function's dialog over the current table: the molecule
  column, and a checkbox per rule set (PAINS ticked, BMS, SureChEMBL, MLSMR, Dundee, Inpharmatica, LINT and
  Glaxo not). OK adds one boolean column per ticked rule set, named "<rule set> (<molecule column>)", true in
  the rows whose molecule raises an alert of that set (help `datagrok/solutions/domains/chem/info-panels/structural-alerts.md`
  names the rule sets; the dialog itself has no help page, so the claims follow `structuralAlertsTopMenu` in
  Chem `src/package.ts`). The rows are ethanol (no alert), catechol (PAINS and BMS, not Glaxo) and
  p-benzoquinone (PAINS, BMS and Glaxo), so the two scenarios tell the rule sets apart.

  Left out on purpose:
  - which molecule raises which alert across all 1,251 SMARTS is Chem's `src/tests/screening-tools.ts` and
    `chemical-engine-tests.ts` (both engines); one row per outcome is enough here;
  - the Structural Alerts info pane and its "Calculate for whole dataset" button are in
    `panels/info-panels.feature`.
  Nothing is put on the server.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens a table "alerts" with:
      | name     | smiles           |
      | ethanol  | CCO              |
      | catechol | Oc1ccccc1O       |
      | quinone  | O=C1C=CC(=O)C=C1 |
    And the semantic types of the current table are detected

  Scenario: The dialog opens with PAINS ticked, and OK adds a PAINS column true for the molecules that raise one
    When user picks "Chem > Analyze > Structural Alerts..." from the top menu
    Then "Structural Alerts" dialog should be visible
    And Molecules input in "Structural Alerts" dialog should contain text "smiles"
    And PAINS input in "Structural Alerts" dialog should be checked
    And BMS input in "Structural Alerts" dialog should not be checked
    And Glaxo input in "Structural Alerts" dialog should not be checked
    When user clicks on OK button in "Structural Alerts" dialog
    Then "Structural Alerts" dialog should be hidden
    And a new column "PAINS (smiles)" should have been added
    And 1 new column should have been added
    And "PAINS (smiles)" column should have type "bool"
    And the value of "PAINS (smiles)" column in row 1 should be "false"
    And the value of "PAINS (smiles)" column in row 2 should be "true"
    And the value of "PAINS (smiles)" column in row 3 should be "true"
    And the table should have 3 rows
    And no errors should have been logged
    When user removes "PAINS (smiles)" column
    Then the table should have the columns "name, smiles"

  Scenario: Only the rule sets ticked get a column, each with its own alerts
    When user picks "Chem > Analyze > Structural Alerts..." from the top menu
    And user unchecks PAINS input in "Structural Alerts" dialog
    And user checks BMS input in "Structural Alerts" dialog
    And user checks Glaxo input in "Structural Alerts" dialog
    And user clicks on OK button in "Structural Alerts" dialog
    Then a new column "Glaxo (smiles)" should have been added
    And a new column "BMS (smiles)" should have been added
    And 2 new columns should have been added
    And the value of "BMS (smiles)" column in row 1 should be "false"
    And the value of "BMS (smiles)" column in row 2 should be "true"
    And the value of "BMS (smiles)" column in row 3 should be "true"
    And the value of "Glaxo (smiles)" column in row 1 should be "false"
    And the value of "Glaxo (smiles)" column in row 2 should be "false"
    And the value of "Glaxo (smiles)" column in row 3 should be "true"
    And no errors should have been logged
    When user removes "BMS (smiles)" column
    And user removes "Glaxo (smiles)" column
    Then the table should have the columns "name, smiles"
