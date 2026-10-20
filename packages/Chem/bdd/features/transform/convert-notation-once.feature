@journey @realizes:GROK-17964
Feature: Convert Notation is offered once in a molecule column's Actions
  The context panel of canonical_smiles in smiles-50 lists the Convert Notation... action once, and
  still once after the action's dialog is cancelled, after a conversion to molblock has run from it,
  and after the dialog has been opened and cancelled twice more (GROK-17964).

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens smiles-50 dataset
    And the context panel is open

  Scenario: The action is listed once for the molecule column
    Given the "canonical_smiles" column is the current object
    Then the context panel should show "canonical_smiles"
    When user expands Actions accordion header in context panel
    Then there should be 1 visible "Convert Notation..." action in context panel
    And no errors should have been logged

  Scenario: Cancelling the action's dialog leaves one action
    When user clicks on "Convert Notation..." action in context panel
    Then "Convert Notation" dialog should be visible
    And Overwrite input in "Convert Notation" dialog should not be checked
    And Join input in "Convert Notation" dialog should be checked
    And Kekulize input in "Convert Notation" dialog should not be checked
    When user clicks on CANCEL button in "Convert Notation" dialog
    And user clicks on the "header canonical_smiles" area of grid
    Then the context panel should show "canonical_smiles"
    When user expands Actions accordion header in context panel
    Then there should be 1 visible "Convert Notation..." action in context panel
    And no errors should have been logged

  Scenario: Running a conversion from the action leaves one action
    When user clicks on "Convert Notation..." action in context panel
    And user selects "molblock" in "Target Notation" input in "Convert Notation" dialog
    And user clicks on OK button in "Convert Notation" dialog
    Then "Convert Notation" dialog should be hidden
    And the table should have a column "canonical_smiles_molblock"
    Given the "canonical_smiles" column is the current object
    Then the context panel should show "canonical_smiles"
    When user expands Actions accordion header in context panel
    Then there should be 1 visible "Convert Notation..." action in context panel
    And no errors should have been logged

  Scenario: Opening and cancelling the dialog twice more leaves one action
    When user clicks on "Convert Notation..." action in context panel
    And user clicks on CANCEL button in "Convert Notation" dialog
    And user clicks on "Convert Notation..." action in context panel
    And user clicks on CANCEL button in "Convert Notation" dialog
    # the conversion scrolled the grid to its new column: the header of canonical_smiles is out of view
    Given the "canonical_smiles" column is the current object
    Then the context panel should show "canonical_smiles"
    When user expands Actions accordion header in context panel
    Then there should be 1 visible "Convert Notation..." action in context panel
    And no errors should have been logged
