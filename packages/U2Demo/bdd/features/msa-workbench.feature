@demo @realizes:u2.dialog @realizes:u2.tabs
Feature: MSA workbench
  The sample app: a toolbar, an alignment form, an MSA dialog and a results panel. Every element
  below resolves through the generic kinds — nothing in this feature is registered by name except
  the two panels and the parts of the results panel.

  Background:
    Given user opens the MSA workbench

  Scenario: Running an alignment from the dialog
    When user types "Kinase panel" into name input in alignment panel
    And user selects "peptide" in sequence column input in alignment panel
    And user clicks on run msa button in alignment panel
    Then MSA dialog should be visible
    And sequence column input in MSA dialog should have value "peptide"
    When user selects "muscle" in method input in MSA dialog
    And user clicks on OK button in MSA dialog
    Then MSA dialog should be hidden
    And results should be visible
    And status line should contain text "Aligned 5 sequences with muscle"
    And aligned sequences list should have 5 items
    And notification should contain text "Alignment finished"

  Scenario: Cancelling the dialog keeps the results untouched
    When user clicks on run msa button in alignment panel
    And user clicks on cancel button in MSA dialog
    Then MSA dialog should be hidden
    And status line should have text "No alignment yet"
    And aligned sequences list should have 0 items

  Scenario: Switching result tabs
    When user selects "muscle" in method input in alignment panel
    And user clicks on run msa button in alignment panel
    And user clicks on OK button in MSA dialog
    And user clicks on Log tab in results
    Then Log tab in results should be selected
    And log should contain text "muscle"
    And Alignment tab in results should not be selected
