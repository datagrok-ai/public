@journey @realizes:helm.panel.properties
Feature: The Properties pane of a HELM cell
  The context panel of the current HELM cell has a Properties pane with the sequence's molecular
  formula, molecular weight and extinction coefficient; it follows the current cell, and a
  sequence over 1000 characters gets a warning instead of a calculation.

  Not translated, and why: nothing of the manual cases is left out; the old spec's direct
  Helm:propertiesWidget call is replaced by the pane itself.

  Background:
    Given user is logged in
    And the Helm package is initialized
    And user opens helm-showcase dataset
    And the context panel is open
    Then "HELM" column should have units "helm"

  Scenario: The current cell's formula, weight and extinction coefficient
    When user clicks on the "cell 1 of HELM" area of grid
    Then row 1 should be current
    And "Properties" section in context panel should be visible
    Given "Properties" section in context panel is expanded
    Then "Properties" section in context panel should contain text "C6H12N2O3S"
    And "Properties" section in context panel should contain text "192.23"
    And "Properties" section in context panel should contain text "0.06"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Another cell shows its own values
    When user clicks on the "cell 2 of HELM" area of grid
    Then row 2 should be current
    And "Properties" section in context panel should contain text "C50H77N13O15S"
    And "Properties" section in context panel should contain text "1132.30"
    And "Properties" section in context panel should not contain text "C6H12N2O3S"
    And no errors should have been logged

  Scenario: A sequence over 1000 characters gets a warning, not a calculation
    When user sets "HELM" column in row 5 to a peptide of 600 alanines
    And user clicks on the "cell 5 of HELM" area of grid
    Then row 5 should be current
    And "Properties" section in context panel should contain text "Too long sequence"
    And "Properties" section in context panel should not contain text "formula"
    And no error or warning balloon should have been shown
    And no errors should have been logged
