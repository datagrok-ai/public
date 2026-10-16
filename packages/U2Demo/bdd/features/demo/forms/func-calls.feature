@demo @realizes:u2.func-form
Feature: FuncCalls
  One FuncCall edited by the platform's DG.InputForm.forFuncCall in the Platform section and by u2's
  funcForm in the u2 section. An edit on either side writes the call, the other side follows and the
  sync line counts it once; the required parameter of the second function gates isValid on both
  forms alike.

  Background:
    Given user opens the "FuncCalls" demo page

  Scenario: The two sides never lie on one another, side by side or stacked
    Then Stage input in "u2" section should not overlap Stage input in "Platform" section
    And editor of Stage input in "u2" section should be as wide as editor of Stage input in "Platform" section
    And demo page should have no layout problems

  Scenario: A stage picked on either side writes the call and shows on the other
    When user selects "Phase II" in Stage input in "u2" section
    Then Stage input in "Platform" section should have value "Phase II"
    And value of first inputs readout should contain text "stage: 'Phase II'"
    And value of first sync readout should have text "u2 → Dart 1 · Dart → u2 0"
    When user selects "Discovery" in Stage input in "Platform" section
    Then Stage input in "u2" section should have value "Discovery"
    And value of first inputs readout should contain text "stage: 'Discovery'"
    And value of first sync readout should have text "u2 → Dart 1 · Dart → u2 1"

  Scenario: Emptying the required parameter makes both forms invalid
    Then value of isValid readout should have text "Dart true · u2 false"
    When user enters "B-17" into "Batch Id" input in "u2" section
    Then "Batch Id" input in "Platform" section should have value "B-17"
    And "Batch Id" input in "Platform" section should be valid
    And value of isValid readout should have text "Dart true · u2 true"
    When user clears "Batch Id" input in "u2" section
    Then "Batch Id" input in "Platform" section should be invalid
    And "Batch Id" input in "u2" section should have the error "Value can't be empty"
    And value of isValid readout should have text "Dart false · u2 false"
    When user enters "B-18" into "Batch Id" input in "Platform" section
    Then "Batch Id" input in "u2" section should have value "B-18"
    And "Batch Id" input in "u2" section should be valid
    And value of isValid readout should have text "Dart true · u2 true"
