@demo @realizes:u2.inputs
Feature: Chained choices
  Five choices, each loading its options from every pick above it, from a mock server that answers in
  500 ms (pages/chained-choices.ts). A choice is off, saying why, while a pick above it is missing; a
  change reloads everything below it, and a pick the new list lacks is cleared.

  Background:
    Given user opens the "Chained choices" demo page

  Scenario: Each pick loads the next choice's options
    Then Letter input should be disabled because "Pick Name first"
    And Shape input should be disabled because "Pick Name first"
    When user selects "Ada" in Name input
    And user selects "Ada a" in Letter input
    And user selects "Ada a 1" in Digit input
    Then value of picked readout should have text "Ada → Ada a → Ada a 1 → – → –"
    And Sign input should list the items "Ada a 1 +, Ada a 1 =" in this order
    And Shape input should be disabled because "Pick Sign first"

  Scenario: Changing a pick clears what was built on it and reloads the lists below
    When user selects "Ada" in Name input
    And user selects "Ada b" in Letter input
    And user selects "Ada b 3" in Digit input
    And user selects "Cleo" in Name input
    Then value of picked readout should have text "Cleo → – → – → – → –"
    And Letter input should list the items "Cleo a, Cleo e" in this order
    And Digit input should be disabled because "Pick Letter first"
