@demo @realizes:u2.inputs
Feature: Basic inputs
  Every u2 input resolves by its label through the generic kinds — nothing on this page carries
  a name. The "name = value" lines are readouts, a kind this package registers for itself in
  bindings/demo.ts, so "value of search readout" is the page's own signal state.

  Background:
    Given user opens the "Basic inputs" demo page

  Scenario: Text inputs and a bound search
    Then Name input should have value "Aspirin"
    And label of Name input should have text "Name"
    And value of search readout should have text "(empty)"
    When user types "Kinase" into Search input
    Then value of search readout should have text "Kinase"
    And Preview input should contain text "Kinase"
    When user clears Search input
    Then value of search readout should have text "(empty)"

  Scenario: A text area
    When user types "first line" into Notes text area
    Then Notes text area should have value "first line"

  Scenario: Checkboxes and switches
    Then Enabled checkbox should be checked
    And Notifications checkbox should be checked
    When user unchecks Notifications checkbox
    Then Notifications checkbox should be unchecked
    And value of notifications readout should have text "false"
    When user toggles Enabled checkbox
    Then Enabled checkbox should be unchecked

  Scenario: Numbers feed a computed readout
    Then value of "dose * replicates" readout should have text "750"
    When user enters "4" into Replicates input
    Then value of "dose * replicates" readout should have text "1000"
    When user enters "100" into "Dose, mg" input
    Then value of "dose * replicates" readout should have text "400"

  Scenario: Single and multiple choices
    Then Stage input should have value "Discovery"
    When user selects "Phase II" in Stage input
    Then Stage input should have value "Phase II"
    And Kinase checkbox in Targets input should be checked
    When user checks Protease checkbox in Targets input
    And user unchecks Kinase checkbox in Targets input
    Then Protease checkbox in Targets input should be checked
    And Kinase checkbox in Targets input should be unchecked

  Scenario: Validation shows on the input and aggregates into a readout
    Then error of Code input should have text "Value is required"
    And value of "code validity" readout should have text "Value is required"
    When user types "ABC-123" into Code input
    Then error of Code input should be hidden
    And value of "code validity" readout should have text "valid"
    When user types "ABCDEFGHIJKLMNOP" into Code input
    Then error of Code input should have text "At most 10 characters"

  Scenario: Pickers that open a popup
    When user selects "acorn" in Icon input
    Then value of icon readout should have text "acorn"
    When user types "Ky" into City input
    And user selects "Kyiv" in City input
    Then City input should have value "Kyiv"
    And value of city readout should have text "Kyiv"
    When user selects "Abs" in Scorer input
    Then value of scorer readout should contain text "Abs"
