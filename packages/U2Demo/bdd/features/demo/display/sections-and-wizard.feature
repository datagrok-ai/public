@demo @realizes:u2.section @realizes:u2.wizard
Feature: Sections and the wizard

  Background:
    Given user opens the "Sections & wizard" demo page

  Scenario: A collapsible section keeps its content
    Then Advanced section should be collapsed
    And Threshold input should be hidden
    When user expands Advanced section
    Then Threshold input should have value "0.85"
    And value of expanded readout should have text "true"
    When user collapses Advanced section
    Then Threshold input should be hidden

  Scenario: The wizard gates its steps and keeps their state
    Then Name wizard step should be selected
    And NEXT button in wizard should be disabled
    And wizard should contain text "Enter a project name"
    When user types "Atlas" into Project input in wizard
    Then NEXT button in wizard should be enabled
    When user clicks on NEXT button in wizard
    Then Data wizard step should be selected
    And value of step readout should have text "data"
    When user selects "Protease" in Targets input in wizard
    And user closes Targets input in wizard
    And user clicks on NEXT button in wizard
    Then wizard should contain text "Creating \"Atlas\" for Kinase, Protease."
    When user clicks on BACK button in wizard
    Then Data wizard step should be selected
    And Targets input in wizard should contain text "Protease"
    When user clicks on NEXT button in wizard
    And user clicks on FINISH button in wizard
    Then value of wizard readout should have text "finished"
