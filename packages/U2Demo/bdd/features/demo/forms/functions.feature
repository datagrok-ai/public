@demo @realizes:u2.functions-browser @realizes:u2.func-form
Feature: Functions
  The functions browser is a kind with parts (search, list); the u2 funcForm is a "function form".

  Scenario: Picking a platform function and running it
    Given user opens the "Functions" demo page
    When user types "Abs" into search of functions browser
    And user clicks on Abs item in list of functions browser
    Then Abs heading should be visible
    And function form should be visible
    When user enters "7.5" into x input in function form
    And user clicks on Run button
    Then value of result readout should have text "7.5"
